# =====================================================================================
# UPWELLING SUBGLACIAL-DISCHARGE-PLUME LES — GEOMETRIC terminus, Oceananigans 0.109
#
# Newer-Oceananigans build of ../newLES/iceplume.jl. Instead of tilting gravity, the terminus
# is represented GEOMETRICALLY with an immersed ice face and VERTICAL gravity — the faithful
# setup (matches the Ovall et al. 2025 OLEM approach and keeps overcut plumes ice-attached,
# unlike the gravity-tilt idealisation). Built on the 0.109.2 API of the sibling IceShelfCavity
# project: ImmersedBoundaryGrid + GridFittedBoundary + ConjugateGradientPoissonSolver (clean
# immersed-boundary pressure — no FFT near-wall divergence).
#
# Coordinates (all TRUE, gravity = (0,0,-1)):
#   x = distance into the fjord from the glacier side (x=0 west, x=Lx east/ocean)
#   y = along-glacier (centred at 0)
#   z = height above the grounding line: z=0 flat sea floor, z=Lz flat surface (both are just
#       the domain boundaries — no cavity mask needed).
#
# Ice terminus = solid immersed region x < x_face(z), x_face(z) = xf_a + xf_b·z:
#   VERTICAL  (--face_angle=90)            : no immersed ice; ice = west wall at x=0.
#   OVERCUT   (--terminus=overcut, default): x_face(z) = (Lz−z)·tanβ  (base sticks out to
#             x=Lz·tanβ, recedes to x=0 at the surface — slopes AWAY from the ocean going up).
#   UNDERCUT  (--terminus=undercut)        : x_face(z) = z·tanβ  (top overhangs toward ocean).
#   β = 90 − face_angle.  Discharge enters at the grounding-line corner (x_face(0), 0).
#
# >>> FIRST-RUN CHECKS (CPU) <<<
#   1. Ice geometry: quick-look shows the immersed ice blanked and leaning the intended way.
#   2. Plume attaches to the ice face and rises; discharge is fresh at the base.
#   3. CFL < 0.5, no NaNs; CG solver converges (watch the progress log).
# =====================================================================================
using Oceananigans
using Oceananigans.Units
using Oceananigans.Solvers: ConjugateGradientPoissonSolver
using Oceananigans.BoundaryConditions: OpenBoundaryCondition
using NCDatasets            # required so NetCDFWriter (an extension) is available
using Printf
using Statistics: mean
using Oceanostics
using CUDA: has_cuda_gpu

#+++ Command-line arguments (minimal --key=value parser)
function parse_cli()
    cli = Dict{String,Any}(
        "simname" => "iceplume_cg", "arch" => "auto", "terminus" => "overcut",
        "face_angle" => 90.0, "discharge" => 150.0, "outlet_w" => 24.0, "outlet_h" => 6.0,
        "Lz" => 150.0, "Ly" => 192.0, "Lx" => 500.0, "dz" => 0.75, "fine_x" => 375.0,
        "dx_max" => 18.6, "stop_time" => 45.0, "output_interval" => 300.0,   # dz/fine_x/dx_max = Ovall 2025
        "fine_y" => 100.0, "dy_max" => 8.0,   # y: uniform dz within |y|<fine_y (middle 2·fine_y), then stretch to dy_max
        "fine_z" => 120.0, "dz_surf" => 4.0,  # z: uniform dz to fine_z, then stretch to dz_surf toward the surface
        "checkpoint_interval" => 5.0, "wall_time_limit" => Inf,
        "cg_reltol" => 1e-5, "cg_maxiter" => 30.0,   # CG Poisson solver: looser reltol = fewer iters
        "channel_len" => 10.0,   # minimum subglacial-channel length in x [m]; face shifts out if shorter
        "sig_src" => 2.0,        # discharge-source relaxation time [s]; larger = gentler (tune vs Δt)
        "outdir" => "")   # output + checkpoints dir; empty ⇒ <rundir>/output
    provided = Set{String}()
    for a in ARGS
        startswith(a, "--") || continue
        kv = split(a[3:end], "="; limit = 2)
        length(kv) == 2 || error("bad argument '$a' (use --key=value)")
        k, v = kv[1], kv[2]
        haskey(cli, k) || error("unknown argument --$k")
        cli[k] = cli[k] isa AbstractString ? String(v) : parse(Float64, v)
        push!(provided, k)
    end
    return cli, provided
end
cli, provided = parse_cli()
rundir = @__DIR__
include(joinpath(rundir, "ambient_profile.jl"))   # T_ambient(z), S_ambient(z): 2024 cast
#---

#+++ Architecture
arch = cli["arch"] == "cpu" ? CPU() : cli["arch"] == "gpu" ? GPU() : (has_cuda_gpu() ? GPU() : CPU())
@info "Simulation $(cli["simname"]) on $arch  ($(cli["terminus"]), face_angle=$(cli["face_angle"])°, Q=$(cli["discharge"]) m³/s)"
#---

#+++ Terminus geometry: ice = immersed solid  x < x_face(z) = xf_a + xf_b·z.
# ALWAYS immersed (even vertical), so every case uses the GPU-safe immersed-drag path. A
# continuous drag flux BC on the west *domain* boundary does not compile on the GPU in
# Oceananigans 0.109 ("dynamic call to τv"); an immersed drag BC does.
β = 90.0 - cli["face_angle"]              # tilt magnitude from vertical [deg]
tanβ = tand(β)
x_wall = 3.0                               # thickness of the vertical ice slab [m] (a few cells)
if β ≤ 0                                    # vertical: thin immersed wall (axis-aligned ⇒ no staircase)
    xf_a0, xf_b = x_wall, 0.0
elseif cli["terminus"] == "undercut"        # top overhangs toward the ocean
    xf_a0, xf_b = 0.0, tanβ
else                                        # overcut (default): base sticks out, recedes up
    xf_a0, xf_b = cli["Lz"] * tanβ, -tanβ
end
# Shift the ice face out (in +x) so the subglacial channel is at least channel_len long at its
# narrowest point over the channel height. This matters for the ~vertical face (slab only a few m
# thick); the overcut base already extends tens of m, so the shift is a no-op there.
min_face0 = xf_a0 + min(xf_b * cli["outlet_h"], 0.0)     # smallest x_face(z) for z ∈ [0, H]
x_shift   = max(0.0, cli["channel_len"] - min_face0)
xf_a      = xf_a0 + x_shift
immersed_ice = true                         # always immerse the ice
x_gl = xf_a                                 # grounding-line x-position (base of ice)
#---

#+++ Grid — fine near the terminus/plume, stretched away in x, y AND z (vertical gravity ⇒ z = height)
# one-sided stretch 0→L: uniform dx_fine out to `fine`, then geometric growth capped at d_max.
function build_faces(dx_fine, fine, d_max, L; growth = 1.03)
    f = Float64[0.0]
    while f[end] < min(fine, L); push!(f, f[end] + dx_fine); end
    d = dx_fine
    while f[end] < L; d = min(d * growth, d_max); push!(f, f[end] + d); end
    f[end] = L
    # Clipping to L can leave a sliver last cell; a tiny cell next to fast flow (e.g. the surface
    # plume-impingement, or the open outflow) is a stability hazard — merge it into its neighbor.
    if length(f) > 2 && (f[end] - f[end-1]) < 0.5 * dx_fine
        deleteat!(f, length(f) - 1)
    end
    return f
end
# symmetric stretch on [-L/2,+L/2]: uniform within |y|<half_fine, stretch to d_max toward both edges.
function build_sym_faces(dx_fine, half_fine, d_max, L; growth = 1.03)
    pos = build_faces(dx_fine, half_fine, d_max, L/2; growth = growth)   # 0 .. L/2
    return vcat(-reverse(pos[2:end]), pos)                               # -L/2 .. 0 .. +L/2
end

Lz = cli["Lz"]; Ly = cli["Ly"]; Lx = cli["Lx"]
dx_fine = cli["dz"]; fine_x = cli["fine_x"]; dx_max = cli["dx_max"]
fine_y = cli["fine_y"]; dy_max = cli["dy_max"]; fine_z = cli["fine_z"]; dz_surf = cli["dz_surf"]
if arch == CPU() && !("dz" in provided)     # laptop smoke test: coarse but ice-resolving
    @warn "CPU with no --dz: reduced smoke-test resolution (2 m; not for science)."
    Lx = min(Lx, 220.0); Ly = min(Ly, 120.0); dx_fine = 2.0; fine_x = 80.0; dx_max = 12.0
    fine_y = 40.0; dy_max = 8.0; fine_z = 90.0; dz_surf = 6.0
elseif arch == CPU()
    @warn "CPU at user resolution dz=$dx_fine on Lx=$Lx, Ly=$Ly — intermediate/overnight run."
end
# x: fine to fine_x then stretch to dx_max (as Ovall 2025).  y: fine middle 2·fine_y, stretch to dy_max.
# z: fine to fine_z, stretch to dz_surf near the surface (relaxes the surface-impingement CFL).
x_faces = build_faces(dx_fine, fine_x, dx_max, Lx)
y_faces = build_sym_faces(dx_fine, fine_y, dy_max, Ly)
z_faces = build_faces(dx_fine, fine_z, dz_surf, Lz)
Nx = length(x_faces) - 1; Ny = length(y_faces) - 1; Nz = length(z_faces) - 1

grid_base = RectilinearGrid(arch; size = (Nx, Ny, Nz),
                            x = x_faces, y = y_faces, z = z_faces,
                            topology = (Bounded, Bounded, Bounded))

if immersed_ice
    # Ice is solid where x < x_face(z), EXCEPT the subglacial channel: a W×H (= outlet cross-
    # section) conduit extruded in x through the ice base, from the back wall out to the terminus
    # face. Carving it to ocean (0) turns the nudged source region (|y|<W/2, z<H) into a real
    # channel, so the discharge exits as a coherent, developed jet through an ice-face opening —
    # instead of materializing in a 2 m sliver at the face and spraying off (the old behavior).
    ice_mask = let a = xf_a, b = xf_b, Wc = cli["outlet_w"], Hc = cli["outlet_h"]
        (x, y, z) -> ifelse((x < a + b * z) & !((abs(y) < Wc/2) & (z < Hc)), 1, 0)  # 1 ice, 0 ocean
    end
    grid = ImmersedBoundaryGrid(grid_base, GridFittedBoundary(ice_mask))
    @info "Immersed $(cli["terminus"]) ice face: x_face(z)=$(round(xf_a))$(xf_b<0 ? "" : "+")$(round(xf_b,digits=3))·z; grounding line x=$(round(x_gl)) m"
else
    grid = grid_base
    @info "Vertical terminus: ice = west wall (no immersed ice)."
end
@info "Grid" grid Nx Ny Nz total_cells=Nx*Ny*Nz
#---

#+++ Parameters
U_in = cli["discharge"] / (cli["outlet_w"] * cli["outlet_h"])   # channel inflow speed [m/s]
# SMOOTH east-outflow profile g(z)=½(1+tanh((z−z_out)/δ_out)): outflow concentrated above 60 m depth
# but WITHOUT the hard jump (the discontinuity stalls the CG preconditioner). Peak amplitude is
# normalized so the discrete flux ∫∫ u_out dy dz = Q exactly (mass balance ⇒ CG stays well-posed).
z_out = Lz - 60.0; δ_out = 15.0
zC = 0.5 .* (z_faces[1:end-1] .+ z_faces[2:end]); Δzc = diff(z_faces)
I_s = sum(0.5 .* (1 .+ tanh.((zC .- z_out) ./ δ_out)) .* Δzc)      # ∫ g(z) dz  [m]
U_out_peak = cli["discharge"] / (Ly * I_s)                        # Ly · ∫g · U_peak = Q
params = (; Lz, Lx, Ly, xf_a, xf_b,
          W = cli["outlet_w"], H = cli["outlet_h"], U_in,
          U_out_peak, z_out, δ_out,
          σ_src = cli["sig_src"], σ_spg = 10.0, x_src = 3 * dx_fine,
          L_sponge = 60.0,
          Cᴰ = 2.5e-3)
@info "Derived" U_in U_out_peak I_s x_gl
#---

#+++ Ambient (vertical gravity ⇒ z is true height; no rotation)
@inline T∞(z) = T_ambient(z)
@inline S∞(z) = S_ambient(z)
#---

#+++ Discharge freshwater source + fjord-side tracer sponge. The u VOLUME (and the one-way
#    subglacial-channel RUNWAY) is set by the OPEN boundaries below — a momentum nudge against the
#    closed back wall cannot drive net through-flow, so open boundaries are required for the runway.
@inline x_face(z, p) = p.xf_a + p.xf_b * z
# Keep the carved channel fresh/cold (T=0,S=0) so the west inflow is buoyant. The along-channel jet
# comes from the west-inflow BC, not an interior u-nudge.
@inline in_outlet(x, y, z, p) = (x < x_face(z, p) + p.x_src) & (abs(y) < p.W/2) & (z < p.H)
@inline src_T(x,y,z,t,T,p) = ifelse(in_outlet(x,y,z,p), -(T - 0.0)/p.σ_src, zero(T))
@inline src_S(x,y,z,t,S,p) = ifelse(in_outlet(x,y,z,p), -(S - 0.0)/p.σ_src, zero(S))

# Fjord (east) tracer sponge: relax T,S→ambient near the outflow so exiting water is reasonable.
@inline east_frac(x, p) = clamp((x - (p.Lx - p.L_sponge)) / p.L_sponge, 0.0, 1.0)
@inline spg_T(x,y,z,t,T,p) = -east_frac(x,p)/p.σ_spg * (T - T∞(z))
@inline spg_S(x,y,z,t,S,p) = -east_frac(x,p)/p.σ_spg * (S - S∞(z))

@inline f_T(x,y,z,t,T,p) = src_T(x,y,z,t,T,p) + spg_T(x,y,z,t,T,p)
@inline f_S(x,y,z,t,S,p) = src_S(x,y,z,t,S,p) + spg_S(x,y,z,t,S,p)
FT = Forcing(f_T, field_dependencies = :T, parameters = params)
FS = Forcing(f_S, field_dependencies = :S, parameters = params)
forcing = (T = FT, S = FS)
#---

# Quadratic drag on the ICE, always via the IMMERSED boundary (GPU-safe for all cases; a
# continuous drag BC on the west domain boundary won't compile on GPU in 0.109).
@inline _spd(u,v,w) = √(u^2 + v^2 + w^2)
@inline τu(x,y,z,t,u,v,w,p) = -p.Cᴰ * u * _spd(u,v,w)
@inline τv(x,y,z,t,u,v,w,p) = -p.Cᴰ * v * _spd(u,v,w)
@inline τw(x,y,z,t,u,v,w,p) = -p.Cᴰ * w * _spd(u,v,w)
bcpar = (; Cᴰ = params.Cᴰ)
τu_bc = FluxBoundaryCondition(τu, field_dependencies=(:u,:v,:w), parameters=bcpar)
τv_bc = FluxBoundaryCondition(τv, field_dependencies=(:u,:v,:w), parameters=bcpar)
τw_bc = FluxBoundaryCondition(τw, field_dependencies=(:u,:v,:w), parameters=bcpar)

# Discharge VOLUME + runway via prescribed OPEN boundaries: west inflow over the channel opening
# (u=U_in ⇒ Q = U_in·W·H drives the along-channel jet — the runway), balanced by a SMOOTH east
# outflow g(z)=½(1+tanh((z−z_out)/δ_out)) scaled to U_out_peak so ∫∫ = Q. The smooth profile avoids
# the hard 60 m-depth jump that stalls the CG preconditioner. Pure arithmetic ⇒ GPU-safe.
@inline u_west_in(y,z,t,p)  = ifelse((abs(y) < p.W/2) & (z < p.H), p.U_in, zero(p.U_in))
@inline u_east_out(y,z,t,p) = p.U_out_peak * 0.5*(1 + tanh((z - p.z_out)/p.δ_out))
u_bcs = FieldBoundaryConditions(immersed = τu_bc,
                                west = OpenBoundaryCondition(u_west_in,  parameters = params),
                                east = OpenBoundaryCondition(u_east_out, parameters = params))
v_bcs = FieldBoundaryConditions(immersed = τv_bc)
w_bcs = FieldBoundaryConditions(immersed = τw_bc)
T_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(0))
S_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(0))
boundary_conditions = (u=u_bcs, v=v_bcs, w=w_bcs, T=T_bcs, S=S_bcs)
#---

#+++ Model (vertical gravity; CG Poisson solver for the immersed boundary)
eos = LinearEquationOfState(thermal_expansion = 3.87e-5, haline_contraction = 7.86e-4)
# QuasiAdamsBashforth2: ONE pressure solve per step (vs THREE for the RK3 default) — ~3× fewer
# CG solves. The CG solver already uses the FFT solver as a preconditioner, so a modest maxiter
# and a looser reltol keep the iteration count low; tune via --cg_reltol / --cg_maxiter.
model = NonhydrostaticModel(grid;
                            timestepper = :QuasiAdamsBashforth2,
                            buoyancy = SeawaterBuoyancy(equation_of_state = eos),
                            coriolis = FPlane(f = 1.22e-4),
                            advection = WENO(order = 5),
                            tracers = (:T, :S),
                            closure = AnisotropicMinimumDissipation(),
                            forcing = forcing,
                            pressure_solver = ConjugateGradientPoissonSolver(grid;
                                reltol = cli["cg_reltol"], maxiter = round(Int, cli["cg_maxiter"])),
                            boundary_conditions = boundary_conditions)
@info "Model built" model
#---

#+++ Initial condition: 2024 ambient + small noise
u, v, w = model.velocities
Tᵢ(x,y,z) = T∞(z); Sᵢ(x,y,z) = S∞(z)
uᵢ = 5e-3 .* (rand(size(u)...) .- 0.5)
vᵢ = 5e-3 .* (rand(size(v)...) .- 0.5)
wᵢ = 5e-3 .* (rand(size(w)...) .- 0.5)
uᵢ .-= mean(uᵢ); vᵢ .-= mean(vᵢ); wᵢ .-= mean(wᵢ)
set!(model, u = uᵢ, v = vᵢ, w = wᵢ, T = Tᵢ, S = Sᵢ)
#---

#+++ Simulation
Δt₀ = 0.2 * minimum_zspacing(grid_base) / max(U_in, 0.1)
wtl = cli["wall_time_limit"]
simulation = Simulation(model; Δt = Δt₀, stop_time = cli["stop_time"] * minutes,
                        wall_time_limit = isfinite(wtl) ? wtl * 3600 : Inf)
simulation.callbacks[:wizard] = Callback(TimeStepWizard(cfl=0.35, max_change=1.05, min_change=0.2, max_Δt=1.0),
                                         IterationInterval(5))
using Oceanostics.ProgressMessengers: BasicTimeMessenger
simulation.callbacks[:progress] = Callback(BasicTimeMessenger(), IterationInterval(50))

# Watch what sets the CFL: print max|u|,|v|,|w| (the surface-impingement velocities dominate Δt).
using Oceananigans.Fields: interior
@inline _vmax(f) = maximum(abs, interior(f))
function log_vmax(sim)
    u, v, w = sim.model.velocities
    @info @sprintf("      max|u|=%.2f  max|v|=%.2f  max|w|=%.2f m/s   Δt=%.3f s",
                   _vmax(u), _vmax(v), _vmax(w), sim.Δt)
    return nothing
end
simulation.callbacks[:vmax] = Callback(log_vmax, IterationInterval(50))
#---

#+++ Outputs
T, S = model.tracers
ω_y = Field(∂z(u) - ∂x(w))
outputs = (; u, v, w, T, S, ω_y)
prefix = cli["simname"]
ckpt = "checkpoint_" * prefix
# Output + checkpoints go to --outdir (default <rundir>/output). Keep this OFF the git repo
# and on /glade/work or scratch for cluster runs.
outdir = isempty(cli["outdir"]) ? joinpath(rundir, "output") : cli["outdir"]
mkpath(outdir)
pickup = any(startswith("$(ckpt)_iteration"), readdir(outdir))
overwrite = !pickup
pickup && @warn "Checkpoint found for $prefix in $outdir — resuming."

# Metadata so the (non-rotating) quick-look can blank the immersed ice and label geometry.
# y-z "face" slice: a few cells into the fluid off the grounding line (a fixed near-wall x=1
# slice would sit inside the ice wedge for overcut, so slice through the plume near the base).
face_ix = clamp(round(Int, x_gl / dx_fine) + 3, 1, Nx)
face_x = x_faces[face_ix]
gattrs = Dict("terminus" => cli["terminus"], "face_angle_deg" => cli["face_angle"],
              "geometry" => "immersed_vertical_gravity", "theta_tilt_deg" => 0.0,
              "xf_a" => xf_a, "xf_b" => xf_b, "water_depth_m" => Lz,
              "face_x_m" => face_x, "discharge_m3s" => cli["discharge"],
              "channel_w_m" => cli["outlet_w"], "channel_h_m" => cli["outlet_h"],
              "channel_len_m" => cli["channel_len"])

simulation.output_writers[:fields] = NetCDFWriter(model, outputs;
    filename = joinpath(outdir, "$(prefix)_fields.nc"),
    schedule = TimeInterval(cli["output_interval"] * seconds),
    global_attributes = gattrs, overwrite_existing = overwrite)
simulation.output_writers[:face] = NetCDFWriter(model, outputs;
    filename = joinpath(outdir, "$(prefix)_face.nc"),
    schedule = TimeInterval(10seconds), indices = (face_ix, :, :),
    global_attributes = gattrs, overwrite_existing = overwrite)
simulation.output_writers[:midy] = NetCDFWriter(model, outputs;
    filename = joinpath(outdir, "$(prefix)_midy.nc"),
    schedule = TimeInterval(10seconds), indices = (:, round(Int, Ny/2), :),
    global_attributes = gattrs, overwrite_existing = overwrite)

uw = Field(@at (Center, Center, Center) u * w)
wT = Field(@at (Center, Center, Center) w * T)
wS = Field(@at (Center, Center, Center) w * S)
# 15-min averaging window (Ovall 2025 average the last 15 min of a 45-min run; the final
# record here covers minutes 30–45).
simulation.output_writers[:avg] = NetCDFWriter(model, (; u,v,w,T,S,uw,wT,wS);
    filename = joinpath(outdir, "$(prefix)_timeavg.nc"),
    schedule = AveragedTimeInterval(15minutes, window = 15minutes),
    global_attributes = gattrs, overwrite_existing = overwrite)

simulation.output_writers[:checkpointer] = Checkpointer(model;
    schedule = TimeInterval(cli["checkpoint_interval"] * minutes),
    dir = outdir, prefix = ckpt, cleanup = true)
#---

@info "Starting run..." pickup stop_minutes=cli["stop_time"]
run!(simulation; pickup)
@info "Done: $(prefix)"
