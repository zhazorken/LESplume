# =====================================================================================
#  2024 ambient T,S profile for the LeConte / Xeitl Sít' subglacial-discharge-plume LES.
#
#  Source: AmbientProfile_15JUL2024_4modeling.mat  (fields: CT, SA, depth), the 15-Jul-2024
#  cast used for the Ovall et al. (2025, JGR-Oceans) study, MODIFIED so that properties are
#  held CONSTANT below the lowest point on the crest of the moraine (~125 m depth), because
#  deeper fjord water cannot cross the sill to reach the glacier.
#
#  Model vertical coordinate:  z = height above the grounding line, 0 (grounding line, the
#  deepest point) .. Lz (free surface).  With the grounding line at 150 m depth,
#        depth = Lz - z,   Lz = 150 m,   moraine crest at z = z_crest = 25 m (depth 125 m).
#
#  The measured cast is represented by a GPU-safe closed form: a degree-8 polynomial in the
#  normalised height  ζ = clamp((z - z_crest)/(Lz - z_crest), 0, 1)  fitted (endpoint-
#  weighted) to the data over the variable region z ≥ z_crest, and CLAMPED to the crest
#  value for z < z_crest (ζ = 0) so the "constant-below-moraine" condition is exact.
#
#  Fit quality vs. the 800-point cast:  T RMSE 0.031 °C (max 0.11),  S RMSE 0.043 g/kg
#  (max 0.12); salinity is monotonic. See ambient_2024.csv (raw + fitted) and
#  ambient_fit.png. Coefficients were generated from the .mat by scripts/fit_ambient.py.
#
#  NOTE ON VARIABLES: CT (conservative temperature) and SA (absolute salinity) from the .mat
#  are used directly as the model's T and S with the LinearEquationOfState below. For this
#  temperate-fjord water the CT↔potential-temperature and SA↔practical-salinity offsets are
#  <~0.02 °C and ~0.16 g/kg — negligible for an idealised linear-EOS plume study. Swap in a
#  TEOS-10 conversion here if you need absolute density accuracy.
# =====================================================================================

const Z_CREST = 25.0      # m   height above grounding line of the moraine-crest (depth 125 m)
const LZ_AMB  = 150.0     # m   grounding-line depth used when the profile was built

# Degree-8 polynomial coefficients in ζ, highest power first (Horner form below).
const CT_COEFFS = (831.95325, -2965.7683, 4129.1725, -2849.3947, 1013.7086,
                   -157.58479, -11.772606, 7.4299025, 6.9838877)
const CS_COEFFS = (1615.4451, -6402.6066, 10206.931, -8425.0683, 3889.064,
                   -1028.8631, 155.89357, -15.977659, 29.956475)

@inline function _horner(coeffs, ζ)
    v = coeffs[1]
    @inbounds for i in 2:length(coeffs)
        v = v * ζ + coeffs[i]
    end
    return v
end

# ζ ∈ [0,1]; 0 at/below the moraine crest (⇒ constant deep values), 1 at the surface.
@inline _zeta(z) = clamp((z - Z_CREST) / (LZ_AMB - Z_CREST), 0.0, 1.0)

"Ambient (conservative) temperature [°C] as a function of height above the grounding line."
@inline T_ambient(z) = _horner(CT_COEFFS, _zeta(z))

"Ambient (absolute) salinity [g/kg] as a function of height above the grounding line."
@inline S_ambient(z) = _horner(CS_COEFFS, _zeta(z))
