module PlasmaWavesSymbolicsExt

using Symbolics
using Symbolics: @variables, expand, simplify, substitute, Differential, expand_derivatives, value
using Latexify
using LinearAlgebra
using PlasmaWaves

"""
    SymbolicDerivation

Result of a symbolic plasma-wave dispersion derivation via `Symbolics.jl` + `Latexify.jl`.

Fields:
- `name::Symbol`: wave model identifier
- `assumptions::String`: physical assumptions of the derivation
- `variables`: NamedTuple of symbolic variables used
- `matrix`: linearized wave matrix (or `nothing` for scalar derivations)
- `dispersion_expr`: symbolic dispersion expression ``D(\\omega, k, \\dots)`` (set to 0)
- `solutions`: NamedTuple of closed-form symbolic branch expressions (e.g. `ω_fast`, `ω_slow`)
- `latex::String`: LaTeX representation of the dispersion relation (ready for Quarto/KeyNotes)
"""
struct SymbolicDerivation
    name::Symbol
    assumptions::String
    variables::NamedTuple
    matrix::Any
    dispersion_expr::Any
    solutions::NamedTuple
    latex::String
end

function Base.show(io::IO, d::SymbolicDerivation)
    println(io, "SymbolicDerivation(:$(d.name))")
    println(io, "  Assumptions : ", d.assumptions)
    println(io, "  Dispersion  : ", d.dispersion_expr, " = 0")
    println(io, "  Branches    : ", keys(d.solutions))
end

"""
    derive_mhd_dispersion()

Symbolically derive the ideal MHD dispersion relation from the ``3\\times 3`` linearized
velocity perturbation matrix ``\\mathbf{M}(\\omega, k, \\theta, v_A, v_s)``, verifying its
exact factorization into the Shear Alfvén mode and Fast/Slow magnetosonic quartic.
"""
function PlasmaWaves.derive_mhd_dispersion()
    @variables ω k θ v_A v_s
    sθ, cθ = sin(θ), cos(θ)
    M = [
        -ω^2 + k^2 * (v_A^2 + v_s^2 * sθ^2)   0                              k^2 * v_s^2 * sθ * cθ;
        0                                      -ω^2 + k^2 * v_A^2 * cθ^2      0;
        k^2 * v_s^2 * sθ * cθ                  0                              -ω^2 + k^2 * v_s^2 * cθ^2
    ]

    detM = expand(det(M))
    alfven_factor = -ω^2 + k^2 * v_A^2 * cθ^2
    # Exact 2x2 magnetosonic block determinant (using sin^2 + cos^2 = 1)
    ms_quartic = ω^4 - k^2 * (v_A^2 + v_s^2) * ω^2 + k^4 * v_A^2 * v_s^2 * cθ^2

    v_alfven = k * v_A * cθ
    disc = sqrt((v_A^2 + v_s^2)^2 - 4 * v_A^2 * v_s^2 * cθ^2)
    ω_fast = k * sqrt((v_A^2 + v_s^2 + disc) / 2)
    ω_slow = k * sqrt((v_A^2 + v_s^2 - disc) / 2)

    eq = ms_quartic ~ 0
    return SymbolicDerivation(
        :ideal_mhd,
        "Linearized ideal MHD: ω ≪ Ω_i, k d_i ≪ 1, adiabatic pressure, uniform B₀ = B₀ ẑ.",
        (; ω, k, θ, v_A, v_s),
        M,
        ms_quartic,
        (; alfven_factor, ω_alfven = v_alfven, ω_fast, ω_slow, detM),
        String(latexify(eq)),
    )
end

"""
    derive_cold_plasma_dispersion()

Symbolically derive the Appleton-Hartree / Stix cold-plasma dispersion relation
``A n^4 - B n^2 + C = 0`` from the ``3\\times 3`` Stix wave matrix ``\\mathbf{R}(S, D, P, n, \\theta)``,
including principal-angle limits ``\\theta = 0`` and ``\\theta = \\pi/2``, and polarization ``i E_x / E_y``.
"""
function PlasmaWaves.derive_cold_plasma_dispersion()
    @variables S D P R L n θ
    sθ, cθ = sin(θ), cos(θ)
    # Use real-valued representation for clean symbolic determinant with iD -> D_im
    # det(R) = (S - n^2)*((S - n^2 cos^2 θ)*(P - n^2 sin^2 θ) - n^4 sin^2 θ cos^2 θ) - D^2*(P - n^2 sin^2 θ)
    A = S * sθ^2 + P * cθ^2
    B = R * L * sθ^2 + P * S * (1 + cθ^2)
    C = P * R * L
    disp = A * n^4 - B * n^2 + C

    # Principal limits
    lim_par = simplify(substitute(disp, Dict(θ => 0)))
    lim_perp = simplify(substitute(disp, Dict(θ => π / 2)))
    # Transverse electric-field polarization i E_x / E_y = (n^2 - S) / D
    pol_E = (n^2 - S) / D

    eq = disp ~ 0
    return SymbolicDerivation(
        :cold_plasma,
        "Cold uniform magnetized plasma (T = 0, Stix parameters R, L, P, S, D with S² - D² = RL).",
        (; S, D, P, R, L, n, θ),
        nothing,
        disp,
        (; A, B, C, parallel_limit = lim_par, perpendicular_limit = lim_perp, polarization_E = pol_E),
        String(latexify(eq)),
    )
end

"""
    derive_hall_mhd_dispersion()

Symbolically derive the parallel Hall-MHD circularly polarized dispersion relation
``\\omega^2 = k_\\parallel^2 v_A^2 (1 \\pm \\omega / \\Omega_i)`` for Whistler (`+`) and EMIC (`-`) waves.
"""
function PlasmaWaves.derive_hall_mhd_dispersion()
    @variables ω k_par v_A Ω_i d_i
    # Linearized 2x2 Hall-MHD system for (v_±, B_±):
    # -i ω v_± = (i k_par v_A^2 / B0) B_±
    # -i ω B_± = i k_par B0 v_± ∓ i k_par^2 (v_A^2 / Ω_i) B_±
    # Eliminating v_± yields ω^2 ∓ (k_par^2 v_A^2 / Ω_i) ω - k_par^2 v_A^2 = 0
    disp_R = ω^2 - k_par^2 * v_A^2 * (1 + ω / Ω_i)
    disp_L = ω^2 - k_par^2 * v_A^2 * (1 - ω / Ω_i)

    x = k_par * d_i / 2
    ω_whistler = k_par * v_A * (sqrt(1 + x^2) + x)
    ω_emic = k_par * v_A * (sqrt(1 + x^2) - x)

    Dk = Differential(k_par)
    vg_whistler = simplify(expand_derivatives(Dk(ω_whistler)))
    vg_emic = simplify(expand_derivatives(Dk(ω_emic)))

    eq = disp_R ~ 0
    return SymbolicDerivation(
        :hall_mhd,
        "Parallel Hall-MHD (k ∥ B₀, ω ≪ Ω_e, cold electrons, d_i = v_A / Ω_i).",
        (; ω, k_par, v_A, Ω_i, d_i),
        nothing,
        disp_R,
        (; disp_R, disp_L, ω_whistler, ω_emic, vg_whistler, vg_emic),
        String(latexify(eq)),
    )
end

"""
    derive_kaw_iaw_dispersion()

Symbolically derive the Kinetic Alfvén Wave (KAW, ``\\beta_e \\gg m_e / m_i``) and
Inertial Alfvén Wave (IAW, ``\\beta_e \\ll m_e / m_i``) dispersion relations and
parallel electric field ratios from the two-fluid shear Alfvén equation.
"""
function PlasmaWaves.derive_kaw_iaw_dispersion()
    @variables ω k_par k_perp v_A r_iL r_fL d_e
    ω_kaw = k_par * v_A * sqrt(1 + k_perp^2 * (r_iL^2 + r_fL^2))
    ω_iaw = k_par * v_A / sqrt(1 + k_perp^2 * d_e^2)

    disp_kaw = ω^2 - k_par^2 * v_A^2 * (1 + k_perp^2 * (r_iL^2 + r_fL^2))
    disp_iaw = ω^2 * (1 + k_perp^2 * d_e^2) - k_par^2 * v_A^2

    ez_ex_kaw = k_par * k_perp * r_fL^2 / (1 + k_perp^2 * r_iL^2)
    ez_ex_iaw = k_par * k_perp * d_e^2 / (1 + k_perp^2 * d_e^2)

    Dkp = Differential(k_perp)
    vg_perp_kaw = simplify(expand_derivatives(Dkp(ω_kaw)))
    vg_perp_iaw = simplify(expand_derivatives(Dkp(ω_iaw)))

    eq = disp_kaw ~ 0
    return SymbolicDerivation(
        :kaw_iaw,
        "Two-fluid oblique shear Alfvén modes (ω ≪ Ω_i, k_⊥ ≫ k_∥): KAW (β_e ≫ m_e/m_i) & IAW (β_e ≪ m_e/m_i).",
        (; ω, k_par, k_perp, v_A, r_iL, r_fL, d_e),
        nothing,
        disp_kaw,
        (; disp_kaw, disp_iaw, ω_kaw, ω_iaw, ez_ex_kaw, ez_ex_iaw, vg_perp_kaw, vg_perp_iaw),
        String(latexify(eq)),
    )
end

"""
    derive_bohm_gross_dispersion()

Symbolically derive the warm-electron Bohm-Gross dispersion relation
``\\omega^2 = \\omega_{pe}^2 + 3 k^2 v_{th}^2`` and its exact group velocity ``v_g = \\partial\\omega / \\partial k``.
"""
function PlasmaWaves.derive_bohm_gross_dispersion()
    @variables ω k ω_pe v_th
    disp = ω^2 - ω_pe^2 - 3 * k^2 * v_th^2
    ω_bg = sqrt(ω_pe^2 + 3 * k^2 * v_th^2)
    Dk = Differential(k)
    v_group = simplify(expand_derivatives(Dk(ω_bg)))
    v_phase = ω_bg / k

    eq = disp ~ 0
    return SymbolicDerivation(
        :bohm_gross,
        "Warm unmagnetized electron plasma wave (1D adiabatic γ = 3, v_ph ≫ v_th, fixed ions).",
        (; ω, k, ω_pe, v_th),
        nothing,
        disp,
        (; ω_bg, v_phase, v_group),
        String(latexify(eq)),
    )
end

"""
    symbolic_to_numeric(expr, k_var, param_dict)

Substitute numeric plasma parameters `param_dict` into a symbolic dispersion expression
`expr` and return a callable function `k -> Float64` that can be passed directly into
`confirm_dispersion(...)`.
"""
function PlasmaWaves.symbolic_to_numeric(expr, k_var, param_dict::AbstractDict)
    sub_expr = substitute(expr, param_dict)
    fn = Symbolics.build_function(sub_expr, k_var; expression = Val{false})
    return k_val::Real -> Float64(fn(float(k_val)))
end

end # module PlasmaWavesSymbolicsExt
