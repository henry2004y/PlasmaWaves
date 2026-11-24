"""
Cold uniform plasma dispersion relation (Appleton-Hartree / Stix),
Hall-MHD dispersive waves (Whistler & Ion Cyclotron / EMIC),
O/X modes, Upper/Lower Hybrid resonances, Faraday rotation, and CMA diagram boundaries.
"""

"""
    stix_parameters(plasma::PlasmaParams, ω::Number)

Compute Stix cold-plasma parameters `(; R, L, P, S, D)` at frequency `ω`.
Uses signed cyclotron frequencies ``\\Omega_s = q_s B_0 / m_s`` so that
``R = 1 - \\sum_s \\frac{\\omega_{ps}^2}{\\omega(\\omega + \\Omega_s)}`` (resonant at ``|\\Omega_e|``)
and ``L = 1 - \\sum_s \\frac{\\omega_{ps}^2}{\\omega(\\omega - \\Omega_s)}`` (resonant at ``|\\Omega_i|``).
"""
function stix_parameters(plasma::PlasmaParams, ω::Number)
    ωc = complex(float(ω))
    R = one(ωc)
    L = one(ωc)
    P = one(ωc)
    for s in plasma.species
        ωps2 = plasma_frequency(s, plasma)^2
        Ωs = cyclotron_frequency(s, plasma)
        R -= ωps2 / (ωc * (ωc + Ωs))
        L -= ωps2 / (ωc * (ωc - Ωs))
        P -= ωps2 / ωc^2
    end
    S = (R + L) / 2
    D = (R - L) / 2
    return (; R, L, P, S, D)
end

"Cold-plasma relative dielectric tensor ``\\boldsymbol{\\epsilon}_R(\\omega)`` in Stix form."
function cold_dielectric_tensor(plasma::PlasmaParams, ω::Number)
    (; S, D, P) = stix_parameters(plasma, ω)
    return SMatrix{3,3,typeof(S)}(
        S,     im * D, zero(S),
       -im * D, S,     zero(S),
        zero(S), zero(S), P,
    )
end

"""
    cold_cutoffs(plasma::PlasmaParams)

Compute the three cold-plasma cutoff frequencies `(; ω_pe, ω_R, ω_L)` for an
electron-ion plasma (solving `P = 0`, `R = 0`, `L = 0`).
"""
function cold_cutoffs(plasma::PlasmaParams)
    e = get_electron(plasma)
    i = get_ion(plasma)
    ωpe = plasma_frequency(e, plasma)
    ωpi = plasma_frequency(i, plasma)
    Ωe = abs(cyclotron_frequency(e, plasma))
    Ωi = abs(cyclotron_frequency(i, plasma))

    # Exact electron-ion cutoff approximations (including ion mass ratio correction)
    ω_p = sqrt(ωpe^2 + ωpi^2)
    ω_R = 0.5 * ((Ωe - Ωi) + sqrt((Ωe + Ωi)^2 + 4 * (ωpe^2 + ωpi^2)))
    ω_L = 0.5 * (-(Ωe - Ωi) + sqrt((Ωe + Ωi)^2 + 4 * (ωpe^2 + ωpi^2)))
    return (; ω_pe = ω_p, ω_R, ω_L)
end

"""
    cold_resonances(plasma::PlasmaParams)

Compute principal cold-plasma resonance frequencies `(; Ω_i, Ω_e, ω_LH, ω_UH)`.
"""
function cold_resonances(plasma::PlasmaParams)
    e = get_electron(plasma)
    i = get_ion(plasma)
    ωpe = plasma_frequency(e, plasma)
    ωpi = plasma_frequency(i, plasma)
    Ωe = abs(cyclotron_frequency(e, plasma))
    Ωi = abs(cyclotron_frequency(i, plasma))

    # Exact two-species quadratic solution of S(ω) = 0:
    # ω⁴ - (ωpe² + ωpi² + Ωe² + Ωi²) ω² + (Ωe² Ωi² + ωpe² Ωi² + ωpi² Ωe²) = 0
    B_s = ωpe^2 + ωpi^2 + Ωe^2 + Ωi^2
    C_s = Ωe^2 * Ωi^2 + ωpe^2 * Ωi^2 + ωpi^2 * Ωe^2
    disc = sqrt(max(0.0, B_s^2 - 4 * C_s))
    ω_UH = sqrt(0.5 * (B_s + disc))
    ω_LH = sqrt(2 * C_s / (B_s + disc))
    return (; Ω_i = Ωi, Ω_e = Ωe, ω_LH, ω_UH)
end

"""
    cold_plasma_n2(plasma::PlasmaParams, ω::Real, θ::Real)

Solve the cold-plasma dispersion relation ``A n^4 - B n^2 + C = 0`` at frequency `ω`
and propagation angle `θ`. Returns `(n2_plus, n2_minus)` corresponding to ``(B \\pm F)/(2A)``.
"""
function cold_plasma_n2(plasma::PlasmaParams, ω::Real, θ::Real)
    (; R, L, P, S, D) = stix_parameters(plasma, ω)
    Rr, Lr, Pr, Sr, Dr = real(R), real(L), real(P), real(S), real(D)
    s2 = sin(θ)^2
    c2 = cos(θ)^2

    A = Sr * s2 + Pr * c2
    B = Rr * Lr * s2 + Pr * Sr * (1 + c2)
    C = Pr * Rr * Lr
    F = sqrt(max(0.0, (Rr * Lr - Pr * Sr)^2 * s2^2 + 4 * Pr^2 * Dr^2 * c2))

    if abs(A) < 1e-14
        return (Inf, C / B)
    end
    return ((B + F) / (2A), (B - F) / (2A))
end

"""
    cold_plasma_eigenmode(plasma::PlasmaParams, ω::Real, θ::Real, n2::Real; mode::Symbol=:cold_plasma, amp::Real=0.15)

Given a propagating branch ``n^2 > 0``, construct the full `WaveSolution` by finding the
exact complex electric-field eigenvector from the nullspace of the Stix wave matrix.
"""
function cold_plasma_eigenmode(
    plasma::PlasmaParams,
    ω::Real,
    θ::Real,
    n2::Real;
    mode::Symbol = :cold_plasma,
    amp::Real = 0.15,
)
    n2 <= 0 && error("Non-propagating (evanescent) branch: n² = $n2 <= 0 at ω = $ω, θ = $θ.")
    n = sqrt(n2)
    k = ω * n / plasma.c
    sθ, cθ = sin(θ), cos(θ)
    kvec = SVector(Float64(k * sθ), 0.0, Float64(k * cθ))

    (; S, D, P) = stix_parameters(plasma, ω)
    M = Matrix{ComplexF64}([
        S - n2 * cθ^2      -im * D        n2 * sθ * cθ;
        im * D              S - n2        0.0;
        n2 * sθ * cθ        0.0           P - n2 * sθ^2
    ])
    F = svd(M)
    e_unit = SVector{3,ComplexF64}(F.V[:, 3])
    # Fix global phase so the dominant real component is positive
    idx_max = argmax(abs.(e_unit))
    phase_fix = abs(e_unit[idx_max]) > 0 ? conj(e_unit[idx_max]) / abs(e_unit[idx_max]) : 1.0 + 0im
    e_unit = e_unit * phase_fix

    # Scale so |δB| ≈ amp * B0 (or |δE| = amp * c * B0 if purely electrostatic)
    b_unit = cross(kvec, e_unit) / ω
    bnorm = norm(b_unit)
    if bnorm > 1e-12 * norm(e_unit) / plasma.c
        scale = amp * plasma.B0 / bnorm
    else
        scale = amp * alfven_speed(plasma) * plasma.B0
    end
    δE = scale * e_unit
    δB = scale * b_unit

    # Ion fluid velocity from cold equation of motion
    ion = get_ion(plasma)
    Ωi = cyclotron_frequency(ion, plasma)
    denom = 1.0 - (Ωi / ω)^2
    if abs(denom) > 1e-10
        fac = im * ion.q / (ion.m * ω)
        vix = fac * (δE[1] + im * (Ωi / ω) * δE[2]) / denom
        viy = fac * (δE[2] - im * (Ωi / ω) * δE[1]) / denom
        viz = fac * δE[3]
        δv = SVector(vix, viy, viz)
    else
        δv = SVector(ComplexF64(0), ComplexF64(0), ComplexF64(0))
    end
    δn = ion.n * dot(kvec, δv) / ω
    δp = ComplexF64(0)

    return WaveSolution(mode, ComplexF64(ω), kvec, δE, δB, δv, δn, δp, plasma)
end

"""
    hall_mhd_dispersion(plasma::PlasmaParams, k_par::Real; branch::Symbol=:whistler)

Exact parallel Hall-MHD dispersion relation:
``\\omega / (k_\\parallel v_A) = \\sqrt{1 + (k_\\parallel d_i / 2)^2} \\pm k_\\parallel d_i / 2``.
`branch = :whistler` (`+`, right-hand) or `:emic` / `:ion_cyclotron` (`-`, left-hand).
"""
function hall_mhd_dispersion(plasma::PlasmaParams, k_par::Real; branch::Symbol = :whistler)
    vA = alfven_speed(plasma)
    di = inertial_length(get_ion(plasma), plasma)
    x = abs(k_par) * di / 2
    if branch === :whistler || branch === :right
        return abs(k_par) * vA * (sqrt(1 + x^2) + x)
    elseif branch === :emic || branch === :ion_cyclotron || branch === :left
        return abs(k_par) * vA * (sqrt(1 + x^2) - x)
    else
        error("Unknown Hall-MHD branch: $branch. Use :whistler or :emic.")
    end
end

function _whistler_n2(plasma::PlasmaParams, ω::Real, θ::Real)
    if abs(θ) < 1e-12
        return real(stix_parameters(plasma, ω).R)
    end
    n2p, n2m = cold_plasma_n2(plasma, ω, θ)
    if n2p > 0 && n2m > 0
        # Below Ω_i both branches propagate; the smaller n² is the fast/whistler (R) mode
        return min(n2p, n2m)
    else
        return max(n2p, n2m)
    end
end

"""
    whistler_wave(plasma::PlasmaParams, k::Real, θ::Real=0.0; model::Symbol=:hall, amp::Real=0.15)

Construct a right-hand polarized Whistler wave.
- `model = :hall`: exact Hall-MHD whistler (for `θ = 0`) or oblique EMHD whistler
  ``\\omega = \\Omega_e (k d_e)^2 \\cos\\theta / (1 + k^2 d_e^2)`` when `k * d_i >> 1`.
- `model = :cold`: solves the full cold-plasma dispersion relation for the whistler branch
  below the electron cyclotron resonance cone ``\\omega < |\\Omega_e| \\cos\\theta``.
"""
function whistler_wave(
    plasma::PlasmaParams,
    k::Real,
    θ::Real = 0.0;
    model::Symbol = :hall,
    amp::Real = 0.15,
)
    Ωe = abs(cyclotron_frequency(get_electron(plasma), plasma))
    ω_max = 0.999 * Ωe * abs(cos(θ))
    ω_max <= 0 && error("Whistler cannot propagate at θ = π/2 (outside resonance cone).")

    if model === :hall && abs(θ) < 1e-12
        ω = min(hall_mhd_dispersion(plasma, k; branch = :whistler), 0.95 * Ωe)
        vA = alfven_speed(plasma)
        kvec = SVector(0.0, 0.0, Float64(k))
        B1 = amp * plasma.B0
        # Right-hand circular polarization around +B0
        δB = SVector(ComplexF64(B1), ComplexF64(0, B1), ComplexF64(0))
        # Faraday: k × δE = ω δB -> δE = (ω/k) (δB × ẑ) = (ω/k) (δBy, -δBx, 0)
        vph = ω / k
        δE = SVector(vph * δB[2], -vph * δB[1], ComplexF64(0))
        # Hall-MHD Walén relation: δv = -(k vA² / (ω B0)) δB
        δv = -(k * vA^2 / (ω * plasma.B0)) * δB
        return WaveSolution(:whistler, ComplexF64(ω), kvec, δE, δB, δv, ComplexF64(0), ComplexF64(0), plasma)
    end

    # Solve full cold plasma R-whistler branch for given (k, θ) via bisection on n²(ω) - (ck/ω)² = 0
    f_res(ω) = _whistler_n2(plasma, ω, θ) - (plasma.c * k / ω)^2
    ω_lo = 1e-4 * abs(cyclotron_frequency(get_ion(plasma), plasma))
    ω_hi = ω_max
    ω = _bisect(f_res, ω_lo, ω_hi)
    n2 = _whistler_n2(plasma, ω, θ)
    return cold_plasma_eigenmode(plasma, ω, θ, n2; mode = :whistler, amp)
end

"""
    ion_cyclotron_wave(plasma::PlasmaParams, k::Real; model::Symbol=:hall, amp::Real=0.15)

Construct a left-hand polarized Alfvén–Ion Cyclotron Wave (EMIC) propagating along ``\\mathbf{B}_0``.
"""
function ion_cyclotron_wave(plasma::PlasmaParams, k::Real; model::Symbol = :hall, amp::Real = 0.15)
    Ωi = abs(cyclotron_frequency(get_ion(plasma), plasma))
    if model === :hall
        ω = hall_mhd_dispersion(plasma, k; branch = :emic)
        vA = alfven_speed(plasma)
        kvec = SVector(0.0, 0.0, Float64(k))
        B1 = amp * plasma.B0
        # Left-hand circular polarization around +B0
        δB = SVector(ComplexF64(B1), ComplexF64(0, -B1), ComplexF64(0))
        vph = ω / k
        δE = SVector(vph * δB[2], -vph * δB[1], ComplexF64(0))
        δv = -(k * vA^2 / (ω * plasma.B0)) * δB
        return WaveSolution(:ion_cyclotron, ComplexF64(ω), kvec, δE, δB, δv, ComplexF64(0), ComplexF64(0), plasma)
    else
        f_res(ω) = real(stix_parameters(plasma, ω).L) - (plasma.c * k / ω)^2
        ω = _bisect(f_res, 1e-5 * Ωi, 0.9999 * Ωi)
        n2 = real(stix_parameters(plasma, ω).L)
        return cold_plasma_eigenmode(plasma, ω, 0.0, n2; mode = :ion_cyclotron, amp)
    end
end

"Ordinary (O-mode) electromagnetic wave propagating perpendicular to ``\\mathbf{B}_0`` (``n^2 = P``)."
function o_mode_wave(plasma::PlasmaParams, k::Real; amp::Real = 0.15)
    ωpe = plasma_frequency(get_electron(plasma), plasma)
    ω = sqrt(ωpe^2 + (plasma.c * k)^2)
    n2 = (plasma.c * k / ω)^2
    return cold_plasma_eigenmode(plasma, ω, π / 2, n2; mode = :o_mode, amp)
end

"Extraordinary (X-mode) electromagnetic wave propagating perpendicular to ``\\mathbf{B}_0`` (``n^2 = RL/S``)."
function x_mode_wave(plasma::PlasmaParams, ω::Real; amp::Real = 0.15)
    (; R, L, S) = stix_parameters(plasma, ω)
    n2 = real(R * L / S)
    return cold_plasma_eigenmode(plasma, ω, π / 2, n2; mode = :x_mode, amp)
end

"Upper hybrid electrostatic oscillation (``\\mathbf{k} \\perp \\mathbf{B}_0``, ``\\omega^2 = \\omega_{UH}^2 + k^2 v_{th,e}^2``)."
function upper_hybrid_wave(plasma::PlasmaParams, k_perp::Real; amp::Real = 0.15)
    e = get_electron(plasma)
    ω_UH = cold_resonances(plasma).ω_UH
    vte = thermal_speed_1d(e, plasma; perp = true)
    ω = sqrt(ω_UH^2 + (k_perp * vte)^2)
    kvec = SVector(Float64(k_perp), 0.0, 0.0)
    E1 = amp * plasma.B0 * alfven_speed(plasma)
    δE = SVector(ComplexF64(E1), ComplexF64(0), ComplexF64(0))
    δB = SVector(ComplexF64(0), ComplexF64(0), ComplexF64(0))
    Ωe = cyclotron_frequency(e, plasma)
    # Electron velocity from cold equation of motion
    fac = im * e.q / (e.m * ω)
    denom = 1 - (Ωe / ω)^2
    vex = fac * E1 / denom
    vey = fac * (-im * (Ωe / ω) * E1) / denom
    δv = SVector(vex, vey, ComplexF64(0))
    δn = e.n * k_perp * vex / ω
    return WaveSolution(:upper_hybrid, ComplexF64(ω), kvec, δE, δB, δv, δn, ComplexF64(0), plasma)
end

"Lower hybrid electrostatic wave (``\\mathbf{k} \\perp \\mathbf{B}_0``, ``\\omega \\approx \\omega_{LH}``)."
function lower_hybrid_wave(plasma::PlasmaParams, k_perp::Real; amp::Real = 0.15)
    ω_LH = cold_resonances(plasma).ω_LH
    kvec = SVector(Float64(k_perp), 0.0, 0.0)
    E1 = amp * plasma.B0 * alfven_speed(plasma)
    δE = SVector(ComplexF64(E1), ComplexF64(0), ComplexF64(0))
    δB = SVector(ComplexF64(0), ComplexF64(0), ComplexF64(0))
    ion = get_ion(plasma)
    # Ions nearly unmagnetized along x (ω >> Ωi), electrons E×B drifting along y
    vix = im * ion.q * E1 / (ion.m * ω_LH)
    δv = SVector(vix, ComplexF64(0), ComplexF64(0))
    δn = ion.n * k_perp * vix / ω_LH
    return WaveSolution(:lower_hybrid, ComplexF64(ω_LH), kvec, δE, δB, δv, δn, ComplexF64(0), plasma)
end

"""
    faraday_rotation_fields(plasma::PlasmaParams, ω::Real, z::AbstractVector, t::Real; E0=1.0)

Superposition of equal-amplitude R- and L-mode waves at frequency ``\\omega > \\omega_R``
propagating along ``+\\mathbf{B}_0``: demonstrates Faraday rotation of the linear
polarization plane with rotation rate ``d\\phi/dz = (k_L - k_R)/2``.
Returns `(; Ex, Ey, kR, kL, rotation_rate)`.
"""
function faraday_rotation_fields(plasma::PlasmaParams, ω::Real, z::AbstractVector, t::Real; E0::Real = 1.0)
    (; R, L) = stix_parameters(plasma, ω)
    Rr, Lr = real(R), real(L)
    (Rr <= 0 || Lr <= 0) && error("Faraday rotation requires ω above both R and L cutoffs (R=$Rr, L=$Lr).")
    kR = ω * sqrt(Rr) / plasma.c
    kL = ω * sqrt(Lr) / plasma.c

    # E_R = (E0/2)(1, i) e^{i(kR z - ωt)}, E_L = (E0/2)(1, -i) e^{i(kL z - ωt)}
    Ex = [real(0.5 * E0 * (exp(im * (kR * zi - ω * t)) + exp(im * (kL * zi - ω * t)))) for zi in z]
    Ey = [real(0.5 * E0 * im * (exp(im * (kR * zi - ω * t)) - exp(im * (kL * zi - ω * t)))) for zi in z]
    return (; Ex, Ey, kR, kL, rotation_rate = 0.5 * (kL - kR))
end

# Internal bisection helper for monotone dispersion brackets
function _bisect(f, a::Real, b::Real; tol::Real = 1e-11, maxiter::Int = 100)
    fa = f(a)
    fb = f(b)
    if fa * fb > 0
        # Scan logarithmically for a sign change inside [a, b]
        xs = exp.(range(log(a), log(b), length = 200))
        prev_x, prev_f = xs[1], f(xs[1])
        found = false
        for xi in xs[2:end]
            fi = f(xi)
            if isfinite(prev_f) && isfinite(fi) && prev_f * fi <= 0
                a, b, fa, fb = prev_x, xi, prev_f, fi
                found = true
                break
            end
            prev_x, prev_f = xi, fi
        end
        found || error("Could not bracket root in [$a, $b].")
    end
    lo, hi = float(a), float(b)
    flo = fa
    for _ in 1:maxiter
        mid = 0.5 * (lo + hi)
        fmid = f(mid)
        if abs(hi - lo) <= tol * max(1.0, abs(mid)) || abs(fmid) <= tol
            return mid
        end
        if flo * fmid <= 0
            hi = mid
        else
            lo = mid
            flo = fmid
        end
    end
    return 0.5 * (lo + hi)
end
