"""
Warm-fluid & analytical kinetic electrostatic plasma wave theory:
- Plasma dispersion function Z(ζ) and Z'(ζ) via Faddeeva function erfcx
- 1D unmagnetized Vlasov electrostatic dispersion (Langmuir/Bohm-Gross + Landau damping, Ion-Acoustic)
- Landau damping (x, v) phase-space test-particle trajectories
- Perpendicular Electrostatic Bernstein Waves (EBW / IBW)

For general hot magnetized electromagnetic kinetic dispersion relations (including arbitrary
harmonics, anisotropies, and non-Maxwellian VDFs), use `PlasmaBO.jl` via `confirm_dispersion`.
"""

using SpecialFunctions: erfcx, besselix

"""
    plasma_dispersion_Z(ζ::Number)

Fried-Conte plasma dispersion function
``Z(\\zeta) = \\frac{1}{\\sqrt{\\pi}} \\int_{-\\infty}^{\\infty} \\frac{e^{-t^2}}{t - \\zeta}\\,dt = i\\sqrt{\\pi}\\,\\mathrm{erfcx}(-i\\zeta)``,
valid across the entire complex plane.
"""
plasma_dispersion_Z(ζ::Number) = im * sqrt(π) * erfcx(-im * complex(float(ζ)))

"Derivative of the plasma dispersion function: ``Z'(\\zeta) = -2[1 + \\zeta Z(\\zeta)]``."
plasma_dispersion_dZ(ζ::Number) = -2 * (1 + complex(float(ζ)) * plasma_dispersion_Z(ζ))

# Lightweight complex secant solver for 1D electrostatic Vlasov dispersion
function _complex_secant(f, z0::ComplexF64; tol::Real = 1e-11, maxiter::Int = 80)
    z1 = z0 * (1.0 + 1e-4) + 1e-6 * (1.0 + im)
    f0 = complex(f(z0))
    f1 = complex(f(z1))
    abs(f0) < tol && return z0
    abs(f1) < tol && return z1

    for _ in 1:maxiter
        df = f1 - f0
        abs(df) < 1e-30 && break
        dz = -f1 * (z1 - z0) / df
        maxstep = max(0.5 * abs(z1), 0.5)
        abs(dz) > maxstep && (dz *= maxstep / abs(dz))
        z2 = z1 + dz
        f2 = complex(f(z2))
        if abs(dz) <= tol * max(1.0, abs(z2)) && abs(f2) <= 1e-8
            return z2
        end
        z0, f0 = z1, f1
        z1, f1 = z2, f2
    end
    return z1
end

"""
    electrostatic_dielectric(plasma::PlasmaParams, ω::Number, k::Real; species_list=plasma.species)

Unmagnetized Vlasov electrostatic dielectric function:
``\\epsilon(\\omega, k) = 1 + \\sum_s \\frac{1}{k^2 \\lambda_{Ds}^2} [1 + \\zeta_s Z(\\zeta_s)]``,
where ``\\zeta_s = \\omega / (|k| v_{th,s})`` and ``v_{th,s} = \\sqrt{2 k_B T_{\\parallel s} / m_s}``.
"""
function electrostatic_dielectric(
    plasma::PlasmaParams,
    ω::Number,
    k::Real;
    species_list = plasma.species,
)
    ωc = complex(float(ω))
    absk = abs(k)
    ϵ = one(ComplexF64)
    for s in species_list
        s.T_par <= 0 && continue
        λD = debye_length(s, plasma)
        vth = thermal_speed(s, plasma)
        ζ = ωc / (absk * vth)
        ϵ += (1.0 / (absk * λD)^2) * (1.0 + ζ * plasma_dispersion_Z(ζ))
    end
    return ϵ
end

"""
    langmuir_wave(plasma::PlasmaParams, k::Real; model::Symbol=:vlasov, amp::Real=0.15)

Construct a Langmuir (electron plasma) wave along ``\\hat{\\mathbf{x}}``.
- `model = :vlasov`: solves the 1D Vlasov electrostatic dispersion relation for complex ``\\omega = \\omega_r + i\\gamma``.
- `model = :bohm_gross`: warm-fluid Bohm-Gross approximation ``\\omega_r = \\sqrt{\\omega_{pe}^2 + 3 k^2 v_{th,1D}^2}``
  with asymptotic weak Landau damping rate.
"""
function langmuir_wave(plasma::PlasmaParams, k::Real; model::Symbol = :vlasov, amp::Real = 0.15)
    e = get_electron(plasma)
    ωpe = plasma_frequency(e, plasma)
    λDe = debye_length(e, plasma)
    vte_1d = thermal_speed_1d(e, plasma)
    kλ = abs(k) * λDe

    ω_bg = sqrt(ωpe^2 + 3 * (k * vte_1d)^2)
    γ_weak = kλ > 1e-3 ? -ωpe * sqrt(π / 8) / kλ^3 * exp(-1 / (2 * kλ^2) - 1.5) : 0.0

    ω = if model === :bohm_gross
        ComplexF64(ω_bg, γ_weak)
    else
        f(w) = electrostatic_dielectric(plasma, w, k; species_list = [e])
        guess = ComplexF64(ω_bg, min(γ_weak, -1e-4 * ωpe))
        _complex_secant(f, guess)
    end

    kvec = SVector(Float64(k), 0.0, 0.0)
    E1 = ComplexF64(amp * plasma.B0 * alfven_speed(plasma))
    δE = SVector(E1, ComplexF64(0), ComplexF64(0))
    δB = SVector(ComplexF64(0), ComplexF64(0), ComplexF64(0))
    δn = -im * plasma.ε0 * k * E1 / abs(e.q)
    δv = SVector((ω / (k * e.n)) * δn, ComplexF64(0), ComplexF64(0))
    δp = 3 * plasma.kB * e.T_par * δn

    return WaveSolution(:langmuir, ω, kvec, δE, δB, δv, δn, δp, plasma)
end

"""
    ion_acoustic_wave(plasma::PlasmaParams, k::Real; model::Symbol=:vlasov, amp::Real=0.15)

Construct an electrostatic Ion-Acoustic Wave along ``\\hat{\\mathbf{x}}``.
- `model = :vlasov`: solves the 1D two-species Vlasov electrostatic dispersion relation.
- `model = :fluid`: warm-fluid dispersion ``\\omega^2 = k^2 c_s^2 / (1 + k^2 \\lambda_{De}^2) + 3 k^2 v_{th,i,1D}^2``.
"""
function ion_acoustic_wave(plasma::PlasmaParams, k::Real; model::Symbol = :vlasov, amp::Real = 0.15)
    e = get_electron(plasma)
    ion = get_ion(plasma)
    λDe = debye_length(e, plasma)
    cs = sqrt(plasma.kB * e.T_par / ion.m)
    vti_1d = thermal_speed_1d(ion, plasma)
    ω_fluid = abs(k) * sqrt(cs^2 / (1 + (k * λDe)^2) + 3 * vti_1d^2)

    ω = if model === :fluid
        ComplexF64(ω_fluid, 0.0)
    else
        f(w) = electrostatic_dielectric(plasma, w, k)
        τ = e.T_par / max(ion.T_par, 1e-30)
        γ0 = -ω_fluid * (τ > 5 ? 0.03 : 0.25)
        _complex_secant(f, ComplexF64(ω_fluid, γ0))
    end

    kvec = SVector(Float64(k), 0.0, 0.0)
    E1 = ComplexF64(amp * plasma.B0 * alfven_speed(plasma))
    δE = SVector(E1, ComplexF64(0), ComplexF64(0))
    δB = SVector(ComplexF64(0), ComplexF64(0), ComplexF64(0))
    δn = im * e.n * abs(e.q) * E1 / (k * plasma.kB * e.T_par)
    δv = SVector((ω / (k * ion.n)) * δn, ComplexF64(0), ComplexF64(0))
    δp = plasma.kB * (e.T_par + 3 * ion.T_par) * δn

    return WaveSolution(:ion_acoustic, ω, kvec, δE, δB, δv, δn, δp, plasma)
end

"""
    landau_phase_space(sol::WaveSolution, x0s, v0s, ts; q_over_m=nothing, include_damping=true)

Integrate test-particle trajectories ``\\ddot{x} = (q/m) E_x(x, t)`` in an electrostatic wave
to illustrate Landau damping and cat's-eye phase-space trapping around ``v_\\phi = \\Re(\\omega)/k``.
"""
function landau_phase_space(
    sol::WaveSolution,
    x0s::AbstractVector,
    v0s::AbstractVector,
    ts::AbstractVector;
    q_over_m::Union{Nothing,Real} = nothing,
    include_damping::Bool = true,
)
    qm = if q_over_m === nothing
        e = get_electron(sol.plasma)
        e.q / e.m
    else
        Float64(q_over_m)
    end
    kx = sol.k[1]
    ω_eff = include_damping ? sol.ω : ComplexF64(real(sol.ω), 0.0)
    Ex_amp = sol.δE[1]

    accel(x, t) = qm * real(Ex_amp * exp(im * (kx * x - ω_eff * t)))

    np = length(x0s)
    nt = length(ts)
    X = Matrix{Float64}(undef, np, nt)
    V = Matrix{Float64}(undef, np, nt)
    X[:, 1] .= x0s
    V[:, 1] .= v0s

    for it in 1:(nt - 1)
        t = ts[it]
        dt = ts[it + 1] - t
        for ip in 1:np
            x = X[ip, it]
            v = V[ip, it]

            k1x = v
            k1v = accel(x, t)
            k2x = v + 0.5 * dt * k1v
            k2v = accel(x + 0.5 * dt * k1x, t + 0.5 * dt)
            k3x = v + 0.5 * dt * k2v
            k3v = accel(x + 0.5 * dt * k2x, t + 0.5 * dt)
            k4x = v + dt * k3v
            k4v = accel(x + dt * k3x, t + dt)

            X[ip, it + 1] = x + (dt / 6) * (k1x + 2k2x + 2k3x + k4x)
            V[ip, it + 1] = v + (dt / 6) * (k1v + 2k2v + 2k3v + k4v)
        end
    end
    return (; x = X, v = V)
end

# --- Electrostatic Bernstein Waves (EBW / IBW) ---

"""
    bernstein_dielectric(plasma::PlasmaParams, ω::Real, k_perp::Real;
                         species::Symbol=:electron, n_max::Int=25)

Exact perpendicular electrostatic Bernstein dispersion function (``k_\\parallel = 0``):
``\\epsilon(\\omega, k_\\perp) = 1 - \\sum_s \\frac{2\\omega_{ps}^2}{\\lambda_s \\Omega_s^2}
\\sum_{n=1}^{n_{\\max}} \\frac{n^2 \\Omega_s^2}{\\omega^2 - n^2 \\Omega_s^2} \\Gamma_n(\\lambda_s)``.
"""
function bernstein_dielectric(
    plasma::PlasmaParams,
    ω::Real,
    k_perp::Real;
    species::Symbol = :electron,
    n_max::Int = 25,
)
    sp_list = if species === :all
        plasma.species
    elseif species === :electron
        [get_electron(plasma)]
    elseif species === :ion
        [get_ion(plasma)]
    else
        error("species must be :electron, :ion, or :all")
    end

    ϵ = 1.0
    ω2 = ω^2
    for s in sp_list
        ωps2 = plasma_frequency(s, plasma)^2
        Ωs = abs(cyclotron_frequency(s, plasma))
        rL = larmor_radius(s, plasma)
        λ = (k_perp * rL)^2
        if λ < 1e-10
            ϵ -= ωps2 / (ω2 - Ωs^2)
        else
            acc = 0.0
            for n in 1:n_max
                Γn = besselix(n, λ)
                nΩ2 = (n * Ωs)^2
                acc += (nΩ2 / (ω2 - nΩ2)) * Γn
            end
            ϵ -= (2 * ωps2 / (λ * Ωs^2)) * acc
        end
    end
    return ϵ
end

"""
    bernstein_wave(plasma::PlasmaParams, k_perp::Real;
                   species::Symbol=:electron, harmonic::Int=1, n_max::Int=25, amp::Real=0.15)

Solve for the `harmonic`-th perpendicular electrostatic Bernstein wave (EBW or IBW)
in the frequency band ``(n\\Omega_s, (n+1)\\Omega_s)``.
"""
function bernstein_wave(
    plasma::PlasmaParams,
    k_perp::Real;
    species::Symbol = :electron,
    harmonic::Int = 1,
    n_max::Int = 25,
    amp::Real = 0.15,
)
    sp = species === :electron ? get_electron(plasma) : get_ion(plasma)
    Ωs = abs(cyclotron_frequency(sp, plasma))

    f(ω) = bernstein_dielectric(plasma, ω, k_perp; species, n_max)
    δ = 1e-12 * Ωs
    ω_lo = harmonic * Ωs + δ
    ω_hi = (harmonic + 1) * Ωs - δ
    f_lo = f(ω_lo)
    f_hi = f(ω_hi)
    ω = if f_lo * f_hi <= 0
        _bisect(f, ω_lo, ω_hi; tol = 1e-12)
    elseif f_lo > 0 && f_hi > 0
        float(harmonic * Ωs)
    else
        ω_uh = cold_resonances(plasma).ω_UH
        ω_top = max(ω_hi, 1.2 * ω_uh)
        f(ω_top) >= 0 ? _bisect(f, ω_lo, ω_top; tol = 1e-12) : float((harmonic + 1) * Ωs)
    end

    kvec = SVector(Float64(k_perp), 0.0, 0.0)
    E1 = ComplexF64(amp * plasma.B0 * alfven_speed(plasma))
    δE = SVector(E1, ComplexF64(0), ComplexF64(0))
    δB = SVector(ComplexF64(0), ComplexF64(0), ComplexF64(0))
    δn = -im * plasma.ε0 * k_perp * E1 / sp.q
    δv = SVector((ω / (k_perp * sp.n)) * δn, ComplexF64(0), ComplexF64(0))

    mode = species === :electron ? :ebw : :ibw
    return WaveSolution(mode, ComplexF64(ω), kvec, δE, δB, δv, δn, ComplexF64(0), plasma)
end
