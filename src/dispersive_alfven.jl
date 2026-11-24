"""
Dispersive Alfvén Waves:
- Kinetic Alfvén Wave (KAW, β_e >> m_e/m_i): two-fluid Padé & gyrokinetic Bessel dispersion
- Inertial Alfvén Wave (IAW, β_e << m_e/m_i): electron-inertia dispersion
- Parallel test-particle trapping in KAW's finite E_parallel vs. free streaming in ideal SAW
"""

using SpecialFunctions: besselix

"""
    kinetic_alfven_wave(plasma::PlasmaParams, k_par::Real, k_perp::Real;
                        model::Symbol=:bessel, amp::Real=0.15)

Construct a Kinetic Alfvén Wave (KAW) in the warm-electron regime ``\\beta_e \\gg m_e / m_i``.
- `model = :bessel`: gyrokinetic modified-Bessel dispersion (Johnson & Cheng 1997)
  ``\\omega^2 = k_\\parallel^2 v_A^2 \\left[\\frac{b_i}{1 - I_0(b_i)e^{-b_i}} + \\frac{T_e}{T_i} b_i\\right]``
  with ``b_i = k_\\perp^2 r_{iL}^2``.
- `model = :twofluid`: Padé / two-fluid approximation
  ``\\omega^2 = k_\\parallel^2 v_A^2 [1 + k_\\perp^2 (r_{iL}^2 + r_{fL}^2)]``.

Computes the finite parallel electric field ``\\delta E_z``, modified ``|\\delta E_\\perp / \\delta B_\\perp|``,
and Boltzmann density perturbation ``\\delta n``.
"""
function kinetic_alfven_wave(
    plasma::PlasmaParams,
    k_par::Real,
    k_perp::Real;
    model::Symbol = :bessel,
    amp::Real = 0.15,
)
    vA = alfven_speed(plasma)
    ion = get_ion(plasma)
    elec = get_electron(plasma)
    Ωi = abs(cyclotron_frequency(ion, plasma))

    riL = larmor_radius(ion, plasma)
    rfL = sqrt(plasma.kB * elec.T_par / ion.m) / Ωi # ion-acoustic gyroradius ρ_s
    bi = (k_perp * riL)^2
    τ = elec.T_par / max(ion.T_perp, 1e-30)

    if model === :bessel && bi > 1e-7
        Γ0 = besselix(0, bi) # I_0(bi) * exp(-bi), overflow-safe
        disp_factor = bi / (1.0 - Γ0) + τ * bi
    else
        disp_factor = 1.0 + k_perp^2 * (riL^2 + rfL^2)
    end
    ω = abs(k_par) * vA * sqrt(disp_factor)

    kvec = SVector(Float64(k_perp), 0.0, Float64(k_par))
    # Perturbed magnetic field in y (shear mode)
    By = ComplexF64(amp * plasma.B0)
    δB = SVector(ComplexF64(0), By, ComplexF64(0))

    # Polarization: |δE_x / δB_y| = v_A (1 + k_perp^2 r_iL^2) / sqrt(disp_factor)
    sgn = k_par >= 0 ? 1.0 : -1.0
    Ex = sgn * By * vA * (1.0 + bi) / sqrt(disp_factor)
    # Parallel electric field (distinctive KAW feature!):
    # δE_z / δE_x = - k_par k_perp r_fL^2 / (1 + k_perp^2 r_iL^2)
    Ez = -Ex * (k_par * k_perp * rfL^2) / (1.0 + bi)
    δE = SVector(Ex, ComplexF64(0), Ez)

    # E×B drift + ion polarization contribution (leading order v_y = -E_x / B0)
    δv = SVector(ComplexF64(0), -Ex / plasma.B0, ComplexF64(0))
    # Boltzmann electron density response: δn / n0 = i e δE_z / (k_par k_B T_e)
    δn = if abs(k_par) > 1e-14 && elec.T_par > 0
        elec.n * (im * abs(elec.q) * Ez) / (k_par * plasma.kB * elec.T_par)
    else
        ComplexF64(0)
    end
    δp = plasma.kB * (elec.T_par + ion.T_par) * δn

    return WaveSolution(:kaw, ComplexF64(ω), kvec, δE, δB, δv, δn, δp, plasma)
end

"""
    inertial_alfven_wave(plasma::PlasmaParams, k_par::Real, k_perp::Real; amp::Real=0.15)

Construct an Inertial Alfvén Wave (IAW) in the cold-electron regime ``\\beta_e \\ll m_e / m_i``:
``\\omega^2 = k_\\parallel^2 v_A^2 / (1 + k_\\perp^2 d_e^2)``.
"""
function inertial_alfven_wave(plasma::PlasmaParams, k_par::Real, k_perp::Real; amp::Real = 0.15)
    vA = alfven_speed(plasma)
    elec = get_electron(plasma)
    de = inertial_length(elec, plasma)
    kperp_de2 = (k_perp * de)^2

    ω = abs(k_par) * vA / sqrt(1.0 + kperp_de2)
    kvec = SVector(Float64(k_perp), 0.0, Float64(k_par))
    By = ComplexF64(amp * plasma.B0)
    δB = SVector(ComplexF64(0), By, ComplexF64(0))

    sgn = k_par >= 0 ? 1.0 : -1.0
    Ex = sgn * By * vA * sqrt(1.0 + kperp_de2)
    Ez = -Ex * (k_par * k_perp * de^2) / (1.0 + kperp_de2)
    δE = SVector(Ex, ComplexF64(0), Ez)
    δv = SVector(ComplexF64(0), -Ex / plasma.B0, ComplexF64(0))

    return WaveSolution(:iaw, ComplexF64(ω), kvec, δE, δB, δv, ComplexF64(0), ComplexF64(0), plasma)
end

"""
    kaw_particle_trajectories(sol::WaveSolution, z0s, vz0s, ts; q_over_m=nothing)

Integrate 1D parallel electron motion ``\\ddot{z} = (q/m) E_z(z, t)`` in the wave's
parallel electric field using 4th-order Runge-Kutta.
For ideal SAW (``E_z = 0``) particles stream freely; for KAW (``E_z \\neq 0``) resonant
particles near ``v_z \\approx \\omega / k_\\parallel`` become trapped in the wave potential troughs.

Returns `(; z, vz)` matrices of size `(length(z0s), length(ts))`.
"""
function kaw_particle_trajectories(
    sol::WaveSolution,
    z0s::AbstractVector,
    vz0s::AbstractVector,
    ts::AbstractVector;
    q_over_m::Union{Nothing,Real} = nothing,
)
    qm = if q_over_m === nothing
        e = get_electron(sol.plasma)
        e.q / e.m
    else
        Float64(q_over_m)
    end
    kz = sol.k[3]
    ωr = real(sol.ω)
    Ez_amp = sol.δE[3]

    accel(z, t) = qm * real(Ez_amp * exp(im * (kz * z - ωr * t)))

    np = length(z0s)
    nt = length(ts)
    Z = Matrix{Float64}(undef, np, nt)
    VZ = Matrix{Float64}(undef, np, nt)
    Z[:, 1] .= z0s
    VZ[:, 1] .= vz0s

    for it in 1:(nt - 1)
        t = ts[it]
        dt = ts[it + 1] - t
        for ip in 1:np
            z = Z[ip, it]
            v = VZ[ip, it]

            k1z = v
            k1v = accel(z, t)

            k2z = v + 0.5 * dt * k1v
            k2v = accel(z + 0.5 * dt * k1z, t + 0.5 * dt)

            k3z = v + 0.5 * dt * k2v
            k3v = accel(z + 0.5 * dt * k2z, t + 0.5 * dt)

            k4z = v + dt * k3v
            k4v = accel(z + dt * k3z, t + dt)

            Z[ip, it + 1] = z + (dt / 6) * (k1z + 2k2z + 2k3z + k4z)
            VZ[ip, it + 1] = v + (dt / 6) * (k1v + 2k2v + 2k3v + k4v)
        end
    end
    return (; z = Z, vz = VZ)
end
