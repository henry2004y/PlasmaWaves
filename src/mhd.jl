"""
Ideal MHD wave modes:
- Linear Shear Alfvén Wave (SAW)
- Exact Nonlinear Circularly Polarized Alfvén Wave (CPAW)
- Fast and Slow Magnetosonic Waves
- Friedrich phase & group velocity diagrams
- Non-uniform slab SAW phase mixing
"""

"""
    friedrich_phase_speeds(vA, vs, θ)

Return `(; v_fast, v_alfven, v_slow)` ideal MHD phase speeds at propagation angle `θ`
relative to ``\\mathbf{B}_0``.
"""
function friedrich_phase_speeds(vA::Real, vs::Real, θ::Real)
    cθ = cos(θ)
    sum2 = vA^2 + vs^2
    disc = sqrt(max(0.0, sum2^2 - 4 * vA^2 * vs^2 * cθ^2))
    v_fast = sqrt(0.5 * (sum2 + disc))
    v_slow = sqrt(max(0.0, 0.5 * (sum2 - disc)))
    v_alfven = vA * abs(cθ)
    return (; v_fast, v_alfven, v_slow)
end

"""
    friedrich_group_velocities(vA, vs, θ)

Return `(; vg_fast, vg_alfven, vg_slow)` as 2-vectors `(v_gx, v_gz)` in the ``(x, z)`` plane
using the exact angular derivative ``\\mathbf{v}_g = v_{\\mathrm{ph}}\\hat{\\mathbf{k}} + \\frac{\\partial v_{\\mathrm{ph}}}{\\partial \\theta}\\hat{\\boldsymbol{\\theta}}``.
"""
function friedrich_group_velocities(vA::Real, vs::Real, θ::Real)
    sθ, cθ = sin(θ), cos(θ)
    khat = SVector(sθ, cθ)
    θhat = SVector(cθ, -sθ)

    # Central finite difference for smooth derivative, exact for Alfvén
    δ = 1e-6
    sp0 = friedrich_phase_speeds(vA, vs, θ)
    sp_p = friedrich_phase_speeds(vA, vs, θ + δ)
    sp_m = friedrich_phase_speeds(vA, vs, θ - δ)

    dvf = (sp_p.v_fast - sp_m.v_fast) / (2δ)
    dvs = (sp_p.v_slow - sp_m.v_slow) / (2δ)

    vg_fast = sp0.v_fast * khat + dvf * θhat
    vg_slow = sp0.v_slow * khat + dvs * θhat
    vg_alfven = SVector(0.0, sign(cθ == 0 ? 1.0 : cθ) * vA)
    return (; vg_fast, vg_alfven, vg_slow)
end

"""
    shear_alfven_wave(plasma::PlasmaParams, k::Real, θ::Real=0.0; amp::Real=0.15)

Construct a linearly polarized Shear Alfvén Wave (SAW) with ``\\mathbf{k} = k(\\sin\\theta, 0, \\cos\\theta)``,
``\\omega = |k_\\parallel| v_A``, transverse ``\\delta v_y, \\delta B_y, \\delta E_x``, and zero compression.
"""
function shear_alfven_wave(plasma::PlasmaParams, k::Real, θ::Real = 0.0; amp::Real = 0.15)
    vA = alfven_speed(plasma)
    kvec = SVector(Float64(k * sin(θ)), 0.0, Float64(k * cos(θ)))
    kz = kvec[3]
    ω = abs(kz) * vA
    sgn = kz >= 0 ? 1.0 : -1.0

    vy = ComplexF64(amp * vA)
    δv = SVector(ComplexF64(0), vy, ComplexF64(0))
    δE = SVector(-plasma.B0 * vy, ComplexF64(0), ComplexF64(0))
    δB = SVector(ComplexF64(0), -sgn * plasma.B0 * (vy / vA), ComplexF64(0))

    return WaveSolution(
        :shear_alfven,
        ComplexF64(ω),
        kvec,
        δE,
        δB,
        δv,
        ComplexF64(0),
        ComplexF64(0),
        plasma,
    )
end

"""
    circular_alfven_wave(plasma::PlasmaParams, k::Real; hand::Symbol=:right, amp::Real=0.2)

Construct an exact arbitrary-amplitude circularly polarized Alfvén wave (CPAW)
propagating along ``+\\mathbf{B}_0`` with ``|\\mathbf{B}| = \\text{const}``.
`hand` can be `:right` (electron gyration sense) or `:left` (ion gyration sense).
"""
function circular_alfven_wave(plasma::PlasmaParams, k::Real; hand::Symbol = :right, amp::Real = 0.2)
    vA = alfven_speed(plasma)
    kvec = SVector(0.0, 0.0, Float64(k))
    ω = abs(k) * vA
    h = hand === :right ? 1.0 : hand === :left ? -1.0 : error("hand must be :right or :left")

    B1 = amp * plasma.B0
    # Re[(1, h*i) exp(i(kz - ωt))] = (cos(kz - ωt), -h*sin(kz - ωt))
    # At z=0: (cos(ωt), h*sin(ωt)) -> right-handed around +z for h = +1
    δB = SVector(ComplexF64(B1), ComplexF64(0, h * B1), ComplexF64(0))
    # Parallel propagation along +B0: δv = -(vA/B0) δB
    sgn = k >= 0 ? 1.0 : -1.0
    δv = -sgn * (vA / plasma.B0) * δB
    # δE = -δv × B0
    δE = SVector(-δv[2] * plasma.B0, δv[1] * plasma.B0, ComplexF64(0))

    return WaveSolution(
        hand === :right ? :cpaw_right : :cpaw_left,
        ComplexF64(ω),
        kvec,
        δE,
        δB,
        δv,
        ComplexF64(0),
        ComplexF64(0),
        plasma,
    )
end

function _magnetosonic_wave(plasma::PlasmaParams, k::Real, θ::Real, branch::Symbol; amp::Real = 0.15)
    vA = alfven_speed(plasma)
    vs = sound_speed(plasma)
    speeds = friedrich_phase_speeds(vA, vs, θ)
    vph = branch === :fast ? speeds.v_fast : speeds.v_slow
    ω = k * vph

    sθ, cθ = sin(θ), cos(θ)
    kvec = SVector(Float64(k * sθ), 0.0, Float64(k * cθ))

    # Eigenvector (vx, vz) of the 2x2 x-z ideal MHD block:
    # [ -vph^2 + vA^2 + vs^2 sin^2 θ,  vs^2 sin θ cos θ ] [vx] = 0
    # [ vs^2 sin θ cos θ,              -vph^2 + vs^2 cos^2 θ ] [vz]
    a11 = -vph^2 + vA^2 + vs^2 * sθ^2
    a12 = vs^2 * sθ * cθ
    a22 = -vph^2 + vs^2 * cθ^2

    if abs(a12) < 1e-12 * max(vA^2, vs^2, 1e-30)
        # Parallel (θ=0) or perpendicular (θ=π/2)
        if abs(a11) <= abs(a22)
            vx_unit, vz_unit = 1.0, 0.0
        else
            vx_unit, vz_unit = 0.0, 1.0
        end
    else
        # Choose row with larger norm for numerical stability
        if hypot(a11, a12) >= hypot(a12, a22)
            vx_unit, vz_unit = -a12, a11
        else
            vx_unit, vz_unit = a22, -a12
        end
        nrm = hypot(vx_unit, vz_unit)
        vx_unit /= nrm
        vz_unit /= nrm
    end

    # Normalize amplitude by amp * max(vA, vs)
    vscale = amp * max(vA, vs)
    vx = ComplexF64(vscale * vx_unit)
    vz = ComplexF64(vscale * vz_unit)
    δv = SVector(vx, ComplexF64(0), vz)

    ρ0 = mass_density(plasma)
    n0 = get_ion(plasma).n
    kdotv = kvec[1] * vx + kvec[3] * vz

    if ω > 1e-14
        δρ = (ρ0 / ω) * kdotv
        δn = (n0 / ω) * kdotv
        δp = vs^2 * δρ
        δE = SVector(ComplexF64(0), plasma.B0 * vx, ComplexF64(0))
        δB = SVector(-plasma.B0 * kvec[3] * vx / ω, ComplexF64(0), plasma.B0 * kvec[1] * vx / ω)
    else
        # Perpendicular slow / mirror limit (vph -> 0): finite δp + δB_z B0/μ0 = 0
        δp = ComplexF64(amp * thermal_pressure(plasma))
        δn = ComplexF64(amp * n0)
        δE = SVector(ComplexF64(0), ComplexF64(0), ComplexF64(0))
        δB = SVector(ComplexF64(0), ComplexF64(0), -plasma.μ0 * δp / plasma.B0)
    end

    return WaveSolution(
        branch === :fast ? :fast_magnetosonic : :slow_magnetosonic,
        ComplexF64(ω),
        kvec,
        δE,
        δB,
        δv,
        δn,
        δp,
        plasma,
    )
end

"Construct a Fast Magnetosonic wave (`δp` and `δB_z` in phase)."
fast_magnetosonic_wave(plasma::PlasmaParams, k::Real, θ::Real = π / 3; amp::Real = 0.15) =
    _magnetosonic_wave(plasma, k, θ, :fast; amp)

"Construct a Slow Magnetosonic wave (`δp` and `δB_z` anti-phase)."
slow_magnetosonic_wave(plasma::PlasmaParams, k::Real, θ::Real = π / 4; amp::Real = 0.15) =
    _magnetosonic_wave(plasma, k, θ, :slow; amp)

"""
    saw_phase_mixing(x, t; kz=1.0, vA_profile = x -> 1.0 + 0.4 * tanh(x),
                     B0_profile = x -> exp(-x^2 / 4))

Exact 1D non-uniform slab Shear Alfvén Wave continuous-spectrum solution:
``B_{1y}(x, t) = \\hat{B}_{1y}(x, 0) \\cos[\\omega_A(x) t]`` where ``\\omega_A(x) = k_z v_A(x)``.
Demonstrates phase mixing where ``|k_x(x, t)| \\simeq |\\omega_A'(x)| t`` grows linearly in time.
"""
function saw_phase_mixing(
    x::AbstractVector,
    t::Real;
    kz::Real = 1.0,
    vA_profile = x -> 1.0 + 0.5 * tanh(x),
    B0_profile = x -> exp(-x^2 / 4.0),
)
    return [B0_profile(xi) * cos(kz * vA_profile(xi) * t) for xi in x]
end
