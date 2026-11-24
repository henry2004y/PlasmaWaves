"""
Core physical types, plasma parameter helpers, wave eigenmode representation,
transport ratios (Gary 1993), and space-time field evaluators.
"""

using LinearAlgebra
using StaticArrays

# Physical constants (SI)
const C_SI  = 2.99792458e8        # speed of light [m/s]
const EPS0_SI = 8.8541878128e-12  # vacuum permittivity [F/m]
const MU0_SI  = 4π * 1e-7         # vacuum permeability [H/m]
const KB_SI   = 1.380649e-23      # Boltzmann constant [J/K]
const QE_SI   = 1.602176634e-19   # elementary charge [C]
const ME_SI   = 9.1093837015e-31  # electron mass [kg]
const MP_SI   = 1.67262192369e-27 # proton mass [kg]
const EV_TO_J = QE_SI             # 1 eV in Joules

"""
    Species{T<:AbstractFloat}

A single charged particle species in a magnetized plasma.

Fields:
- `name`: species label (e.g. `:electron`, `:proton`)
- `m`: particle mass
- `q`: signed particle charge (`< 0` for electrons, `> 0` for ions)
- `n`: equilibrium number density
- `T_par`: parallel temperature (in Kelvin for SI, or energy units when `kB == 1`)
- `T_perp`: perpendicular temperature
"""
struct Species{T<:AbstractFloat}
    name::Symbol
    m::T
    q::T
    n::T
    T_par::T
    T_perp::T
end

function Species(name::Symbol, m::Real, q::Real, n::Real, T_par::Real, T_perp::Real=T_par)
    T = float(promote_type(typeof(m), typeof(q), typeof(n), typeof(T_par), typeof(T_perp)))
    return Species{T}(name, T(m), T(q), T(n), T(T_par), T(T_perp))
end

"""
    PlasmaParams{T<:AbstractFloat}

Uniform magnetized plasma with background magnetic field ``\\mathbf{B}_0 = B_0 \\hat{\\mathbf{z}}``.
Works seamlessly in either dimensionless normalized units ([`normalized_plasma`](@ref))
or SI units ([`si_plasma`](@ref)).
"""
struct PlasmaParams{T<:AbstractFloat}
    B0::T
    species::Vector{Species{T}}
    c::T
    ε0::T
    μ0::T
    kB::T
    γ::T
end

"""
    normalized_plasma(; mi_me=1836.0, β_i=0.1, β_e=0.1, vA_c=0.01,
                        anisotropy_i=1.0, anisotropy_e=1.0, γ=5/3)

Construct a two-species (proton + electron) plasma in standard Alfvén-normalized units:
``m_i = q_i = n_0 = B_0 = \\mu_0 = k_B = 1``, so that ``v_A = 1``, ``\\Omega_i = 1``, and ``d_i = 1``.

Plasma beta is defined as ``\\beta_s = 2\\mu_0 n_s k_B T_{\\parallel s} / B_0^2``,
and `anisotropy_s = T_perp / T_par`.
"""
function normalized_plasma(;
    mi_me::Real = 1836.0,
    β_i::Real = 0.1,
    β_e::Real = 0.1,
    vA_c::Real = 0.01,
    anisotropy_i::Real = 1.0,
    anisotropy_e::Real = 1.0,
    γ::Real = 5 / 3,
)
    T = Float64
    B0 = one(T)
    μ0 = one(T)
    kB = one(T)
    n0 = one(T)
    c = T(1 / vA_c)
    ε0 = T(1 / (μ0 * c^2))

    mi = one(T)
    me = T(1 / mi_me)
    Ti_par = T(β_i / 2)
    Te_par = T(β_e / 2)

    ions = Species(:proton, mi, one(T), n0, Ti_par, T(anisotropy_i * Ti_par))
    elec = Species(:electron, me, -one(T), n0, Te_par, T(anisotropy_e * Te_par))
    return PlasmaParams{T}(B0, [ions, elec], c, ε0, μ0, kB, T(γ))
end

"""
    si_plasma(; B0=1e-8, n0=5e6, Ti_eV=10.0, Te_eV=10.0, mi_me=MP_SI/ME_SI,
                anisotropy_i=1.0, anisotropy_e=1.0, γ=5/3)

Construct a proton-electron plasma in SI units (defaults approximate solar wind at 1 AU:
``B_0 = 10\\,\\mathrm{nT}``, ``n_0 = 5\\,\\mathrm{cm}^{-3}``, ``T = 10\\,\\mathrm{eV}``).
"""
function si_plasma(;
    B0::Real = 1e-8,
    n0::Real = 5e6,
    Ti_eV::Real = 10.0,
    Te_eV::Real = 10.0,
    mi_me::Real = MP_SI / ME_SI,
    anisotropy_i::Real = 1.0,
    anisotropy_e::Real = 1.0,
    γ::Real = 5 / 3,
)
    Ti_K = Ti_eV * EV_TO_J / KB_SI
    Te_K = Te_eV * EV_TO_J / KB_SI
    ions = Species(:proton, ME_SI * mi_me, QE_SI, Float64(n0), Ti_K, anisotropy_i * Ti_K)
    elec = Species(:electron, ME_SI, -QE_SI, Float64(n0), Te_K, anisotropy_e * Te_K)
    return PlasmaParams{Float64}(Float64(B0), [ions, elec], C_SI, EPS0_SI, MU0_SI, KB_SI, Float64(γ))
end

# --- Species lookup helpers ---

function get_electron(p::PlasmaParams)
    idx = findfirst(s -> s.q < 0, p.species)
    idx === nothing && error("No negatively charged (electron) species found in PlasmaParams.")
    return p.species[idx]
end

function get_ion(p::PlasmaParams)
    idx = findfirst(s -> s.q > 0, p.species)
    idx === nothing && error("No positively charged (ion) species found in PlasmaParams.")
    return p.species[idx]
end

# --- Derived plasma quantities ---

"Plasma frequency ``\\omega_{ps} = \\sqrt{n_s q_s^2 / (\\varepsilon_0 m_s)}``."
plasma_frequency(s::Species, p::PlasmaParams) = sqrt(s.n * s.q^2 / (p.ε0 * s.m))

"Signed cyclotron frequency ``\\Omega_s = q_s B_0 / m_s`` (negative for electrons)."
cyclotron_frequency(s::Species, p::PlasmaParams) = s.q * p.B0 / s.m

"Total mass density ``\\rho_0 = \\sum_s n_s m_s``."
mass_density(p::PlasmaParams) = sum(s -> s.n * s.m, p.species)

"Total equilibrium thermal pressure ``p_0 = \\sum_s n_s k_B T_{\\parallel s}``."
thermal_pressure(p::PlasmaParams) = sum(s -> s.n * p.kB * s.T_par, p.species)

"Alfvén speed ``v_A = B_0 / \\sqrt{\\mu_0 \\rho_0}``."
alfven_speed(p::PlasmaParams) = p.B0 / sqrt(p.μ0 * mass_density(p))

"MHD adiabatic sound speed ``v_s = \\sqrt{\\gamma p_0 / \\rho_0}``."
sound_speed(p::PlasmaParams) = sqrt(p.γ * thermal_pressure(p) / mass_density(p))

"""
    thermal_speed(s::Species, p::PlasmaParams; perp=false)

Kinetic thermal speed ``v_{th,s} = \\sqrt{2 k_B T_s / m_s}`` (Stix / Gary convention).
Use [`thermal_speed_1d`](@ref) for the RMS speed ``\\sqrt{k_B T_s / m_s}``.
"""
function thermal_speed(s::Species, p::PlasmaParams; perp::Bool = false)
    T = perp ? s.T_perp : s.T_par
    return sqrt(2 * p.kB * T / s.m)
end

"1D RMS thermal speed ``\\sqrt{k_B T_s / m_s}``."
function thermal_speed_1d(s::Species, p::PlasmaParams; perp::Bool = false)
    T = perp ? s.T_perp : s.T_par
    return sqrt(p.kB * T / s.m)
end

"Inertial length (skin depth) ``d_s = c / \\omega_{ps}``."
inertial_length(s::Species, p::PlasmaParams) = p.c / plasma_frequency(s, p)

"Thermal Larmor radius ``r_{Ls} = \\sqrt{k_B T_{\\perp s} / m_s} / |\\Omega_s|``."
larmor_radius(s::Species, p::PlasmaParams) = thermal_speed_1d(s, p; perp = true) / abs(cyclotron_frequency(s, p))

"Debye length ``\\lambda_{Ds} = \\sqrt{\\varepsilon_0 k_B T_{\\parallel s} / (n_s q_s^2)}``."
debye_length(s::Species, p::PlasmaParams) = sqrt(p.ε0 * p.kB * s.T_par / (s.n * s.q^2))

"Species parallel plasma beta ``\\beta_s = 2\\mu_0 n_s k_B T_{\\parallel s} / B_0^2``."
plasma_beta(s::Species, p::PlasmaParams; perp::Bool = false) =
    2 * p.μ0 * s.n * p.kB * (perp ? s.T_perp : s.T_par) / p.B0^2

"Total plasma beta ``\\beta = \\sum_s \\beta_s``."
plasma_beta(p::PlasmaParams) = sum(s -> plasma_beta(s, p), p.species)

# --- Wave eigenmode solution & diagnostics ---

"""
    WaveSolution{T<:AbstractFloat}

Linear (or exact circularly polarized) plane-wave eigenmode solution
``\\delta A(\\mathbf{r}, t) = \\Re[\\tilde{A} \\exp(i(\\mathbf{k}\\cdot\\mathbf{r} - \\omega t))]``
in a uniform plasma with ``\\mathbf{B}_0 = B_0 \\hat{\\mathbf{z}}``.
"""
struct WaveSolution{T<:AbstractFloat}
    mode::Symbol
    ω::Complex{T}
    k::SVector{3,T}
    δE::SVector{3,Complex{T}}
    δB::SVector{3,Complex{T}}
    δv::SVector{3,Complex{T}}
    δn::Complex{T}
    δp::Complex{T}
    plasma::PlasmaParams{T}
end

"Real phase speed ``v_{\\mathrm{ph}} = \\Re(\\omega) / |\\mathbf{k}|``."
phase_speed(sol::WaveSolution) = real(sol.ω) / norm(sol.k)

"Temporal growth/damping rate ``\\gamma = \\Im(\\omega)`` (negative = damped)."
damping_rate(sol::WaveSolution) = imag(sol.ω)

"Propagation angle ``\\theta`` between ``\\mathbf{k}`` and ``\\mathbf{B}_0 = B_0 \\hat{\\mathbf{z}}`` in radians."
propagation_angle(sol::WaveSolution) = atan(hypot(sol.k[1], sol.k[2]), sol.k[3])

"""
    polarization_E(sol::WaveSolution)

Electric-field polarization in the plane perpendicular to ``\\mathbf{B}_0``:
``P_E = i \\delta E_x / \\delta E_y``.
- ``\\Re(P_E) > 0``: right-hand polarized (electron gyration sense).
- ``\\Re(P_E) < 0``: left-hand polarized (ion gyration sense).
- ``P_E = +1`` / ``-1``: right / left circular polarization.
"""
function polarization_E(sol::WaveSolution)
    abs(sol.δE[2]) < 1e-14 * norm(sol.δE) && return ComplexF64(NaN)
    return im * sol.δE[1] / sol.δE[2]
end

"""
    polarization_B(sol::WaveSolution)

Magnetic-field polarization (Gary 1993) ``P_B = i b_s / b_A``, where
``\\hat{\\mathbf{A}} = \\hat{\\mathbf{k}}\\times\\hat{\\mathbf{z}}`` (out of ``(\\mathbf{k}, \\mathbf{B}_0)`` plane)
and ``\\hat{\\mathbf{S}} = \\hat{\\mathbf{A}}\\times\\hat{\\mathbf{k}}`` (in-plane, ``\\perp \\mathbf{k}``).
"""
function polarization_B(sol::WaveSolution)
    bnorm = norm(sol.δB)
    bnorm == 0 && return ComplexF64(NaN)
    khat = sol.k / norm(sol.k)
    zhat = SVector(0.0, 0.0, 1.0)
    avec = cross(khat, zhat)
    if norm(avec) < 1e-12
        # Parallel propagation: fall back to (x, y) basis
        bA = sol.δB[2]
        bs = sol.δB[1]
    else
        Ahat = avec / norm(avec)
        Shat = cross(Ahat, khat)
        bA = dot(Ahat, sol.δB)
        bs = dot(Shat, sol.δB)
    end
    abs(bA) < 1e-14 * bnorm && return ComplexF64(NaN)
    return im * bs / bA
end

"Magnetic compressibility ``C_B = |\\delta B_z|^2 / |\\delta\\mathbf{B}|^2``."
function magnetic_compressibility(sol::WaveSolution)
    b2 = sum(abs2, sol.δB)
    b2 == 0 && return 0.0
    return abs2(sol.δB[3]) / b2
end

"Alfvén ratio ``R_A = \\mu_0 \\rho_0 |\\delta\\mathbf{v}|^2 / |\\delta\\mathbf{B}|^2``."
function alfven_ratio(sol::WaveSolution)
    b2 = sum(abs2, sol.δB)
    b2 == 0 && return Inf
    v2 = sum(abs2, sol.δv)
    return sol.plasma.μ0 * mass_density(sol.plasma) * v2 / b2
end

"Electrostatic fraction ``|\\hat{\\mathbf{k}}\\cdot\\delta\\mathbf{E}| / |\\delta\\mathbf{E}| \\in [0, 1]``."
function electrostatic_ratio(sol::WaveSolution)
    enorm = norm(sol.δE)
    enorm == 0 && return 0.0
    khat = sol.k / norm(sol.k)
    return abs(dot(khat, sol.δE)) / enorm
end

"Time-averaged Poynting flux ``\\langle\\mathbf{S}\\rangle = \\frac{1}{2\\mu_0}\\Re(\\delta\\mathbf{E} \\times \\delta\\mathbf{B}^*)``."
function poynting_flux(sol::WaveSolution)
    return real.(cross(sol.δE, conj.(sol.δB))) / (2 * sol.plasma.μ0)
end

"""
    evaluate_fields(sol::WaveSolution, r, t; include_damping=true)

Evaluate real-valued wave perturbations `(E, B, v, n, p)` at position `r` (3-vector) and time `t`.
Set `include_damping=false` to visualize the undamped eigenmode structure with ``\\Re(\\omega)`` only.
"""
function evaluate_fields(sol::WaveSolution, r, t::Real; include_damping::Bool = true)
    ω_eff = include_damping ? sol.ω : Complex(real(sol.ω), 0.0)
    phase = exp(im * (sol.k[1] * r[1] + sol.k[2] * r[2] + sol.k[3] * r[3] - ω_eff * t))
    return (
        E = real.(sol.δE .* phase),
        B = real.(sol.δB .* phase),
        v = real.(sol.δv .* phase),
        n = real(sol.δn * phase),
        p = real(sol.δp * phase),
    )
end

# --- Interoperability with PlasmaPopulations.jl & PlasmaDispersionBase.jl ---

"""
    to_bo_plasma(p::PlasmaParams; Ω_ref=nothing)

Convert a `PlasmaParams` into a `(norm_plasma::NormalizedPlasma, scales)` pair suitable
for `PlasmaDispersionBase.DispersionProblem` and `PlasmaBO.solve(..., BOHH())`.

In `PlasmaDispersionBase`, `NormalizedPlasma` uses frequency unit `Ω_ref` and wavevector
unit `Ω_ref / c`, so `k_bo = k * p.c / Ω_ref` and `ω_phys = ω_bo * Ω_ref`.
"""
function to_bo_plasma(p::PlasmaParams; Ω_ref::Union{Nothing,Real} = nothing)
    Ωref = Ω_ref === nothing ? abs(cyclotron_frequency(get_ion(p), p)) : float(Ω_ref)
    Ωref == 0 && (Ωref = plasma_frequency(get_electron(p), p))

    pops = Tuple(
        map(p.species) do s
            wc = cyclotron_frequency(s, p) / Ωref
            wp2 = (plasma_frequency(s, p) / Ωref)^2
            vdf = if s.T_par <= 0 && s.T_perp <= 0
                PlasmaPopulations.Cold()
            else
                vpar = thermal_speed(s, p; perp = false)
                vperp = thermal_speed(s, p; perp = true)
                PlasmaPopulations.BiMaxwellian(vperp, vpar, 0.0)
            end
            PlasmaDispersionBase.NormalizedPopulation(wc, wp2, vdf, p.c)
        end
    )
    norm_plasma = PlasmaDispersionBase.NormalizedPlasma(pops)
    return norm_plasma, (; omega = Ωref, k = Ωref / p.c)
end

"""
    PlasmaParams(p::PlasmaPopulations.Plasma; γ=5/3)

Convert a physical SI `PlasmaPopulations.Plasma` into a `PlasmaParams` for analytical
theory calculations, transport ratios, and animations.
"""
function PlasmaParams(p::PlasmaPopulations.Plasma; γ::Real = 5 / 3)
    B0 = Float64(norm(p.B0))
    species = Species{Float64}[]
    for pop in p.populations
        part = pop.particle
        m = Float64(part.m)
        q = Float64(part.q)
        n = Float64(pop.n)
        vperp, vpar = PlasmaPopulations.vth_perp_para(pop.vdf)
        T_perp = 0.5 * m * Float64(vperp)^2 / KB_SI
        T_par  = 0.5 * m * Float64(vpar)^2 / KB_SI
        sname = q < 0 ? :electron : :ion
        push!(species, Species{Float64}(sname, m, q, n, T_par, T_perp))
    end
    return PlasmaParams{Float64}(B0, species, C_SI, EPS0_SI, MU0_SI, KB_SI, Float64(γ))
end

