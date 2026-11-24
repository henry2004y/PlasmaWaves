"""
Numerical dispersion confirmation powered by PlasmaBO.jl and PlasmaDispersionBase.jl.

Compares analytical (or symbolically derived) wave dispersion relations against full
hot magnetized kinetic electromagnetic solutions from `PlasmaBO.solve(..., BOHH())`,
using `PlasmaDispersionBase.track` for smooth branch continuation and quantifying
regime validity boundaries and kinetic damping/growth rates.
"""

"""
    DispersionComparison

Result of confirming a theoretical dispersion relation against `PlasmaBO.jl`.

Fields:
- `mode::Symbol`: wave mode label (e.g., `:alfven`, `:whistler`, `:kaw`, `:langmuir`)
- `ks::Vector{Float64}`: wavevector magnitudes (or `k_perp` for KAW sweeps)
- `θ::Float64`: propagation angle relative to ``\\mathbf{B}_0`` [rad]
- `k_label::String`: axis label for the swept wavevector coordinate (supports LaTeX `\$...\$`)
- `ω_theory::Vector{ComplexF64}`: theoretical complex frequency ``\\omega_{\\text{theory}}(k)``
- `ω_numerical::Vector{ComplexF64}`: primary matched/tracked kinetic root from `PlasmaBO`
- `rel_error::Vector{Float64}`: relative real-frequency discrepancy ``|\\Re\\omega_{\\text{th}} - \\Re\\omega_{\\text{num}}| / |\\Re\\omega_{\\text{num}}|``
- `damping_ratio::Vector{Float64}`: kinetic damping/growth ratio ``|\\Im\\omega_{\\text{num}}| / |\\Re\\omega_{\\text{num}}|``
- `ω_scale::Float64`: normalization frequency (e.g. ``\\Omega_i`` or ``\\omega_{pe}``) for plotting
- `ω_label::String`: frequency axis label (supports LaTeX `\$...\$`)
- `primary_label::String`: legend label for the primary `PlasmaBO` branch
- `extra_branches`: optional additional tracked kinetic branches `(; label, ω)` (e.g. Whistler R-branch alongside EMIC L-branch)
"""
struct DispersionComparison
    mode::Symbol
    ks::Vector{Float64}
    θ::Float64
    k_label::String
    ω_theory::Vector{ComplexF64}
    ω_numerical::Vector{ComplexF64}
    rel_error::Vector{Float64}
    damping_ratio::Vector{Float64}
    ω_scale::Float64
    ω_label::String
    primary_label::String
    extra_branches::Vector{NamedTuple{(:label, :ω), Tuple{String, Vector{ComplexF64}}}}
end

function DispersionComparison(
    mode::Symbol,
    ks::Vector{Float64},
    θ::Float64,
    k_label::String,
    ω_theory::Vector{ComplexF64},
    ω_numerical::Vector{ComplexF64},
    rel_error::Vector{Float64},
    damping_ratio::Vector{Float64},
    ω_scale::Float64,
    ω_label::String,
)
    return DispersionComparison(
        mode, ks, θ, k_label, ω_theory, ω_numerical, rel_error, damping_ratio,
        ω_scale, ω_label, "PlasmaBO (BOHH)",
        NamedTuple{(:label, :ω), Tuple{String, Vector{ComplexF64}}}[],
    )
end

"Return a boolean mask of `k` points where theoretical and kinetic real frequencies agree within `tol`."
validity_mask(comp::DispersionComparison; tol::Real = 0.05) = comp.rel_error .<= tol

function _eval_errors(ω_th::Vector{ComplexF64}, ω_num::Vector{ComplexF64})
    nk = length(ω_th)
    rel_err = Vector{Float64}(undef, nk)
    damp_rat = Vector{Float64}(undef, nk)
    for i in 1:nk
        wn = ω_num[i]
        if !isfinite(wn)
            rel_err[i] = NaN
            damp_rat[i] = NaN
        else
            denom = max(abs(real(wn)), 1e-12)
            rel_err[i] = abs(real(ω_th[i]) - real(wn)) / denom
            damp_rat[i] = abs(imag(wn)) / denom
        end
    end
    return rel_err, damp_rat
end

"""
    confirm_dispersion(theory_ω_fn, plasma, ks; θ=0.0, alg=BOHH(N=4, J=8),
                       mode=:custom, k_label="k", ω_scale=1.0, ω_label="ω",
                       k_vectors=nothing, Ω_ref=nothing, region_box=nothing,
                       track_seed=nothing, max_damping=0.6,
                       primary_label="PlasmaBO (BOHH)", extra_seeds=())

Confirm a theoretical dispersion function `theory_ω_fn(k) -> Complex` against the full
hot magnetized kinetic solver `PlasmaBO.solve(..., alg)`.

- `region_box`: optional `(lowerleft, upperright)` tuple in `Ω_ref`-normalized units to exclude
  spurious heavily damped Hermite poles when building the `PlasmaDispersionBase` survey.
- `track_seed`: optional `(ω_seed_phys, seed_idx)` tuple to trace a single continuous branch
  across `ks` using `PlasmaDispersionBase.track`.
- `extra_seeds`: optional collection of `(label, ω_seed_phys, seed_idx)` tuples to track
  additional physical branches from the same `SurveySolution`.
"""
function confirm_dispersion(
    theory_ω_fn,
    plasma::PlasmaParams,
    ks::AbstractVector{<:Real};
    θ::Real = 0.0,
    alg = PlasmaBO.BOHH(N = 4, J = 8),
    mode::Symbol = :custom,
    k_label::AbstractString = "\$k\$",
    ω_scale::Real = 1.0,
    ω_label::AbstractString = "\$\\omega\$",
    k_vectors::Union{Nothing,AbstractVector{<:PlasmaDispersionBase.Wavevector}} = nothing,
    Ω_ref::Union{Nothing,Real} = nothing,
    region_box::Union{Nothing,Tuple} = nothing,
    track_seed::Union{Nothing,Tuple{Number,Int}} = nothing,
    max_damping::Real = 0.6,
    primary_label::AbstractString = "PlasmaBO (BOHH)",
    extra_seeds = (),
)
    norm_p, u = to_bo_plasma(plasma; Ω_ref)
    ks_f = Float64.(collect(ks))
    nk = length(ks_f)

    ω_th = ComplexF64[complex(float(theory_ω_fn(k))) for k in ks_f]

    wvs = if k_vectors === nothing
        [PlasmaDispersionBase.Wavevector(; k = k / u.k, theta = float(θ)) for k in ks_f]
    else
        [PlasmaDispersionBase.Wavevector(w.k_perp / u.k, w.k_para / u.k) for w in k_vectors]
    end

    target = region_box === nothing ? PlasmaDispersionBase.Region(nothing) : PlasmaDispersionBase.Region(region_box)
    prob = PlasmaDispersionBase.DispersionProblem(norm_p, target, wvs)
    survey = PlasmaDispersionBase.solve(prob, alg)

    ω_num = if track_seed !== nothing
        w_seed, idx = track_seed
        br = PlasmaDispersionBase.track(survey, ComplexF64(w_seed) / u.omega, wvs[idx])
        ComplexF64.(br.omega .* u.omega)
    else
        out = Vector{ComplexF64}(undef, nk)
        for i in 1:nk
            roots_phys = survey.omega[i] .* u.omega
            cand = filter(roots_phys) do w
                isfinite(w) && real(w) > 0 && imag(w) >= -max_damping * max(real(w), 1e-6)
            end
            if isempty(cand)
                cand = filter(w -> isfinite(w) && real(w) > 0, roots_phys)
            end
            if isempty(cand)
                out[i] = ComplexF64(NaN, NaN)
            else
                wt = ω_th[i]
                out[i] = argmin(w -> abs(w - wt), cand)
            end
        end
        out
    end

    rel_err, damp_rat = _eval_errors(ω_th, ω_num)

    extras = NamedTuple{(:label, :ω), Tuple{String, Vector{ComplexF64}}}[]
    for (lbl, w_seed, idx) in extra_seeds
        br = PlasmaDispersionBase.track(survey, ComplexF64(w_seed) / u.omega, wvs[idx])
        push!(extras, (label = String(lbl), ω = ComplexF64.(br.omega .* u.omega)))
    end

    return DispersionComparison(
        mode,
        ks_f,
        Float64(θ),
        String(k_label),
        ω_th,
        ω_num,
        rel_err,
        damp_rat,
        Float64(ω_scale),
        String(ω_label),
        String(primary_label),
        extras,
    )
end

# Also accept a PlasmaPopulations.Plasma directly!
function confirm_dispersion(theory_ω_fn, plasma::PlasmaPopulations.Plasma, ks::AbstractVector{<:Real}; kwargs...)
    return confirm_dispersion(theory_ω_fn, PlasmaParams(plasma), ks; kwargs...)
end

# --- Built-in Confirmation Presets ---

"""
    confirm_alfven(plasma, ks; θ=0.0, alg=BOHH(N=3, J=6))

Confirm the Ideal MHD Shear Alfvén dispersion ``\\omega = k_\\parallel v_A`` against `PlasmaBO`.
Tracks both the left-hand EMIC (L) and right-hand Whistler (R) kinetic branches backward
from high ``k d_i`` where they are cleanly separated, revealing Hall splitting as ``k d_i \\sim 1``.
"""
function confirm_alfven(
    plasma::PlasmaParams = normalized_plasma(; mi_me = 100.0, β_i = 0.05, β_e = 0.05, vA_c = 0.05),
    ks::AbstractVector{<:Real} = range(0.05, 1.2, length = 30);
    θ::Real = 0.0,
    alg = PlasmaBO.BOHH(N = 3, J = 6),
)
    Ωi = abs(cyclotron_frequency(get_ion(plasma), plasma))
    k_last = last(ks)
    w_emic_end = hall_mhd_dispersion(plasma, k_last * cos(θ); branch = :emic)
    w_whis_end = whistler_wave(plasma, k_last, θ; model = :cold).ω

    box = (0.01 - 0.15im, max(2.5, 1.3 * real(w_whis_end) / Ωi) + 0.05im)
    return confirm_dispersion(
        k -> shear_alfven_wave(plasma, k, θ).ω,
        plasma,
        ks;
        θ,
        alg,
        mode = :alfven,
        k_label = "\$k d_i\$",
        ω_scale = Ωi,
        ω_label = "\$\\omega / \\Omega_i\$",
        region_box = box,
        track_seed = (w_emic_end, length(ks)),
        primary_label = "PlasmaBO EMIC (L)",
        extra_seeds = [("PlasmaBO Whistler (R)", w_whis_end, length(ks))],
    )
end

"""
    confirm_whistler(plasma, ks; θ=0.0, model=:cold, alg=BOHH(N=4, J=8))

Confirm the Cold/Hall Whistler wave dispersion against `PlasmaBO`.
Demonstrates agreement in the fluid whistler band and onset of electron cyclotron
damping as ``\\omega \\to |\\Omega_e|``.
"""
function confirm_whistler(
    plasma::PlasmaParams = normalized_plasma(; mi_me = 100.0, β_i = 0.05, β_e = 0.05, vA_c = 0.05),
    ks::AbstractVector{<:Real} = range(0.2, 4.0, length = 24);
    θ::Real = 0.0,
    model::Symbol = :cold,
    alg = PlasmaBO.BOHH(N = 4, J = 8),
)
    Ωi = abs(cyclotron_frequency(get_ion(plasma), plasma))
    return confirm_dispersion(
        k -> whistler_wave(plasma, k, θ; model).ω,
        plasma,
        ks;
        θ,
        alg,
        mode = :whistler,
        k_label = "\$k d_i\$",
        ω_scale = Ωi,
        ω_label = "\$\\omega / \\Omega_i\$",
        primary_label = "PlasmaBO Whistler",
    )
end

"""
    confirm_kaw(plasma, k_par, k_perps; model=:bessel, alg=BOHH(N=5, J=6))

Confirm the Kinetic Alfvén Wave (KAW) dispersion against `PlasmaBO` across perpendicular
wavenumbers `k_perps` at fixed `k_par`, revealing kinetic Landau damping and finite-β corrections.
"""
function confirm_kaw(
    plasma::PlasmaParams = normalized_plasma(; mi_me = 100.0, β_i = 0.2, β_e = 0.2, vA_c = 0.05),
    k_par::Real = 0.05,
    k_perps::AbstractVector{<:Real} = range(0.2, 4.0, length = 16);
    model::Symbol = :bessel,
    alg = PlasmaBO.BOHH(N = 5, J = 6),
)
    vA = alfven_speed(plasma)
    riL = larmor_radius(get_ion(plasma), plasma)
    wvs = [PlasmaDispersionBase.Wavevector(Float64(kp), Float64(k_par)) for kp in k_perps]
    return confirm_dispersion(
        kp -> kinetic_alfven_wave(plasma, k_par, kp; model).ω,
        plasma,
        k_perps .* riL;
        θ = atan(last(k_perps), k_par),
        alg,
        mode = :kaw,
        k_label = "\$k_\\perp r_{iL}\$",
        ω_scale = k_par * vA,
        ω_label = "\$\\omega / (k_\\parallel v_A)\$",
        k_vectors = wvs,
        primary_label = "PlasmaBO KAW",
    )
end

"""
    confirm_langmuir(plasma, kλDs; alg=BOHH(N=2, J=8))

Confirm the warm-fluid Bohm-Gross dispersion ``\\omega^2 = \\omega_{pe}^2 + 3 k^2 v_{th,e}^2``
against full Vlasov kinetic `PlasmaBO` across ``k \\lambda_{De}``, revealing the onset of
collisionless Landau damping.
"""
function confirm_langmuir(
    plasma::PlasmaParams = normalized_plasma(; mi_me = 100.0, β_i = 0.05, β_e = 0.5, vA_c = 0.05),
    kλDs::AbstractVector{<:Real} = range(0.08, 0.45, length = 16);
    alg = PlasmaBO.BOHH(N = 2, J = 8),
)
    e = get_electron(plasma)
    ωpe = plasma_frequency(e, plasma)
    λDe = debye_length(e, plasma)
    wvs = [PlasmaDispersionBase.Wavevector(0.0, Float64(kλ / λDe)) for kλ in kλDs]
    return confirm_dispersion(
        kλ -> langmuir_wave(plasma, kλ / λDe; model = :bohm_gross).ω,
        plasma,
        kλDs;
        θ = 0.0,
        alg,
        mode = :langmuir,
        k_label = "\$k \\lambda_{De}\$",
        ω_scale = ωpe,
        ω_label = "\$\\omega / \\omega_{pe}\$",
        k_vectors = wvs,
        Ω_ref = ωpe,
        primary_label = "PlasmaBO Langmuir",
    )
end

"""
    confirm_ion_acoustic(plasma, ks; alg=BOHH(N=2, J=8))

Confirm the fluid Ion-Acoustic dispersion against full Vlasov kinetic `PlasmaBO`.
"""
function confirm_ion_acoustic(
    plasma::PlasmaParams = normalized_plasma(; mi_me = 100.0, β_i = 0.02, β_e = 0.30, vA_c = 0.05),
    ks::AbstractVector{<:Real} = range(0.2, 2.0, length = 16);
    alg = PlasmaBO.BOHH(N = 2, J = 8),
)
    Ωi = abs(cyclotron_frequency(get_ion(plasma), plasma))
    return confirm_dispersion(
        k -> ion_acoustic_wave(plasma, k; model = :fluid).ω,
        plasma,
        ks;
        θ = 0.0,
        alg,
        mode = :ion_acoustic,
        k_label = "\$k d_i\$",
        ω_scale = Ωi,
        ω_label = "\$\\omega / \\Omega_i\$",
        primary_label = "PlasmaBO Ion-Acoustic",
    )
end
