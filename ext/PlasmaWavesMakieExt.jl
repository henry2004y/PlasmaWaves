module PlasmaWavesMakieExt

using Makie
using Makie: Figure, Axis, Axis3, PolarAxis, Observable, @lift, Point2f, Point3f
using Makie: lines!, scatter!, linesegments!, heatmap!, band!, hlines!, vlines!, text!, axislegend, Colorbar, Label, record
using Makie: @L_str, LaTeXString, to_colormap
using LinearAlgebra
using StaticArrays
using PlasmaWaves

# ColorBrewer Set2_8 discrete palette
const SET2 = to_colormap(:Set2_8)
const C_FAST     = SET2[2]
const C_ALFVEN   = SET2[3]
const C_SLOW     = SET2[1]
const C_WHISTLER = SET2[6]
const C_EMIC     = SET2[4]
const C_EXTRA    = SET2[5]
const C_EFIELD   = SET2[2]
const C_BFIELD   = SET2[3]
const C_VEL      = SET2[1]
const C_PARTICLE = SET2[6]

function _to_latex(s::AbstractString)
    return startswith(s, '$') ? LaTeXString(s) : LaTeXString("\$" * s * "\$")
end

function _finalize_anim(fig, t_obs, ts, path; fps::Int = 24)
    if path === nothing
        return (; fig, t = t_obs)
    end
    mkpath(dirname(abspath(path)))
    record(fig, path, ts; framerate = fps) do ti
        t_obs[] = ti
    end
    return path
end

# ==============================================================================
# 1. Static Diagnostic & Dispersion Plots
# ==============================================================================

"""
    plot_friedrich(; vA=1.0, vs_low=0.5, vs_high=1.4)

Side-by-side ideal MHD Friedrich phase-velocity and group-velocity polar diagrams
for low-β (`v_s < v_A`) and high-β (`v_s > v_A`) regimes.
"""
function PlasmaWaves.plot_friedrich(; vA::Real = 1.0, vs_low::Real = 0.5, vs_high::Real = 1.4)
    θs = range(0, 2π, length = 361)
    fig = Figure(size = (1020, 490), fontsize = 18)

    for (col, (vs, title_tex)) in enumerate([
        (vs_low, LaTeXString("\\text{Low-}\\beta\\; (v_s / v_A = $(round(vs_low/vA, digits=2)))")),
        (vs_high, LaTeXString("\\text{High-}\\beta\\; (v_s / v_A = $(round(vs_high/vA, digits=2)))")),
    ])
        ax = Axis(
            fig[1, col];
            title = title_tex,
            titlesize = 19,
            xlabel = L"v_x / v_A\; (\perp \mathbf{B}_0)",
            ylabel = L"v_z / v_A\; (\parallel \mathbf{B}_0)",
            aspect = DataAspect(),
        )
        vf = [friedrich_phase_speeds(vA, vs, θ).v_fast for θ in θs]
        va = [friedrich_phase_speeds(vA, vs, θ).v_alfven for θ in θs]
        vslow = [friedrich_phase_speeds(vA, vs, θ).v_slow for θ in θs]

        lines!(ax, vf .* sin.(θs), vf .* cos.(θs); color = C_FAST, linewidth = 3.0, label = L"\text{Fast } (v_{\mathrm{ph}})")
        lines!(ax, va .* sin.(θs), va .* cos.(θs); color = C_ALFVEN, linewidth = 3.0, linestyle = :dash, label = L"\text{Alfvén } (v_{\mathrm{ph}})")
        lines!(ax, vslow .* sin.(θs), vslow .* cos.(θs); color = C_SLOW, linewidth = 3.0, label = L"\text{Slow } (v_{\mathrm{ph}})")

        vg_slow = [friedrich_group_velocities(vA, vs, θ).vg_slow for θ in θs]
        lines!(ax, first.(vg_slow), last.(vg_slow); color = C_SLOW, linewidth = 2.0, linestyle = :dot, label = L"\text{Slow } (v_g)")

        lim = 1.2 * max(vA, vs)
        lines!(ax, [0.0, 0.0], [-lim, lim]; color = (:black, 0.25), linestyle = :dash)
        lines!(ax, [-lim, lim], [0.0, 0.0]; color = (:black, 0.25), linestyle = :dash)
        col == 1 && axislegend(ax; position = :rb, framevisible = false, labelsize = 15)
    end
    return fig
end

"""
    plot_dispersion_parallel(plasma::PlasmaParams)

Log-log cold-plasma dispersion relation for parallel propagation (`θ = 0`):
R-mode (Whistler + free-space), L-mode (EMIC + free-space), and light line `ω = c k`.
"""
function PlasmaWaves.plot_dispersion_parallel(plasma::PlasmaParams = normalized_plasma(; mi_me = 25.0, vA_c = 0.2))
    cut = cold_cutoffs(plasma)
    res = cold_resonances(plasma)
    Ωi, Ωe = res.Ω_i, res.Ω_e
    ωL, ωR = cut.ω_L, cut.ω_R
    ωmax = 3.5 * max(ωR, Ωe)

    ω_emic = exp.(range(log(0.02 * Ωi), log(0.98 * Ωi), length = 200))
    ω_whis = exp.(range(log(0.02 * Ωi), log(0.98 * Ωe), length = 300))
    ω_Lfree = exp.(range(log(1.01 * ωL), log(ωmax), length = 200))
    ω_Rfree = exp.(range(log(1.01 * ωR), log(ωmax), length = 200))

    k_of_L(ωs) = [ω * sqrt(max(real(stix_parameters(plasma, ω).L), 0.0)) / plasma.c for ω in ωs]
    k_of_R(ωs) = [ω * sqrt(max(real(stix_parameters(plasma, ω).R), 0.0)) / plasma.c for ω in ωs]

    k_emic = k_of_L(ω_emic)
    k_whis = k_of_R(ω_whis)
    k_Lfree = k_of_L(ω_Lfree)
    k_Rfree = k_of_R(ω_Rfree)

    fig = Figure(size = (880, 540), fontsize = 18)
    ax = Axis(
        fig[1, 1];
        title = L"\text{Parallel Cold-Plasma Dispersion } (\mathbf{k} \parallel \mathbf{B}_0)",
        titlesize = 19,
        xlabel = L"k d_i",
        ylabel = L"\omega / \Omega_i",
        xscale = log10,
        yscale = log10,
    )
    di = inertial_length(get_ion(plasma), plasma)

    lines!(ax, k_whis .* di, ω_whis ./ Ωi; color = C_WHISTLER, linewidth = 3.2, label = L"\text{R-mode: Whistler / ECW}")
    lines!(ax, k_emic .* di, ω_emic ./ Ωi; color = C_EMIC, linewidth = 3.2, label = L"\text{L-mode: Alfvén / EMIC}")
    lines!(ax, k_Rfree .* di, ω_Rfree ./ Ωi; color = C_FAST, linewidth = 2.6, linestyle = :dash, label = L"\text{R free-space } (\omega > \omega_R)")
    lines!(ax, k_Lfree .* di, ω_Lfree ./ Ωi; color = C_ALFVEN, linewidth = 2.6, linestyle = :dash, label = L"\text{L free-space } (\omega > \omega_L)")

    k_ref = exp.(range(log(max(minimum(k_whis), 1e-3)), log(maximum(k_whis)), length = 100))
    lines!(ax, k_ref .* di, (plasma.c .* k_ref) ./ Ωi; color = (:black, 0.5), linestyle = :dot, linewidth = 2.0, label = L"\text{Light } (\omega = c k)")
    lines!(ax, k_ref .* di, (alfven_speed(plasma) .* k_ref) ./ Ωi; color = (:gray, 0.7), linestyle = :dashdot, linewidth = 2.0, label = L"\text{Alfvén } (\omega = k v_A)")

    hlines!(ax, [1.0, Ωe / Ωi]; color = (:black, 0.35), linestyle = :dot)
    axislegend(ax; position = :lt, framevisible = false, labelsize = 15)
    return fig
end

"""
    plot_cma_diagram(; mi_me=10.0)

Clemmow–Mullaly–Allis (CMA) parameter-space diagram `(X = ω_pe²/ω², Y = Ω_e/ω)`
showing the cutoff boundaries (`P=0`, `R=0`, `L=0`) and resonance boundaries.
"""
function PlasmaWaves.plot_cma_diagram(; mi_me::Real = 10.0)
    μ = float(mi_me)
    fig = Figure(size = (820, 600), fontsize = 18)
    ax = Axis(
        fig[1, 1];
        title = LaTeXString("\\text{CMA Diagram } (m_i/m_e = $(Int(round(μ))))"),
        titlesize = 19,
        xlabel = L"X = \omega_{pe}^2 / \omega^2\; (\propto n_0)",
        ylabel = L"Y = |\Omega_e| / \omega\; (\propto B_0)",
        limits = (0, 2.5, 0, 2.5),
    )

    Y = range(0.001, 2.5, length = 400)
    vlines!(ax, [1.0 / (1 + 1 / μ)]; color = :black, linewidth = 2.6, label = L"P = 0\; (\omega = \omega_p)")
    hlines!(ax, [1.0]; color = C_WHISTLER, linewidth = 2.6, linestyle = :dash, label = L"R \to \infty\; (\omega = \Omega_e)")
    X_R = @. (1 - Y) * (1 + Y / μ) / (1 + 1 / μ)
    mask_R = (X_R .>= 0) .& (X_R .<= 2.5)
    lines!(ax, X_R[mask_R], Y[mask_R]; color = C_FAST, linewidth = 2.8, label = L"R = 0\; (\omega_R)")
    X_L = @. (1 + Y) * (1 - Y / μ) / (1 + 1 / μ)
    mask_L = (X_L .>= 0) .& (X_L .<= 2.5)
    lines!(ax, X_L[mask_L], Y[mask_L]; color = C_ALFVEN, linewidth = 2.8, label = L"L = 0\; (\omega_L)")
    Y_uh = range(0.0, 0.999, length = 200)
    X_uh = @. (1 - Y_uh^2) * (1 - (Y_uh / μ)^2) / (1 + 1 / μ - (Y_uh^2 / μ) * (1 + 1 / μ))
    lines!(ax, X_uh, Y_uh; color = C_SLOW, linewidth = 2.8, linestyle = :dot, label = L"S = 0\; (\omega_{UH})")

    axislegend(ax; position = :rt, framevisible = true, labelsize = 15)
    return fig
end

"""
    plot_kaw_iaw_dispersion(plasma::PlasmaParams)

Two-panel comparison of Shear Alfvén Wave (SAW), Kinetic Alfvén Wave (KAW, Bessel & two-fluid),
and Inertial Alfvén Wave (IAW): phase speed `ω / (k_∥ v_A)` and field ratios `|δE_⊥/δB_⊥|` and `|δE_∥/δE_⊥|`.
"""
function PlasmaWaves.plot_kaw_iaw_dispersion(plasma::PlasmaParams = normalized_plasma(; β_i = 0.2, β_e = 0.2))
    ion = get_ion(plasma)
    riL = larmor_radius(ion, plasma)
    vA = alfven_speed(plasma)

    plasma_cold = normalized_plasma(; β_i = 1e-5, β_e = 1e-5)
    de_cold = inertial_length(get_electron(plasma_cold), plasma_cold)

    xs = range(0.01, 4.0, length = 200)
    k_par = 0.05

    vph_kaw_bessel = [real(kinetic_alfven_wave(plasma, k_par, x / riL; model = :bessel).ω) / (k_par * vA) for x in xs]
    vph_kaw_2f = [real(kinetic_alfven_wave(plasma, k_par, x / riL; model = :twofluid).ω) / (k_par * vA) for x in xs]
    vph_iaw = [real(inertial_alfven_wave(plasma_cold, k_par, x / de_cold).ω) / (k_par * vA) for x in xs]

    ez_ex = [
        let s = kinetic_alfven_wave(plasma, k_par, x / riL; model = :bessel)
            abs(s.δE[3]) / abs(s.δE[1])
        end for x in xs
    ]
    ex_by = [
        let s = kinetic_alfven_wave(plasma, k_par, x / riL; model = :bessel)
            abs(s.δE[1]) / (abs(s.δB[2]) * vA)
        end for x in xs
    ]

    fig = Figure(size = (1040, 460), fontsize = 18)
    ax1 = Axis(
        fig[1, 1];
        title = L"\text{Dispersive Alfvén Phase Speed}",
        titlesize = 19,
        xlabel = L"k_\perp r_{iL}\; (\text{or } k_\perp d_e \text{ for IAW})",
        ylabel = L"\omega / (k_\parallel v_A)",
    )
    lines!(ax1, xs, ones(length(xs)); color = (:black, 0.55), linestyle = :dash, linewidth = 2.0, label = L"\text{Ideal SAW}")
    lines!(ax1, xs, vph_kaw_bessel; color = C_FAST, linewidth = 3.0, label = L"\text{KAW (Bessel)}")
    lines!(ax1, xs, vph_kaw_2f; color = C_WHISTLER, linewidth = 2.5, linestyle = :dot, label = L"\text{KAW (Two-Fluid)}")
    lines!(ax1, xs, vph_iaw; color = C_ALFVEN, linewidth = 3.0, label = L"\text{IAW } (\beta_e \ll m_e/m_i)")
    axislegend(ax1; position = :lt, framevisible = false, labelsize = 15)

    ax2 = Axis(
        fig[1, 2];
        title = L"\text{KAW Polarization \& } E_\parallel",
        titlesize = 19,
        xlabel = L"k_\perp r_{iL}",
        ylabel = L"\text{Normalized Ratio}",
    )
    lines!(ax2, xs, ex_by; color = C_BFIELD, linewidth = 3.0, label = L"|\delta E_\perp / (v_A \delta B_\perp)|")
    lines!(ax2, xs, ez_ex .* 20; color = C_EFIELD, linewidth = 3.0, label = L"20\,|\delta E_\parallel / \delta E_\perp|")
    axislegend(ax2; position = :lt, framevisible = false, labelsize = 15)
    return fig
end

"""
    plot_landau_dispersion(plasma::PlasmaParams)

Two-panel Vlasov kinetic electrostatic dispersion:
Left: Langmuir wave real frequency `ω_r / ω_pe` and Landau damping `-γ / ω_pe`.
Right: Ion-Acoustic wave damping ratio `|γ / ω_r|` vs `T_e / T_i`.
"""
function PlasmaWaves.plot_landau_dispersion(plasma::PlasmaParams = normalized_plasma(; β_i = 0.05, β_e = 0.5))
    e = get_electron(plasma)
    ωpe = plasma_frequency(e, plasma)
    λDe = debye_length(e, plasma)

    kλs = range(0.08, 0.8, length = 80)
    ωr_v = Float64[]
    γ_v = Float64[]
    ωr_bg = Float64[]
    for kλ in kλs
        k = kλ / λDe
        s_v = langmuir_wave(plasma, k; model = :vlasov)
        s_bg = langmuir_wave(plasma, k; model = :bohm_gross)
        push!(ωr_v, real(s_v.ω) / ωpe)
        push!(γ_v, -imag(s_v.ω) / ωpe)
        push!(ωr_bg, real(s_bg.ω) / ωpe)
    end

    τs = exp.(range(log(1.0), log(30.0), length = 50))
    ia_damp = Float64[]
    for τ in τs
        p_ia = normalized_plasma(; β_i = 0.02, β_e = 0.02 * τ)
        λD = debye_length(get_electron(p_ia), p_ia)
        s_ia = ion_acoustic_wave(p_ia, 0.1 / λD; model = :vlasov)
        push!(ia_damp, abs(imag(s_ia.ω)) / max(real(s_ia.ω), 1e-12))
    end

    fig = Figure(size = (1040, 460), fontsize = 18)
    ax1 = Axis(
        fig[1, 1];
        title = L"\text{Langmuir Wave \& Landau Damping}",
        titlesize = 19,
        xlabel = L"k \lambda_{De}",
        ylabel = L"\text{Frequency} / \omega_{pe}",
    )
    lines!(ax1, kλs, ωr_v; color = C_ALFVEN, linewidth = 3.0, label = L"\omega_r / \omega_{pe}\; (\text{Vlasov})")
    lines!(ax1, kλs, ωr_bg; color = :black, linewidth = 2.2, linestyle = :dash, label = L"\omega_r / \omega_{pe}\; (\text{Bohm-Gross})")
    lines!(ax1, kλs, γ_v; color = C_FAST, linewidth = 3.0, label = L"-\gamma / \omega_{pe}\; (\text{Landau})")
    axislegend(ax1; position = :lt, framevisible = false, labelsize = 15)

    ax2 = Axis(
        fig[1, 2];
        title = L"\text{Ion-Acoustic Damping vs. } T_e / T_i",
        titlesize = 19,
        xlabel = L"T_e / T_i",
        ylabel = L"|\gamma / \omega_r|",
        xscale = log10,
    )
    lines!(ax2, τs, ia_damp; color = C_SLOW, linewidth = 3.2, label = L"\text{Ion-Acoustic } |\gamma/\omega_r|")
    hlines!(ax2, [0.1]; color = (:black, 0.4), linestyle = :dot, linewidth = 2.0, label = L"\text{Weak damping } (0.1)")
    axislegend(ax2; position = :rt, framevisible = false, labelsize = 15)
    return fig
end

"""
    plot_bernstein_dispersion(plasma::PlasmaParams; species=:electron, n_harmonics=5)

Plot Electron or Ion Bernstein wave harmonic passbands `ω / Ω_s` vs `k_⊥ r_Ls`.
"""
function PlasmaWaves.plot_bernstein_dispersion(
    plasma::PlasmaParams = normalized_plasma(; mi_me = 1836.0, β_i = 0.05, β_e = 0.05, vA_c = 0.05);
    species::Symbol = :electron,
    n_harmonics::Int = 5,
)
    sp = species === :electron ? get_electron(plasma) : get_ion(plasma)
    Ωs = abs(cyclotron_frequency(sp, plasma))
    rL = larmor_radius(sp, plasma)

    xs = range(0.05, 5.0, length = 120)
    fig = Figure(size = (840, 520), fontsize = 18)
    sname = species === :electron ? "Electron" : "Ion"
    ax = Axis(
        fig[1, 1];
        title = LaTeXString("\\text{$sname Bernstein Waves } (k_\\parallel = 0)"),
        titlesize = 19,
        xlabel = L"k_\perp r_L",
        ylabel = L"\omega / \Omega",
    )

    for n in 1:n_harmonics
        hlines!(ax, [Float64(n)]; color = (:black, 0.25), linestyle = :dot)
        ωs = [real(bernstein_wave(plasma, x / rL; species, harmonic = n).ω) / Ωs for x in xs]
        lines!(ax, xs, ωs; color = SET2[mod1(n, length(SET2))], linewidth = 3.0, label = LaTeXString("n = $n"))
    end
    hlines!(ax, [Float64(n_harmonics + 1)]; color = (:black, 0.25), linestyle = :dot)
    axislegend(ax; position = :rt, framevisible = true, labelsize = 15)
    return fig
end

"""
    plot_transport_ratios(plasma::PlasmaParams)

Plot Gary's observational wave-identification ratios (`C_B` and `R_A`) across propagation angle `θ`.
"""
function PlasmaWaves.plot_transport_ratios(plasma::PlasmaParams = normalized_plasma(; β_i = 0.15, β_e = 0.15))
    θs = range(0.05, π / 2 - 0.05, length = 120)
    θdeg = rad2deg.(θs)

    cb_f = [magnetic_compressibility(fast_magnetosonic_wave(plasma, 1.0, θ)) for θ in θs]
    cb_s = [magnetic_compressibility(slow_magnetosonic_wave(plasma, 1.0, θ)) for θ in θs]
    cb_a = [magnetic_compressibility(shear_alfven_wave(plasma, 1.0, θ)) for θ in θs]

    ra_f = [alfven_ratio(fast_magnetosonic_wave(plasma, 1.0, θ)) for θ in θs]
    ra_s = [min(alfven_ratio(slow_magnetosonic_wave(plasma, 1.0, θ)), 10.0) for θ in θs]
    ra_a = [alfven_ratio(shear_alfven_wave(plasma, 1.0, θ)) for θ in θs]

    fig = Figure(size = (1040, 460), fontsize = 18)
    ax1 = Axis(
        fig[1, 1];
        title = L"\text{Magnetic Compressibility } C_B",
        titlesize = 19,
        xlabel = L"\theta\; [^\circ]",
        ylabel = L"C_B = |\delta B_\parallel|^2 / |\delta\mathbf{B}|^2",
    )
    lines!(ax1, θdeg, cb_f; color = C_FAST, linewidth = 3.0, label = L"\text{Fast}")
    lines!(ax1, θdeg, cb_a; color = C_ALFVEN, linewidth = 3.0, linestyle = :dash, label = L"\text{Alfvén}")
    lines!(ax1, θdeg, cb_s; color = C_SLOW, linewidth = 3.0, label = L"\text{Slow}")
    axislegend(ax1; position = :lt, framevisible = false, labelsize = 15)

    ax2 = Axis(
        fig[1, 2];
        title = L"\text{Alfvén Ratio } R_A",
        titlesize = 19,
        xlabel = L"\theta\; [^\circ]",
        ylabel = L"R_A = \mu_0 \rho_0 |\delta\mathbf{v}|^2 / |\delta\mathbf{B}|^2",
    )
    lines!(ax2, θdeg, ra_f; color = C_FAST, linewidth = 3.0, label = L"\text{Fast}")
    lines!(ax2, θdeg, ra_a; color = C_ALFVEN, linewidth = 3.0, linestyle = :dash, label = L"\text{Alfvén } (R_A = 1)")
    lines!(ax2, θdeg, ra_s; color = C_SLOW, linewidth = 3.0, label = L"\text{Slow (capped at 10)}")
    axislegend(ax2; position = :lt, framevisible = false, labelsize = 15)
    return fig
end

"""
    plot_confirmation(comp::DispersionComparison; tol=0.05)

Three-panel Theory-vs-Numerical Confirmation figure comparing an analytical (or symbolically
derived) wave dispersion relation against the full hot magnetized kinetic `PlasmaBO.jl` solver:
1. Real frequency `ω_r(k)` — Theory (dashed) vs. tracked `PlasmaBO` branches
2. Kinetic damping/growth rate `γ(k)` from `PlasmaBO`
3. Relative real-frequency discrepancy `|Δω_r / ω_r|` with validity threshold `tol`
"""
function PlasmaWaves.plot_confirmation(comp::DispersionComparison; tol::Real = 0.05)
    ks = comp.ks
    s = comp.ω_scale
    ωr_th = real.(comp.ω_theory) ./ s
    ωr_bo = real.(comp.ω_numerical) ./ s
    γ_bo  = imag.(comp.ω_numerical) ./ s
    err_pct = comp.rel_error .* 100

    k_tex = _to_latex(comp.k_label)
    ω_tex = _to_latex(comp.ω_label)
    ω_inner = strip(String(comp.ω_label), '$')
    γ_tex = LaTeXString("\\gamma / ($(ω_inner)\\text{ scale})")

    fig = Figure(size = (1220, 440), fontsize = 18)

    ax1 = Axis(
        fig[1, 1];
        title = LaTeXString("\\text{Theory vs. PlasmaBO: } \\mathrm{$(comp.mode)}"),
        titlesize = 19,
        xlabel = k_tex,
        ylabel = ω_tex,
    )
    lines!(ax1, ks, ωr_th; color = :black, linewidth = 2.8, linestyle = :dash, label = "Analytical Theory")
    lines!(ax1, ks, ωr_bo; color = C_ALFVEN, linewidth = 3.0, label = comp.primary_label)
    scatter!(ax1, ks, ωr_bo; color = C_ALFVEN, markersize = 7)

    ax2 = Axis(
        fig[1, 2];
        title = L"\text{Kinetic Damping / Growth}",
        titlesize = 19,
        xlabel = k_tex,
        ylabel = γ_tex,
    )
    hlines!(ax2, [0.0]; color = (:black, 0.35), linestyle = :dot)
    lines!(ax2, ks, γ_bo; color = C_ALFVEN, linewidth = 3.0, label = comp.primary_label)
    scatter!(ax2, ks, γ_bo; color = C_ALFVEN, markersize = 7)

    ax3 = Axis(
        fig[1, 3];
        title = L"\text{Theory Breakdown Diagnostic}",
        titlesize = 19,
        xlabel = k_tex,
        ylabel = L"|\Delta\omega_r / \omega_r|\; [\%]",
    )
    lines!(ax3, ks, err_pct; color = C_ALFVEN, linewidth = 3.0, label = comp.primary_label)

    # Render any additional tracked kinetic branches (e.g. Whistler alongside EMIC)
    extra_colors = [C_WHISTLER, C_EMIC, C_EXTRA]
    for (idx, br) in enumerate(comp.extra_branches)
        col = extra_colors[mod1(idx, length(extra_colors))]
        ωr_ex = real.(br.ω) ./ s
        γ_ex  = imag.(br.ω) ./ s
        err_ex = [isfinite(br.ω[i]) ? abs(real(comp.ω_theory[i]) - real(br.ω[i])) / max(abs(real(br.ω[i])), 1e-12) * 100 : NaN for i in eachindex(ks)]

        lines!(ax1, ks, ωr_ex; color = col, linewidth = 3.0, label = br.label)
        scatter!(ax1, ks, ωr_ex; color = col, markersize = 7)
        lines!(ax2, ks, γ_ex; color = col, linewidth = 3.0, label = br.label)
        scatter!(ax2, ks, γ_ex; color = col, markersize = 7)
        lines!(ax3, ks, err_ex; color = col, linewidth = 3.0, label = br.label)
    end

    hlines!(ax3, [tol * 100]; color = C_SLOW, linestyle = :dash, linewidth = 2.2, label = "$(round(tol * 100, digits=1))% threshold")
    axislegend(ax1; position = :lt, framevisible = false, labelsize = 14)
    axislegend(ax2; position = :lb, framevisible = false, labelsize = 14)
    axislegend(ax3; position = :lt, framevisible = false, labelsize = 14)

    return fig
end

# ==============================================================================
# 2. Animated Wave Propagation Demos (MP4 / GIF / WebM / Interactive)
# ==============================================================================

"""
    animate_alfven_wave(path=nothing; plasma=normalized_plasma(), k=1.0, nframes=72, fps=24)

3D animation of a Shear Alfvén Wave propagating along `B₀ = B₀ ẑ`.
"""
function PlasmaWaves.animate_alfven_wave(
    path::Union{Nothing,AbstractString} = nothing;
    plasma::PlasmaParams = normalized_plasma(),
    k::Real = 1.0,
    nframes::Int = 72,
    fps::Int = 24,
)
    sol = shear_alfven_wave(plasma, k, 0.0; amp = 0.35)
    ω = real(sol.ω)
    T = 2π / ω
    ts = range(0, T, length = nframes + 1)[1:end-1]
    zs = range(0, 4π / k, length = 160)
    z_vec = range(0.2, 4π / k - 0.2, length = 18)

    t_obs = Observable(0.0)
    fig = Figure(size = (900, 560), fontsize = 18)
    ax = Axis3(
        fig[1, 1];
        title = "Shear Alfvén Wave (ω = k∥ vₐ, δv ⊥ B₀, δρ = 0)",
        titlesize = 19,
        xlabel = "z (∥ B₀)",
        ylabel = "x (δE)",
        zlabel = "y (δB, δv)",
        azimuth = 1.25π,
        elevation = 0.14π,
        limits = (0, 4π / k, -1.2, 1.2, -1.2, 1.2),
    )

    for x0 in (-0.6, 0.0, 0.6)
        pts = @lift([Point3f(z, x0, 0.8 * sin(k * z - ω * $t_obs)) for z in zs])
        lines!(ax, pts; color = C_BFIELD, linewidth = x0 == 0 ? 3.8 : 2.0)
    end

    b_segs = @lift(
        vcat([[Point3f(z, 0.0, 0.0), Point3f(z, 0.0, 0.8 * sin(k * z - ω * $t_obs))] for z in z_vec]...)
    )
    e_segs = @lift(
        vcat([[Point3f(z, 0.0, 0.0), Point3f(z, -0.8 * sin(k * z - ω * $t_obs), 0.0)] for z in z_vec]...)
    )
    linesegments!(ax, b_segs; color = (C_BFIELD, 0.9), linewidth = 2.5)
    linesegments!(ax, e_segs; color = (C_EFIELD, 0.9), linewidth = 2.5)

    lines!(ax, [Point3f(NaN, NaN, NaN)]; color = C_BFIELD, linewidth = 3.2, label = L"\text{Field lines \& }\delta B_y")
    lines!(ax, [Point3f(NaN, NaN, NaN)]; color = C_EFIELD, linewidth = 3.2, label = L"\text{Electric field }\delta E_x")
    axislegend(ax; position = :rt, framevisible = true, labelsize = 15)

    return _finalize_anim(fig, t_obs, ts, path; fps)
end

"""
    animate_magnetosonic_waves(path=nothing; plasma=normalized_plasma(β_i=0.2, β_e=0.2), θ=π/3, nframes=72, fps=24)

Side-by-side 2D animation of Fast vs. Slow Magnetosonic waves showing magnetic field lines
overlaid on plasma density/pressure compression.
"""
function PlasmaWaves.animate_magnetosonic_waves(
    path::Union{Nothing,AbstractString} = nothing;
    plasma::PlasmaParams = normalized_plasma(; β_i = 0.2, β_e = 0.2),
    θ::Real = π / 3,
    k::Real = 1.0,
    nframes::Int = 72,
    fps::Int = 24,
)
    sol_f = fast_magnetosonic_wave(plasma, k, θ; amp = 0.25)
    kx, kz = sol_f.k[1], sol_f.k[3]

    ϕs = range(0, 2π, length = nframes + 1)[1:end-1]
    ϕ_obs = Observable(0.0)

    xs = range(0, 4π / k, length = 80)
    zs = range(0, 4π / k, length = 80)
    x_lines = range(0.4, 4π / k - 0.4, length = 9)
    z_line = range(0, 4π / k, length = 120)

    fig = Figure(size = (1060, 520), fontsize = 18)

    for (col, (title_tex, sgn_p)) in enumerate([
        (L"\text{Fast Mode } (\delta p \text{ \& } \delta B_\parallel \text{ IN PHASE})", 1.0),
        (L"\text{Slow Mode } (\delta p \text{ \& } \delta B_\parallel \text{ ANTI-PHASE})", -1.0),
    ])
        ax = Axis(
            fig[1, col];
            title = title_tex,
            titlesize = 18,
            xlabel = L"x\; (\perp \mathbf{B}_0)",
            ylabel = L"z\; (\parallel \mathbf{B}_0)",
            aspect = DataAspect(),
            limits = (0, 4π / k, 0, 4π / k),
        )
        p_mat = @lift([sgn_p * cos(kx * x + kz * z - $ϕ_obs) for x in xs, z in zs])
        heatmap!(ax, xs, zs, p_mat; colormap = :RdBu, colorrange = (-1.1, 1.1))

        for x0 in x_lines
            fl_x = @lift([x0 - 0.35 * sin(kx * x0 + kz * z - $ϕ_obs) for z in z_line])
            lines!(ax, fl_x, z_line; color = :black, linewidth = 2.0)
        end
    end

    Colorbar(fig[1, 3]; colormap = :RdBu, limits = (-1, 1), label = L"\delta p / |\delta p|")
    Label(
        fig[0, 1:2],
        LaTeXString("\\text{MHD Fast vs. Slow Magnetosonic Waves } (\\theta = $(Int(round(rad2deg(θ))))^\\circ,\\; \\beta = $(round(plasma_beta(plasma), digits=2)))");
        fontsize = 20,
    )

    return _finalize_anim(fig, ϕ_obs, ϕs, path; fps)
end

"""
    animate_whistler_emic(path=nothing; plasma=normalized_plasma(), nframes=72, fps=24)

Side-by-side 3D helical polarization animation contrasting Whistler (Right-Hand) and EMIC (Left-Hand).
"""
function PlasmaWaves.animate_whistler_emic(
    path::Union{Nothing,AbstractString} = nothing;
    plasma::PlasmaParams = normalized_plasma(),
    k::Real = 1.2,
    nframes::Int = 72,
    fps::Int = 24,
)
    sol_w = whistler_wave(plasma, k; model = :hall)
    sol_e = ion_cyclotron_wave(plasma, k; model = :hall)
    ϕs = range(0, 2π, length = nframes + 1)[1:end-1]
    ϕ_obs = Observable(0.0)
    zs = range(0, 4π / k, length = 140)
    z_vec = range(0.2, 4π / k - 0.2, length = 22)

    fig = Figure(size = (1040, 520), fontsize = 16)
    for (col, (hand, colr, ttl)) in enumerate([
        (1.0, C_WHISTLER, "Whistler (Right-Hand, ω/k = $(round(phase_speed(sol_w), digits=2)) vₐ)"),
        (-1.0, C_EMIC, "EMIC / ICW (Left-Hand, ω/k = $(round(phase_speed(sol_e), digits=2)) vₐ)"),
    ])
        ax = Axis3(
            fig[1, col];
            title = ttl,
            titlesize = 17,
            xlabel = "z (∥ B₀)",
            ylabel = "δB_x",
            zlabel = "δB_y",
            azimuth = 1.22π,
            elevation = 0.14π,
            limits = (0, 4π / k, -1.2, 1.2, -1.2, 1.2),
        )
        lines!(ax, [Point3f(0, 0, 0), Point3f(4π / k, 0, 0)]; color = (:black, 0.4), linestyle = :dash)
        helix = @lift([Point3f(z, cos(k * z - $ϕ_obs), -hand * sin(k * z - $ϕ_obs)) for z in zs])
        lines!(ax, helix; color = colr, linewidth = 3.2)
        segs = @lift(
            vcat([[Point3f(z, 0, 0), Point3f(z, cos(k * z - $ϕ_obs), -hand * sin(k * z - $ϕ_obs))] for z in z_vec]...)
        )
        linesegments!(ax, segs; color = (colr, 0.8), linewidth = 2.0)
        θc = range(0, 2π, length = 60)
        lines!(ax, [Point3f(0, cos(a), sin(a)) for a in θc]; color = (:gray, 0.45), linewidth = 1.4)
        tip = @lift([Point3f(0, cos(-$ϕ_obs), -hand * sin(-$ϕ_obs))])
        scatter!(ax, tip; color = C_FAST, markersize = 13)
    end

    return _finalize_anim(fig, ϕ_obs, ϕs, path; fps)
end

"""
    animate_whistler_packet(path=nothing; nframes=72, fps=24)

Animate a dispersive whistler wave packet (`ω ∝ k²`).
"""
function PlasmaWaves.animate_whistler_packet(
    path::Union{Nothing,AbstractString} = nothing;
    nframes::Int = 72,
    fps::Int = 24,
)
    α = 0.35
    ks = range(1.0, 5.0, length = 60)
    k0 = 2.8
    σk = 0.85
    weights = @. exp(-0.5 * ((ks - k0) / σk)^2)
    weights ./= sum(weights)

    zs = range(-2.0, 28.0, length = 400)
    ts = range(0.0, 8.0, length = nframes)
    t_obs = Observable(0.0)

    fig = Figure(size = (920, 460), fontsize = 18)
    ax = Axis(
        fig[1, 1];
        title = L"\text{Dispersive Whistler Wave Packet } (\omega \propto k^2,\; v_g \propto \sqrt{\omega})",
        titlesize = 19,
        xlabel = L"z / d_i\; (\parallel \mathbf{B}_0)",
        ylabel = L"\delta B_x \text{ (normalized)}",
        limits = (-2.0, 28.0, -1.15, 1.15),
    )

    sig = @lift([sum(w * cos(ki * z - α * ki^2 * $t_obs) for (w, ki) in zip(weights, ks)) for z in zs])
    lines!(ax, zs, sig; color = C_WHISTLER, linewidth = 2.6)

    return _finalize_anim(fig, t_obs, ts, path; fps)
end

"""
    animate_faraday_rotation(path=nothing; plasma=normalized_plasma(mi_me=25.0, vA_c=0.2), nframes=72, fps=24)

3D animation of Faraday rotation of linear polarization along `B₀`.
"""
function PlasmaWaves.animate_faraday_rotation(
    path::Union{Nothing,AbstractString} = nothing;
    plasma::PlasmaParams = normalized_plasma(; mi_me = 25.0, vA_c = 0.2),
    nframes::Int = 72,
    fps::Int = 24,
)
    cut = cold_cutoffs(plasma)
    ω = 1.35 * cut.ω_R
    res0 = faraday_rotation_fields(plasma, ω, [0.0], 0.0)
    L_rot = π / max(abs(res0.rotation_rate), 1e-6)
    zs = range(0, L_rot, length = 180)
    z_vec = range(0, L_rot, length = 28)

    T = 2π / ω
    ts = range(0, T, length = nframes + 1)[1:end-1]
    t_obs = Observable(0.0)

    fig = Figure(size = (900, 540), fontsize = 17)
    ax = Axis3(
        fig[1, 1];
        title = "Faraday Rotation (dφ/dz = (k_L - k_R)/2)",
        titlesize = 19,
        xlabel = "z (∥ B₀)",
        ylabel = "E_x",
        zlabel = "E_y",
        azimuth = 1.24π,
        elevation = 0.14π,
        limits = (0, L_rot, -1.2, 1.2, -1.2, 1.2),
    )
    lines!(ax, [Point3f(0, 0, 0), Point3f(L_rot, 0, 0)]; color = (:black, 0.35), linestyle = :dash)

    pts = @lift(
        begin
            f = faraday_rotation_fields(plasma, ω, zs, $t_obs)
            [Point3f(z, ex, ey) for (z, ex, ey) in zip(zs, f.Ex, f.Ey)]
        end
    )
    segs = @lift(
        begin
            f = faraday_rotation_fields(plasma, ω, z_vec, $t_obs)
            vcat([[Point3f(z, 0, 0), Point3f(z, ex, ey)] for (z, ex, ey) in zip(z_vec, f.Ex, f.Ey)]...)
        end
    )
    lines!(ax, pts; color = C_EFIELD, linewidth = 3.0)
    linesegments!(ax, segs; color = (C_EFIELD, 0.7), linewidth = 1.8)

    return _finalize_anim(fig, t_obs, ts, path; fps)
end

"""
    animate_kaw_vs_saw(path=nothing; nframes=72, fps=24)

Side-by-side animation comparing Ideal SAW (`E_∥ = 0`) vs. KAW (`E_∥ ≠ 0` trapping electrons).
"""
function PlasmaWaves.animate_kaw_vs_saw(
    path::Union{Nothing,AbstractString} = nothing;
    plasma::PlasmaParams = normalized_plasma(; mi_me = 100.0, β_i = 0.2, β_e = 0.2),
    nframes::Int = 72,
    fps::Int = 24,
)
    k_par = 1.0
    riL = larmor_radius(get_ion(plasma), plasma)
    sol_saw = shear_alfven_wave(plasma, k_par, 0.0; amp = 0.25)
    sol_kaw = kinetic_alfven_wave(plasma, k_par, 1.5 / riL; amp = 0.25)

    L = 4π / k_par
    zs = range(0, L, length = 160)
    ω_saw = real(sol_saw.ω)
    ω_kaw = real(sol_kaw.ω)
    vph_kaw = ω_kaw / k_par
    T = 2π / ω_kaw
    ts = range(0, 2 * T, length = nframes)

    z0s = collect(range(0.3, L - 0.3, length = 24))
    vz0s = [vph_kaw * (0.75 + 0.5 * mod(i * 0.37, 1.0)) for i in eachindex(z0s)]
    traj_saw = kaw_particle_trajectories(sol_saw, z0s, vz0s, ts; q_over_m = -25.0)
    traj_kaw = kaw_particle_trajectories(sol_kaw, z0s, vz0s, ts; q_over_m = -25.0)

    idx_obs = Observable(1)
    fig = Figure(size = (1060, 480), fontsize = 18)

    ax1 = Axis(
        fig[1, 1];
        title = L"\text{Ideal SAW } (E_\parallel = 0:\text{ Free Electrons})",
        titlesize = 18,
        xlabel = L"z\; (\parallel \mathbf{B}_0)",
        ylabel = L"\text{Field-line displacement } y",
        limits = (0, L, -1.2, 1.2),
    )
    ax2 = Axis(
        fig[1, 2];
        title = L"\text{KAW } (E_\parallel \neq 0:\text{ Trapped Electrons})",
        titlesize = 18,
        xlabel = L"z\; (\parallel \mathbf{B}_0)",
        ylabel = L"\text{Field-line displacement } y",
        limits = (0, L, -1.2, 1.2),
    )

    fl_saw = @lift([0.55 * sin(k_par * z - ω_saw * ts[$idx_obs]) for z in zs])
    lines!(ax1, zs, fl_saw; color = C_BFIELD, linewidth = 3.4, label = L"\text{Magnetic field line}")
    pz_saw = @lift(mod.(traj_saw.z[:, $idx_obs], L))
    py_saw = @lift([0.55 * sin(k_par * zp - ω_saw * ts[$idx_obs]) for zp in $pz_saw])
    scatter!(ax1, pz_saw, py_saw; color = C_PARTICLE, strokecolor = :black, strokewidth = 1.0, markersize = 12, label = L"\text{Electrons } e^-")
    axislegend(ax1; position = :rt, framevisible = false, labelsize = 14)

    fl_kaw = @lift([0.55 * sin(k_par * z - ω_kaw * ts[$idx_obs]) for z in zs])
    ez_kaw = @lift([0.45 * cos(k_par * z - ω_kaw * ts[$idx_obs]) for z in zs])
    lines!(ax2, zs, ez_kaw; color = (C_EFIELD, 0.7), linewidth = 2.4, linestyle = :dash, label = L"E_\parallel(z, t)")
    lines!(ax2, zs, fl_kaw; color = C_BFIELD, linewidth = 3.4, label = L"\text{Magnetic field line}")
    pz_kaw = @lift(mod.(traj_kaw.z[:, $idx_obs], L))
    py_kaw = @lift([0.55 * sin(k_par * zp - ω_kaw * ts[$idx_obs]) for zp in $pz_kaw])
    scatter!(ax2, pz_kaw, py_kaw; color = C_PARTICLE, strokecolor = :black, strokewidth = 1.0, markersize = 12, label = L"\text{Trapped / passing } e^-")
    axislegend(ax2; position = :rt, framevisible = false, labelsize = 14)

    return _finalize_anim(fig, idx_obs, 1:nframes, path; fps)
end

"""
    animate_hybrid_waves(path=nothing; nframes=72, fps=24)

2D perpendicular electrostatic Lower Hybrid wave animation.
"""
function PlasmaWaves.animate_hybrid_waves(
    path::Union{Nothing,AbstractString} = nothing;
    nframes::Int = 72,
    fps::Int = 24,
)
    k = 1.0
    L = 4π / k
    ϕs = range(0, 2π, length = nframes + 1)[1:end-1]
    ϕ_obs = Observable(0.0)

    xg = range(0.4, L - 0.4, length = 14)
    yg = range(0.5, 4.5, length = 5)
    x0 = [x for x in xg for _ in yg]
    y0 = [y for _ in xg for y in yg]
    xs_line = range(0, L, length = 160)

    fig = Figure(size = (940, 540), fontsize = 18)
    ax1 = Axis(
        fig[1, 1];
        title = L"\text{Lower Hybrid Wave } (\mathbf{k}\perp\mathbf{B}_0):\text{ Ions vs. }\mathbf{E}\times\mathbf{B}\text{ Electrons}",
        titlesize = 18,
        xlabel = L"x\; (\parallel \mathbf{k})",
        ylabel = L"y\; (\perp \mathbf{k},\, \perp \mathbf{B}_0)",
        limits = (0, L, 0, 5.0),
    )
    ax2 = Axis(
        fig[2, 1];
        xlabel = L"x\; (\parallel \mathbf{k})",
        ylabel = L"E_x(x, t)",
        limits = (0, L, -1.2, 1.2),
    )
    Makie.rowsize!(fig.layout, 1, Makie.Relative(0.68))

    xi = @lift(x0 .+ 0.32 .* sin.(k .* x0 .- $ϕ_obs))
    yi = y0
    xe = @lift(x0 .+ 0.08 .* sin.(k .* x0 .- $ϕ_obs))
    ye = @lift(y0 .+ 0.35 .* cos.(k .* x0 .- $ϕ_obs))

    scatter!(ax1, xi, yi; color = C_ALFVEN, markersize = 14, label = L"\text{Ions } (\parallel E_x)")
    scatter!(ax1, xe, ye; color = C_FAST, markersize = 9, label = L"\text{Electrons } (\mathbf{E}\times\mathbf{B})")
    axislegend(ax1; position = :rt, framevisible = true, labelsize = 15)

    ex = @lift(sin.(k .* xs_line .- $ϕ_obs))
    lines!(ax2, xs_line, ex; color = C_EFIELD, linewidth = 3.0)

    return _finalize_anim(fig, ϕ_obs, ϕs, path; fps)
end

"""
    animate_landau_damping(path=nothing; kλD=0.35, nframes=72, fps=24)

Dual-panel Langmuir wave + collisionless Landau damping animation.
"""
function PlasmaWaves.animate_landau_damping(
    path::Union{Nothing,AbstractString} = nothing;
    plasma::PlasmaParams = normalized_plasma(; β_i = 0.05, β_e = 0.5),
    kλD::Real = 0.35,
    nframes::Int = 72,
    fps::Int = 24,
)
    e = get_electron(plasma)
    λDe = debye_length(e, plasma)
    k = kλD / λDe
    sol = langmuir_wave(plasma, k; model = :vlasov, amp = 0.25)
    ωr, γ = real(sol.ω), imag(sol.ω)
    vph = ωr / k
    vth = thermal_speed_1d(e, plasma)

    L = 2π / k
    T = 2π / ωr
    ts = range(0, 3.0 * T, length = nframes)
    xs = range(0, L, length = 160)

    nx_p, nv_p = 20, 9
    x0s = Float64[]
    v0s = Float64[]
    for ix in 1:nx_p, iv in 1:nv_p
        push!(x0s, (ix - 0.5) * L / nx_p)
        push!(v0s, vph + (iv - 5) * 0.22 * vth)
    end
    qm_vis = (e.q / e.m) * (0.12 * vth * ωr / max(abs(e.q / e.m * sol.δE[1]), 1e-12))
    traj = landau_phase_space(sol, x0s, v0s, ts; q_over_m = qm_vis, include_damping = true)

    idx_obs = Observable(1)
    fig = Figure(size = (920, 600), fontsize = 18)

    ax1 = Axis(
        fig[1, 1];
        title = LaTeXString("\\text{Langmuir Landau Damping } (k\\lambda_{De} = $kλD,\\; \\gamma/\\omega_r = $(round(γ/ωr, digits=3)))"),
        titlesize = 18,
        xlabel = L"x / \lambda_{De}",
        ylabel = L"E_x(x, t) / E_0",
        limits = (0, L / λDe, -1.15, 1.15),
    )
    ax2 = Axis(
        fig[2, 1];
        title = L"\text{Electron Phase Space } (x, v_x) \text{ around } v_\phi = \omega_r / k",
        titlesize = 18,
        xlabel = L"x / \lambda_{De}",
        ylabel = L"v_x / v_{th}",
        limits = (0, L / λDe, (vph - 1.4 * vth) / vth, (vph + 1.4 * vth) / vth),
    )
    Makie.rowsize!(fig.layout, 1, Makie.Relative(0.38))

    ex_line = @lift([exp(γ * ts[$idx_obs]) * cos(k * x - ωr * ts[$idx_obs]) for x in xs])
    lines!(ax1, xs ./ λDe, ex_line; color = C_EFIELD, linewidth = 3.0)
    hlines!(ax1, [0.0]; color = (:black, 0.3), linestyle = :dot)

    hlines!(ax2, [vph / vth]; color = (:black, 0.55), linestyle = :dash, linewidth = 2.0, label = L"v_\phi = \omega_r / k")
    px = @lift(mod.(traj.x[:, $idx_obs], L) ./ λDe)
    pv = @lift(traj.v[:, $idx_obs] ./ vth)
    pcols = [v0 < vph ? C_ALFVEN : C_FAST for v0 in v0s]
    scatter!(ax2, px, pv; color = pcols, markersize = 8)
    axislegend(ax2; position = :rt, framevisible = true, labelsize = 15)

    return _finalize_anim(fig, idx_obs, 1:nframes, path; fps)
end

"""
    animate_ion_acoustic(path=nothing; nframes=72, fps=24)

Side-by-side Ion-Acoustic Wave animation comparing `T_e/T_i = 15` vs `T_e/T_i = 1`.
"""
function PlasmaWaves.animate_ion_acoustic(
    path::Union{Nothing,AbstractString} = nothing;
    nframes::Int = 72,
    fps::Int = 24,
)
    p_hot = normalized_plasma(; β_i = 0.02, β_e = 0.30)
    p_eq  = normalized_plasma(; β_i = 0.15, β_e = 0.15)

    k = 1.0
    s_hot = ion_acoustic_wave(p_hot, k; model = :vlasov)
    s_eq  = ion_acoustic_wave(p_eq, k; model = :vlasov)

    L = 4π / k
    xs = range(0, L, length = 160)
    T = 2π / real(s_hot.ω)
    ts = range(0, 2.5 * T, length = nframes)
    t_obs = Observable(0.0)

    fig = Figure(size = (1040, 460), fontsize = 18)
    for (col, (sol, ttl, colr)) in enumerate([
        (s_hot, LaTeXString("T_e / T_i = 15\\; (\\gamma/\\omega_r = $(round(imag(s_hot.ω)/real(s_hot.ω), digits=3)))"), C_SLOW),
        (s_eq,  LaTeXString("T_e / T_i = 1\\; (\\gamma/\\omega_r = $(round(imag(s_eq.ω)/real(s_eq.ω), digits=3)))"), C_FAST),
    ])
        ωr, γ = real(sol.ω), imag(sol.ω)
        ax = Axis(
            fig[1, col];
            title = ttl,
            titlesize = 18,
            xlabel = L"x\; (\parallel \mathbf{k})",
            ylabel = L"\delta n_i / n_0",
            limits = (0, L, -1.15, 1.15),
        )
        dn = @lift([exp(γ * $t_obs) * cos(k * x - ωr * $t_obs) for x in xs])
        lines!(ax, xs, dn; color = colr, linewidth = 3.2)
        hlines!(ax, [0.0]; color = (:black, 0.3), linestyle = :dot)
    end
    Label(fig[0, 1:2], L"\text{Ion-Acoustic Wave: Landau Damping vs. } T_e / T_i"; fontsize = 20)

    return _finalize_anim(fig, t_obs, ts, path; fps)
end

"""
    animate_bernstein_wave(path=nothing; nframes=72, fps=24)

Animated Electron Bernstein Wave (`k_⊥ r_Le = 1.5`, `n = 1` band).
"""
function PlasmaWaves.animate_bernstein_wave(
    path::Union{Nothing,AbstractString} = nothing;
    plasma::PlasmaParams = normalized_plasma(; mi_me = 1836.0, β_i = 0.05, β_e = 0.05, vA_c = 0.05),
    nframes::Int = 72,
    fps::Int = 24,
)
    e = get_electron(plasma)
    Ωe = abs(cyclotron_frequency(e, plasma))
    rLe = larmor_radius(e, plasma)
    k_rL = 1.5
    k_perp = k_rL / rLe
    sol = bernstein_wave(plasma, k_perp; species = :electron, harmonic = 1)
    ω = real(sol.ω)

    L = 4π
    xs = range(0, L, length = 160)
    ϕs = range(0, 2π, length = nframes + 1)[1:end-1]
    ϕ_obs = Observable(0.0)

    fig = Figure(size = (1040, 460), fontsize = 18)
    ax1 = Axis(
        fig[1, 1];
        title = LaTeXString("\\text{EBW } (k_\\perp r_{Le} = $k_rL,\\; \\omega/\\Omega_e = $(round(ω/Ωe, digits=2)))"),
        titlesize = 18,
        xlabel = L"x / r_{Le}",
        ylabel = L"y / r_{Le}\; \& \; E_x",
        limits = (0, L, -2.0, 2.0),
    )

    ex = @lift([1.4 * cos(k_rL * x - $ϕ_obs) for x in xs])
    lines!(ax1, xs, ex; color = C_EFIELD, linewidth = 2.8, label = L"E_x(x, t)")

    xc_list = range(1.5, L - 1.5, length = 6)
    for (j, xc) in enumerate(xc_list)
        θc = range(0, 2π, length = 50)
        lines!(ax1, xc .+ cos.(θc), sin.(θc); color = (:gray, 0.35), linestyle = :dot)
        p_e = @lift([Point2f(xc + cos((Ωe / ω) * $ϕ_obs + j), sin((Ωe / ω) * $ϕ_obs + j))])
        scatter!(ax1, p_e; color = C_ALFVEN, markersize = 12)
    end
    axislegend(ax1; position = :rt, framevisible = false, labelsize = 15)

    ax2 = Axis(
        fig[1, 2];
        title = L"\text{EBW Dispersion } (n = 1, 2, 3)",
        titlesize = 18,
        xlabel = L"k_\perp r_{Le}",
        ylabel = L"\omega / \Omega_e",
    )
    kr_grid = range(0.05, 4.0, length = 80)
    for (n, colr) in zip(1:3, [C_ALFVEN, C_FAST, C_SLOW])
        hlines!(ax2, [Float64(n)]; color = (:black, 0.25), linestyle = :dot)
        ωs = [real(bernstein_wave(plasma, kr / rLe; species = :electron, harmonic = n).ω) / Ωe for kr in kr_grid]
        lines!(ax2, kr_grid, ωs; color = colr, linewidth = 2.6, label = LaTeXString("n = $n"))
    end
    scatter!(ax2, [k_rL], [ω / Ωe]; color = :black, markersize = 13, label = L"\text{Animated mode}")
    axislegend(ax2; position = :rt, framevisible = false, labelsize = 14)

    return _finalize_anim(fig, ϕ_obs, ϕs, path; fps)
end

"""
    animate_phase_mixing(path=nothing; nframes=72, fps=24)

Animate non-uniform slab Shear Alfvén Wave phase mixing.
"""
function PlasmaWaves.animate_phase_mixing(
    path::Union{Nothing,AbstractString} = nothing;
    nframes::Int = 72,
    fps::Int = 24,
)
    xs = range(-4.5, 4.5, length = 400)
    ts = range(0.0, 28.0, length = nframes)
    t_obs = Observable(0.0)

    vA_prof(x) = 1.0 + 0.5 * tanh(x)
    fig = Figure(size = (920, 480), fontsize = 18)
    ax = Axis(
        fig[1, 1];
        title = @lift(LaTeXString("\\text{SAW Phase Mixing } (\\omega_{A0} t = $(round($t_obs, digits=1)))")),
        titlesize = 19,
        xlabel = L"x\text{ (non-uniformity direction)}",
        ylabel = L"B_{1y}(x, t)",
        limits = (-4.5, 4.5, -1.15, 1.15),
    )

    by = @lift(saw_phase_mixing(xs, $t_obs; vA_profile = vA_prof))
    lines!(ax, xs, [exp(-x^2 / 4) for x in xs]; color = (:gray, 0.6), linestyle = :dash, linewidth = 2.2, label = L"\text{Envelope } B_{1y}(x, 0)")
    lines!(ax, xs, by; color = C_ALFVEN, linewidth = 2.8, label = L"B_{1y}(x, t)")
    axislegend(ax; position = :rt, framevisible = false, labelsize = 15)

    return _finalize_anim(fig, t_obs, ts, path; fps)
end

end # module PlasmaWavesMakieExt
