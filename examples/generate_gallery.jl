"""
Generate the full gallery of PlasmaWaves.jl static diagnostic figures (.png)
and animated wave propagation demos (.gif).

Usage:
    julia --project=. examples/generate_gallery.jl [output_dir]
"""

using CairoMakie
using PlasmaWaves

function generate_gallery(outdir::AbstractString = joinpath(@__DIR__, "..", "gallery");
                          nframes::Int = 72, fps::Int = 24)
    mkpath(outdir)
    @info "Generating PlasmaWaves.jl gallery into $outdir ..."

    # Static diagnostic figures
    static_figs = [
        ("friedrich.png", plot_friedrich),
        ("dispersion_parallel.png", plot_dispersion_parallel),
        ("cma_diagram.png", plot_cma_diagram),
        ("kaw_iaw_dispersion.png", plot_kaw_iaw_dispersion),
        ("landau_dispersion.png", plot_landau_dispersion),
        ("bernstein_dispersion.png", plot_bernstein_dispersion),
        ("transport_ratios.png", plot_transport_ratios),
        ("confirm_alfven.png", () -> plot_confirmation(confirm_alfven())),
        ("confirm_langmuir.png", () -> plot_confirmation(confirm_langmuir())),
    ]
    for (fname, fn) in static_figs
        p = joinpath(outdir, fname)
        @info "  Rendering $fname ..."
        fig = fn()
        save(p, fig; px_per_unit = 2)
    end

    # Animated wave propagation GIFs
    animations = [
        ("alfven_wave.gif", animate_alfven_wave),
        ("magnetosonic_waves.gif", animate_magnetosonic_waves),
        ("whistler_emic.gif", animate_whistler_emic),
        ("whistler_packet.gif", animate_whistler_packet),
        ("faraday_rotation.gif", animate_faraday_rotation),
        ("kaw_vs_saw.gif", animate_kaw_vs_saw),
        ("hybrid_waves.gif", animate_hybrid_waves),
        ("landau_damping.gif", animate_landau_damping),
        ("ion_acoustic.gif", animate_ion_acoustic),
        ("bernstein_wave.gif", animate_bernstein_wave),
        ("phase_mixing.gif", animate_phase_mixing),
    ]
    for (fname, fn) in animations
        p = joinpath(outdir, fname)
        @info "  Animating $fname ($nframes frames @ $fps fps) ..."
        fn(p; nframes, fps)
    end

    @info "Done! All figures and GIFs saved to $outdir."
    return outdir
end

if abspath(PROGRAM_FILE) == @__FILE__
    outdir = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "..", "gallery")
    generate_gallery(outdir)
end
