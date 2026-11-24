"""
Generate the full gallery of PlasmaWaves.jl static diagnostic figures (.png)
and animated wave propagation demos (.mp4 by default; also supports .webm or .gif).

Usage:
    julia --project=examples examples/generate_gallery.jl [output_dir] [format]

Examples:
    julia --project=examples examples/generate_gallery.jl gallery mp4
    julia --project=examples examples/generate_gallery.jl gallery gif
"""

using CairoMakie
using PlasmaWaves

function generate_gallery(
    outdir::AbstractString = joinpath(@__DIR__, "..", "gallery");
    format::AbstractString = "mp4",
    nframes::Int = 72,
    fps::Int = 24,
)
    ext = lstrip(lowercase(String(format)), '.')
    mkpath(outdir)
    @info "Generating PlasmaWaves.jl gallery (.$ext animations) into $outdir ..."

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
        ("confirm_whistler.png", () -> plot_confirmation(confirm_whistler())),
        ("confirm_langmuir.png", () -> plot_confirmation(confirm_langmuir())),
    ]
    for (fname, fn) in static_figs
        p = joinpath(outdir, fname)
        @info "  Rendering $fname ..."
        fig = fn()
        save(p, fig; px_per_unit = 2)
    end

    # Animated wave propagation demos
    animations = [
        ("alfven_wave", animate_alfven_wave),
        ("magnetosonic_waves", animate_magnetosonic_waves),
        ("whistler_emic", animate_whistler_emic),
        ("whistler_packet", animate_whistler_packet),
        ("faraday_rotation", animate_faraday_rotation),
        ("kaw_vs_saw", animate_kaw_vs_saw),
        ("hybrid_waves", animate_hybrid_waves),
        ("landau_damping", animate_landau_damping),
        ("ion_acoustic", animate_ion_acoustic),
        ("bernstein_wave", animate_bernstein_wave),
        ("phase_mixing", animate_phase_mixing),
    ]
    for (stem, fn) in animations
        fname = "$stem.$ext"
        p = joinpath(outdir, fname)
        @info "  Animating $fname ($nframes frames @ $fps fps) ..."
        fn(p; nframes, fps)
    end

    @info "Done! All figures and animations saved to $outdir."
    return outdir
end

if abspath(PROGRAM_FILE) == @__FILE__
    outdir = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "..", "gallery")
    fmt = length(ARGS) >= 2 ? ARGS[2] : "mp4"
    generate_gallery(outdir; format = fmt)
end
