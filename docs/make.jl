using Documenter
using Bonito
using CairoMakie
using Symbolics
using Latexify
using PlasmaWaves

function _animate_friedrich_sweep(dst::AbstractString; nframes::Int = 72, fps::Int = 24)
    CairoMakie.activate!()
    set2 = CairoMakie.to_colormap(:Set2_8)
    θs = range(0, 2π; length = 181)
    vs_obs = Observable(0.25)

    fig = Figure(size = (680, 520), fontsize = 18)
    ax = Axis(
        fig[1, 1];
        title = L"\text{MHD Friedrich Phase-Velocity Diagram}",
        titlesize = 19,
        xlabel = L"v_x / v_A",
        ylabel = L"v_z / v_A",
        aspect = DataAspect(),
        limits = (-2.2, 2.2, -2.2, 2.2),
    )

    vf_x = @lift([friedrich_phase_speeds(1.0, $vs_obs, θ).v_fast * sin(θ) for θ in θs])
    vf_z = @lift([friedrich_phase_speeds(1.0, $vs_obs, θ).v_fast * cos(θ) for θ in θs])
    va_x = [friedrich_phase_speeds(1.0, 1.0, θ).v_alfven * sin(θ) for θ in θs]
    va_z = [friedrich_phase_speeds(1.0, 1.0, θ).v_alfven * cos(θ) for θ in θs]
    vs_x = @lift([friedrich_phase_speeds(1.0, $vs_obs, θ).v_slow * sin(θ) for θ in θs])
    vs_z = @lift([friedrich_phase_speeds(1.0, $vs_obs, θ).v_slow * cos(θ) for θ in θs])

    lines!(ax, vf_x, vf_z; color = set2[2], linewidth = 3.0, label = L"\text{Fast}")
    lines!(ax, va_x, va_z; color = set2[3], linewidth = 3.0, linestyle = :dash, label = L"\text{Alfvén}")
    lines!(ax, vs_x, vs_z; color = set2[1], linewidth = 3.0, label = L"\text{Slow}")
    axislegend(ax; position = :rb, framevisible = false, labelsize = 15)
    text!(
        ax,
        -2.05,
        1.85;
        text = @lift("vₛ / vₐ = $(round($vs_obs; digits = 2))"),
        fontsize = 17,
        color = :black,
    )

    vs_vals = 1.05 .- 0.80 .* cos.(range(0, 2π; length = nframes + 1)[1:(end - 1)])
    record(fig, dst, vs_vals; framerate = fps) do vs
        vs_obs[] = vs
    end
    return dst
end

# Ensure gallery MP4 animations exist in docs/src/assets/gallery/
const GALLERY_SRC = joinpath(@__DIR__, "..", "gallery")
const GALLERY_DST = joinpath(@__DIR__, "src", "assets", "gallery")
mkpath(GALLERY_DST)

const ANIMATIONS = [
    ("friedrich_sweep.mp4", _animate_friedrich_sweep),
    ("alfven_wave.mp4", animate_alfven_wave),
    ("magnetosonic_waves.mp4", animate_magnetosonic_waves),
    ("whistler_emic.mp4", animate_whistler_emic),
    ("whistler_packet.mp4", animate_whistler_packet),
    ("faraday_rotation.mp4", animate_faraday_rotation),
    ("kaw_vs_saw.mp4", animate_kaw_vs_saw),
    ("hybrid_waves.mp4", animate_hybrid_waves),
    ("landau_damping.mp4", animate_landau_damping),
    ("ion_acoustic.mp4", animate_ion_acoustic),
    ("bernstein_wave.mp4", animate_bernstein_wave),
    ("phase_mixing.mp4", animate_phase_mixing),
]

CairoMakie.activate!()
for (fname, fn) in ANIMATIONS
    src = joinpath(GALLERY_SRC, fname)
    dst = joinpath(GALLERY_DST, fname)
    if isfile(src)
        cp(src, dst; force = true)
    else
        @info "Generating $fname for docs ..."
        fn(dst; nframes = 72, fps = 24)
    end
end

CairoMakie.activate!()

DocMeta.setdocmeta!(PlasmaWaves, :DocTestSetup, :(using PlasmaWaves); recursive = true)

makedocs(
    sitename = "PlasmaWaves.jl",
    modules = [PlasmaWaves],
    authors = "Hongyang Zhou <hyzhou@umich.edu>",
    checkdocs = :exports,
    format = Bonito.DocumenterBonito(
        repo = "github.com/henry2004y/PlasmaWaves",
        devbranch = "main",
        version = "v0.1.0",
        home_page = "index.md",
        home = (
            name = "PlasmaWaves.jl",
            text = "Plasma Wave Theory, Confirmation & Visual Demos",
            tagline = "Analytical plasma wave physics, Symbolics.jl derivations, PlasmaBO.jl kinetic confirmation, and Makie figures & MP4 animations.",
            actions = [
                (text = "Explore Wave Families", link = "mhd.html", theme = "brand"),
                (text = "Theory vs. PlasmaBO", link = "confirm.html", theme = "alt"),
                (text = "GitHub", link = "https://github.com/henry2004y/PlasmaWaves", theme = "alt"),
            ],
            features = [
                (icon = "🌊", title = "MHD to Kinetic Hierarchy", details = "Shear Alfvén, Fast/Slow Magnetosonic, Whistler, EMIC, KAW/IAW, Landau damping, and Bernstein waves.", link = "mhd.html"),
                (icon = "🔬", title = "PlasmaBO.jl Confirmation", details = "Benchmark textbook dispersion relations against tracked hot magnetized kinetic branches.", link = "confirm.html"),
                (icon = "∑", title = "Symbolic Derivations", details = "Derive dispersion matrices, group velocities, and KeyNotes-ready LaTeX with Symbolics.jl.", link = "symbolic.html"),
                (icon = "🎬", title = "Figures & MP4 Animations", details = "Publication-quality CairoMakie dispersion diagrams and 12 looping MP4 wave propagation animations.", link = "mhd.html"),
            ],
        ),
    ),
    pages = [
        "Home" => "index.md",
        "Wave Families" => [
            "MHD & Alfvén Waves" => "mhd.md",
            "Cold Plasma & Whistlers" => "cold.md",
            "Kinetic & Inertial Alfvén" => "kaw.md",
            "Kinetic, Landau & Bernstein" => "kinetic.md",
        ],
        "Theory vs. PlasmaBO" => "confirm.md",
        "Symbolic Derivations" => "symbolic.md",
        "API Reference" => "api.md",
    ],
)

touch(joinpath(@__DIR__, "build", ".nojekyll"))
