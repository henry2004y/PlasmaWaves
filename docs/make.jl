using Documenter
using Bonito
using CairoMakie
using WGLMakie
using Symbolics
using Latexify
using PlasmaWaves

# Ensure gallery MP4 animations exist in docs/src/assets/gallery/
const GALLERY_SRC = joinpath(@__DIR__, "..", "gallery")
const GALLERY_DST = joinpath(@__DIR__, "src", "assets", "gallery")
mkpath(GALLERY_DST)

const ANIMATIONS = [
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
            text = "Plasma Wave Theory, Confirmation & Interactive Demos",
            tagline = "Analytical plasma wave physics, Symbolics.jl derivations, PlasmaBO.jl kinetic confirmation, and interactive Makie/Bonito visualizations.",
            actions = [
                (text = "Explore Wave Families", link = "mhd.html", theme = "brand"),
                (text = "Theory vs. PlasmaBO", link = "confirm.html", theme = "alt"),
                (text = "GitHub", link = "https://github.com/henry2004y/PlasmaWaves", theme = "alt"),
            ],
            features = [
                (icon = "🌊", title = "MHD to Kinetic Hierarchy", details = "Shear Alfvén, Fast/Slow Magnetosonic, Whistler, EMIC, KAW/IAW, Landau damping, and Bernstein waves.", link = "mhd.html"),
                (icon = "🔬", title = "PlasmaBO.jl Confirmation", details = "Benchmark textbook dispersion relations against tracked hot magnetized kinetic branches.", link = "confirm.html"),
                (icon = "∑", title = "Symbolic Derivations", details = "Derive dispersion matrices, group velocities, and KeyNotes-ready LaTeX with Symbolics.jl.", link = "symbolic.html"),
                (icon = "🎬", title = "Interactive 3D & MP4 Demos", details = "Rotatable WebGL WGLMakie scenes, offline pre-recorded sliders, and crisp MP4 wave animations.", link = "cold.html"),
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
