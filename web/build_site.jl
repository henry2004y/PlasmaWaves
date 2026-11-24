"""
Export the PlasmaWaves.jl multi-page Bonito + WGLMakie website to static HTML
(ready for GitHub Pages deployment).

Usage:
    julia --project=web web/build_site.jl [output_dir]
"""

using Bonito
using WGLMakie
using PlasmaWaves

outdir = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "..", "site")
gallery_src = joinpath(@__DIR__, "..", "gallery")
gallery_prefix = isdir(gallery_src) ? "gallery" : nothing

build_website(outdir; gallery_prefix)

if gallery_prefix !== nothing
    cp(gallery_src, joinpath(outdir, "gallery"); force = true)
end

@info "Static PlasmaWaves.jl website exported to $outdir"
