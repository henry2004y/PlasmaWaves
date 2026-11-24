"""
Export the PlasmaWaves.jl multi-page Bonito + WGLMakie website to static HTML
(ready for GitHub Pages deployment).

Usage:
    julia --project=web web/build_site.jl [output_dir]

Environment variables:
    BASE_URL   Subpath prefix (e.g. "/PlasmaWaves" for GitHub Pages; default "")
    MEDIA_EXT  Animation format in gallery/ ("mp4", "webm", or "gif"; auto-detected)
"""

using Bonito
using WGLMakie
using PlasmaWaves

outdir = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "..", "site")
base_url = get(ENV, "BASE_URL", "")
gallery_src = joinpath(@__DIR__, "..", "gallery")
gallery_prefix = isdir(gallery_src) ? "gallery" : nothing

media_ext = get(ENV, "MEDIA_EXT", "")
if isempty(media_ext) && isdir(gallery_src)
    if isfile(joinpath(gallery_src, "alfven_wave.mp4"))
        media_ext = "mp4"
    elseif isfile(joinpath(gallery_src, "alfven_wave.webm"))
        media_ext = "webm"
    else
        media_ext = "gif"
    end
elseif isempty(media_ext)
    media_ext = "mp4"
end

build_website(outdir; gallery_prefix, base_url, media_ext)

if gallery_prefix !== nothing
    cp(gallery_src, joinpath(outdir, "gallery"); force = true)
end

@info "Static PlasmaWaves.jl website exported to $outdir (base_url=$(repr(base_url)), media_ext=$(repr(media_ext)))"
