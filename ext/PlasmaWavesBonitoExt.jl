module PlasmaWavesBonitoExt

using Bonito
using Bonito: App, DOM, Routes, Slider, Styles, export_static, Server
using WGLMakie
using WGLMakie: Figure, Axis, Axis3, Observable, @lift, Point3f, lines!, linesegments!, scatter!, hlines!, axislegend, DataAspect
using PlasmaWaves

const NAV_STYLE = Styles(
    "display" => "flex",
    "gap" => "1.2rem",
    "padding" => "0.9rem 1.5rem",
    "background" => "#0f172a",
    "color" => "#f8fafc",
    "font-family" => "system-ui, -apple-system, sans-serif",
    "align-items" => "center",
    "flex-wrap" => "wrap",
)

const LINK_STYLE = Styles(
    "color" => "#93c5fd",
    "text-decoration" => "none",
    "font-weight" => "600",
    "font-size" => "0.95rem",
)

const PAGE_STYLE = Styles(
    "max-width" => "1080px",
    "margin" => "1.5rem auto",
    "padding" => "0 1.25rem 3rem 1.25rem",
    "font-family" => "system-ui, -apple-system, sans-serif",
    "line-height" => "1.6",
    "color" => "#1e293b",
)

const CARD_STYLE = Styles(
    "border" => "1px solid #e2e8f0",
    "border-radius" => "10px",
    "padding" => "1.25rem",
    "margin" => "1rem 0",
    "background" => "#ffffff",
    "box-shadow" => "0 1px 3px rgba(0,0,0,0.06)",
)

_norm_base(base_url::AbstractString) = rstrip(String(base_url), '/')
_route_href(base::AbstractString, path::AbstractString) = isempty(base) ? path : string(base, path)

function _navbar(base::AbstractString)
    return DOM.nav(
        DOM.span("PlasmaWaves.jl"; style = Styles("font-weight" => "800", "font-size" => "1.1rem", "margin-right" => "1rem", "color" => "#ffffff")),
        DOM.a("Overview"; href = _route_href(base, "/"), style = LINK_STYLE),
        DOM.a("MHD & Alfvén"; href = _route_href(base, "/mhd/"), style = LINK_STYLE),
        DOM.a("Cold & Whistlers"; href = _route_href(base, "/cold/"), style = LINK_STYLE),
        DOM.a("KAW & IAW"; href = _route_href(base, "/kaw/"), style = LINK_STYLE),
        DOM.a("Kinetic & Bernstein"; href = _route_href(base, "/kinetic/"), style = LINK_STYLE),
        DOM.a("Theory vs. PlasmaBO"; href = _route_href(base, "/confirm/"), style = LINK_STYLE);
        style = NAV_STYLE,
    )
end

function _page(base::AbstractString, title_str::AbstractString, children...)
    return DOM.div(
        _navbar(base),
        DOM.main(
            DOM.h1(title_str; style = Styles("margin-bottom" => "0.5rem")),
            children...;
            style = PAGE_STYLE,
        ),
    )
end

function _media_card(title_str::AbstractString, desc::AbstractString, media_url::AbstractString)
    ext = lowercase(splitext(media_url)[2])
    media_el = if ext in (".mp4", ".webm")
        mime = ext == ".mp4" ? "video/mp4" : "video/webm"
        DOM.video(
            DOM.source(; src = media_url, type = mime);
            autoplay = true,
            loop = true,
            muted = true,
            playsinline = true,
            controls = true,
            style = Styles("max-width" => "100%", "border-radius" => "6px"),
        )
    else
        DOM.img(; src = media_url, style = Styles("max-width" => "100%", "border-radius" => "6px"))
    end
    return DOM.div(
        DOM.h3(title_str),
        DOM.p(desc),
        media_el;
        style = CARD_STYLE,
    )
end

function _home_app(base::AbstractString)
    return App() do session
        _page(
            base,
            "Plasma Waves Interactive Explorer",
            DOM.p(
                "Explore the propagation, polarization, dispersion relations, and kinetic wave–particle physics of plasma waves — combining analytical theory, symbolic derivations, PlasmaBO.jl kinetic confirmation, and interactive visualizations.",
            ),
            DOM.div(
                DOM.h2("Wave Families"),
                DOM.ul(
                    DOM.li(DOM.a("MHD & Alfvén Waves"; href = _route_href(base, "/mhd/")), " — Shear Alfvén, Fast & Slow Magnetosonic, Friedrich diagrams, Phase Mixing."),
                    DOM.li(DOM.a("Cold Plasma & Whistlers"; href = _route_href(base, "/cold/")), " — Stix CPDR, Whistler vs. EMIC circular polarization, CMA diagram, Faraday rotation."),
                    DOM.li(DOM.a("Kinetic & Inertial Alfvén Waves"; href = _route_href(base, "/kaw/")), " — Finite Larmor radius / electron inertia dispersion, parallel E_∥, and electron trapping."),
                    DOM.li(DOM.a("Kinetic, Landau & Bernstein Waves"; href = _route_href(base, "/kinetic/")), " — Vlasov plasma dispersion Z(ζ), Langmuir & Ion-Acoustic Landau damping, EBW/IBW harmonics."),
                    DOM.li(DOM.a("Theory vs. PlasmaBO Confirmation"; href = _route_href(base, "/confirm/")), " — Direct quantitative comparison of analytical dispersion relations against full hot magnetized kinetic PlasmaBO.jl solutions."),
                );
                style = CARD_STYLE,
            ),
        )
    end
end

function _mhd_app(base::AbstractString; gallery_url::Union{Nothing,String} = nothing, media_ext::String = "mp4")
    return App() do session
        s_vs = Slider(0.2:0.05:2.0; value = 0.6)
        θs = range(0, 2π, length = 241)

        fig = Figure(size = (560, 440))
        ax = Axis(fig[1, 1]; title = "Interactive MHD Friedrich Phase-Velocity Diagram", xlabel = "v_x / vₐ", ylabel = "v_z / vₐ", aspect = DataAspect())
        vf_x = @lift([friedrich_phase_speeds(1.0, $s_vs, θ).v_fast * sin(θ) for θ in θs])
        vf_z = @lift([friedrich_phase_speeds(1.0, $s_vs, θ).v_fast * cos(θ) for θ in θs])
        va_x = [friedrich_phase_speeds(1.0, 0.6, θ).v_alfven * sin(θ) for θ in θs]
        va_z = [friedrich_phase_speeds(1.0, 0.6, θ).v_alfven * cos(θ) for θ in θs]
        vs_x = @lift([friedrich_phase_speeds(1.0, $s_vs, θ).v_slow * sin(θ) for θ in θs])
        vs_z = @lift([friedrich_phase_speeds(1.0, $s_vs, θ).v_slow * cos(θ) for θ in θs])

        lines!(ax, vf_x, vf_z; color = :crimson, linewidth = 2.5, label = "Fast")
        lines!(ax, va_x, va_z; color = :dodgerblue, linewidth = 2.5, linestyle = :dash, label = "Alfvén")
        lines!(ax, vs_x, vs_z; color = :forestgreen, linewidth = 2.5, label = "Slow")
        axislegend(ax; position = :rb, framevisible = false)

        fig3d = animate_alfven_wave(nothing).fig

        cards = Any[
            DOM.div(
                DOM.h2("Friedrich Diagram Explorer"),
                DOM.p("Adjust vₛ / vₐ (live server mode) to transition between low-β and high-β regimes:"),
                DOM.div("vₛ / vₐ = ", s_vs),
                fig;
                style = CARD_STYLE,
            ),
            DOM.div(
                DOM.h2("Interactive 3D Shear Alfvén Wave (Drag to Rotate)"),
                fig3d;
                style = CARD_STYLE,
            ),
        ]
        if gallery_url !== nothing
            push!(cards, _media_card("Shear Alfvén Wave", "Transverse plucked magnetic field lines.", "$gallery_url/alfven_wave.$media_ext"))
            push!(cards, _media_card("Fast vs. Slow Magnetosonic Waves", "In-phase vs. anti-phase pressure and magnetic compressions.", "$gallery_url/magnetosonic_waves.$media_ext"))
            push!(cards, _media_card("Non-Uniform SAW Phase Mixing", "Perpendicular scale cascading |k_x| ∝ |ω_A'(x)| t.", "$gallery_url/phase_mixing.$media_ext"))
        end
        return _page(base, "MHD & Alfvén Waves", cards...)
    end
end

function _cold_app(base::AbstractString; gallery_url::Union{Nothing,String} = nothing, media_ext::String = "mp4")
    return App() do session
        fig_disp = plot_dispersion_parallel()
        fig_cma = plot_cma_diagram()
        fig_whis = animate_whistler_emic(nothing).fig

        cards = Any[
            DOM.div(DOM.h2("Parallel Cold-Plasma Dispersion (R & L Modes)"), fig_disp; style = CARD_STYLE),
            DOM.div(DOM.h2("Clemmow–Mullaly–Allis (CMA) Diagram"), fig_cma; style = CARD_STYLE),
            DOM.div(DOM.h2("Interactive 3D Whistler vs. EMIC Polarization (Drag to Rotate)"), fig_whis; style = CARD_STYLE),
        ]
        if gallery_url !== nothing
            push!(cards, _media_card("Whistler vs. EMIC Helical Polarization", "Right-hand (electron) vs. left-hand (ion) circular polarization.", "$gallery_url/whistler_emic.$media_ext"))
            push!(cards, _media_card("Dispersive Whistler Wave Packet", "High frequencies travel faster than low frequencies (v_g ∝ √ω).", "$gallery_url/whistler_packet.$media_ext"))
            push!(cards, _media_card("Faraday Rotation", "Rotation of the linear polarization plane along B₀.", "$gallery_url/faraday_rotation.$media_ext"))
            push!(cards, _media_card("Lower Hybrid Wave", "Unmagnetized ions vs. E×B drifting electrons.", "$gallery_url/hybrid_waves.$media_ext"))
        end
        return _page(base, "Cold Plasma, Hall-MHD & Whistler Waves", cards...)
    end
end

function _kaw_app(base::AbstractString; gallery_url::Union{Nothing,String} = nothing, media_ext::String = "mp4")
    return App() do session
        fig_kaw = plot_kaw_iaw_dispersion()
        cards = Any[
            DOM.div(DOM.h2("KAW & IAW Dispersion and Parallel Electric Field"), fig_kaw; style = CARD_STYLE),
        ]
        if gallery_url !== nothing
            push!(cards, _media_card("KAW vs. Ideal SAW Electron Trapping", "Resonant electrons trapped in KAW's finite E_∥ potential troughs.", "$gallery_url/kaw_vs_saw.$media_ext"))
        end
        return _page(base, "Kinetic & Inertial Alfvén Waves", cards...)
    end
end

function _kinetic_app(base::AbstractString; gallery_url::Union{Nothing,String} = nothing, media_ext::String = "mp4")
    return App() do session
        fig_landau = plot_landau_dispersion()
        fig_ebw = plot_bernstein_dispersion()
        cards = Any[
            DOM.div(DOM.h2("Langmuir & Ion-Acoustic Vlasov Dispersion"), fig_landau; style = CARD_STYLE),
            DOM.div(DOM.h2("Electron Bernstein Wave Harmonics"), fig_ebw; style = CARD_STYLE),
        ]
        if gallery_url !== nothing
            push!(cards, _media_card("Langmuir Wave & Landau Damping", "Damped E_x(x, t) and (x, v) phase-space trapping.", "$gallery_url/landau_damping.$media_ext"))
            push!(cards, _media_card("Ion-Acoustic Wave Damping", "Tₑ/Tᵢ = 15 (propagating) vs. Tₑ/Tᵢ = 1 (strongly damped).", "$gallery_url/ion_acoustic.$media_ext"))
            push!(cards, _media_card("Electron Bernstein Wave", "Gyrating electrons coupled to perpendicular electrostatic waves.", "$gallery_url/bernstein_wave.$media_ext"))
        end
        return _page(base, "Kinetic, Landau & Bernstein Waves", cards...)
    end
end

function _confirm_app(base::AbstractString)
    return App() do session
        fig_alf = plot_confirmation(confirm_alfven())
        fig_lang = plot_confirmation(confirm_langmuir())
        return _page(
            base,
            "Theory vs. Kinetic Confirmation (PlasmaBO.jl)",
            DOM.p(
                "Compare analytical textbook dispersion relations directly against full hot magnetized kinetic solutions from PlasmaBO.jl (BOHH), revealing where fluid/asymptotic models hold and where Hall dispersion, cyclotron damping, or Landau damping take over.",
            ),
            DOM.div(DOM.h2("Ideal MHD Alfvén vs. PlasmaBO (Hall Splitting at k dᵢ ∼ 1)"), fig_alf; style = CARD_STYLE),
            DOM.div(DOM.h2("Bohm-Gross Langmuir vs. PlasmaBO (Landau Damping Onset)"), fig_lang; style = CARD_STYLE),
        )
    end
end

"""
    plasma_waves_routes(; gallery_prefix=nothing, base_url="", media_ext="mp4")

Build the multi-page `Bonito.Routes` web application for `PlasmaWaves.jl`.
- `base_url`: subpath prefix for GitHub Pages project sites (e.g. `"/PlasmaWaves"`).
- `gallery_prefix`: folder name of pre-rendered animations (e.g. `"gallery"`).
- `media_ext`: animation file extension (`"mp4"`, `"webm"`, or `"gif"`).
"""
function PlasmaWaves.plasma_waves_routes(;
    gallery_prefix::Union{Nothing,String} = nothing,
    base_url::AbstractString = "",
    media_ext::AbstractString = "mp4",
)
    base = _norm_base(base_url)
    gallery_url = gallery_prefix === nothing ? nothing : _route_href(base, "/" * lstrip(gallery_prefix, '/'))
    ext = String(lstrip(String(media_ext), '.'))
    return Routes(
        "/" => _home_app(base),
        "/mhd" => _mhd_app(base; gallery_url, media_ext = ext),
        "/cold" => _cold_app(base; gallery_url, media_ext = ext),
        "/kaw" => _kaw_app(base; gallery_url, media_ext = ext),
        "/kinetic" => _kinetic_app(base; gallery_url, media_ext = ext),
        "/confirm" => _confirm_app(base),
    )
end

"""
    build_website(output_dir::AbstractString; gallery_prefix=nothing, base_url="", media_ext="mp4")

Export the multi-page `PlasmaWaves.jl` website to static HTML in `output_dir`
using `Bonito.export_static`. Creates `.nojekyll` for clean GitHub Pages hosting.
"""
function PlasmaWaves.build_website(
    output_dir::AbstractString;
    gallery_prefix::Union{Nothing,String} = nothing,
    base_url::AbstractString = "",
    media_ext::AbstractString = "mp4",
)
    mkpath(output_dir)
    routes = plasma_waves_routes(; gallery_prefix, base_url, media_ext)
    export_static(output_dir, routes)
    touch(joinpath(output_dir, ".nojekyll"))
    return output_dir
end

"""
    serve_website(; host="127.0.0.1", port=8080, gallery_prefix=nothing, base_url="", media_ext="mp4")

Launch a live `Bonito.Server` serving the interactive `PlasmaWaves.jl` web application
with full Julia-backed slider reactivity.
"""
function PlasmaWaves.serve_website(;
    host::AbstractString = "127.0.0.1",
    port::Int = 8080,
    gallery_prefix::Union{Nothing,String} = nothing,
    base_url::AbstractString = "",
    media_ext::AbstractString = "mp4",
)
    routes = plasma_waves_routes(; gallery_prefix, base_url, media_ext)
    server = Server(routes, host, port)
    @info "PlasmaWaves.jl website running at http://$host:$port"
    return server
end

end # module PlasmaWavesBonitoExt
