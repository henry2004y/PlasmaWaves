# MHD & Alfvén Waves

In ideal magnetohydrodynamics ($\omega \ll \Omega_i$, $k d_i \ll 1$), linear perturbations around a uniform magnetized plasma $\mathbf{B}_0 = B_0\hat{\mathbf{z}}$ decouple into three eigenmodes:

- **Shear Alfvén Wave (SAW)**: $\omega^2 = k_\parallel^2 v_A^2$, incompressible ($\delta n = \delta p = \delta B_\parallel = 0$), energy guided strictly along $\mathbf{B}_0$ ($\mathbf{v}_g = \pm v_A \hat{\mathbf{z}}$).
- **Fast & Slow Magnetosonic Waves**:
  $$\frac{\omega^2}{k^2} = \frac{1}{2}(v_A^2 + v_s^2) \pm \frac{1}{2}\sqrt{(v_A^2 + v_s^2)^2 - 4 v_A^2 v_s^2 \cos^2\theta}$$
  Thermal and magnetic pressures reinforce ($\delta p \cdot \delta B_\parallel > 0$) in the **fast** mode and oppose ($\delta p \cdot \delta B_\parallel < 0$) in the **slow** mode.

## Static Friedrich Diagrams & Transport Ratios

```@example mhd_static
using CairoMakie, PlasmaWaves
CairoMakie.activate!()
plot_friedrich()
```

```@example mhd_transport
using CairoMakie, PlasmaWaves
CairoMakie.activate!()
plot_transport_ratios()
```

## Interactive Friedrich Explorer (Offline Static Slider)

Drag the $v_s / v_A$ slider below to transition between low-$\beta$ ($v_s < v_A$) and high-$\beta$ ($v_s > v_A$) regimes—this widget is pre-recorded with `Bonito.record_states` so it works completely offline on GitHub Pages!

```@example mhd_slider
using Bonito, WGLMakie, PlasmaWaves
WGLMakie.activate!()

App() do session
    s_vs = Slider(0.2:0.15:2.0; value = 0.65)
    θs = range(0, 2π, length = 181)
    set2 = WGLMakie.to_colormap(:Set2_8)

    fig = Figure(size = (640, 480), fontsize = 18)
    ax = Axis(
        fig[1, 1];
        title = L"\text{MHD Friedrich Phase-Velocity Diagram}",
        titlesize = 19,
        xlabel = L"v_x / v_A",
        ylabel = L"v_z / v_A",
        aspect = DataAspect(),
        limits = (-2.2, 2.2, -2.2, 2.2),
    )
    vf_x = @lift([friedrich_phase_speeds(1.0, $s_vs, θ).v_fast * sin(θ) for θ in θs])
    vf_z = @lift([friedrich_phase_speeds(1.0, $s_vs, θ).v_fast * cos(θ) for θ in θs])
    va_x = [friedrich_phase_speeds(1.0, 0.65, θ).v_alfven * sin(θ) for θ in θs]
    va_z = [friedrich_phase_speeds(1.0, 0.65, θ).v_alfven * cos(θ) for θ in θs]
    vs_x = @lift([friedrich_phase_speeds(1.0, $s_vs, θ).v_slow * sin(θ) for θ in θs])
    vs_z = @lift([friedrich_phase_speeds(1.0, $s_vs, θ).v_slow * cos(θ) for θ in θs])

    lines!(ax, vf_x, vf_z; color = set2[2], linewidth = 3.0, label = L"\text{Fast}")
    lines!(ax, va_x, va_z; color = set2[3], linewidth = 3.0, linestyle = :dash, label = L"\text{Alfvén}")
    lines!(ax, vs_x, vs_z; color = set2[1], linewidth = 3.0, label = L"\text{Slow}")
    axislegend(ax; position = :rb, framevisible = false, labelsize = 15)

    dom = DOM.div(
        DOM.p("vₛ / vₐ = ", s_vs),
        fig,
    )
    return Bonito.record_states(session, dom)
end
```

## Interactive 3D Shear Alfvén Wave (Drag to Rotate)

```@example mhd_3d
using WGLMakie, PlasmaWaves
WGLMakie.activate!()
animate_alfven_wave(nothing).fig
```

## Wave Propagation Animations

### Shear Alfvén Wave

```@raw html
<video src="assets/gallery/alfven_wave.mp4" autoplay loop muted playsinline controls style="max-width:100%; border-radius:8px;"></video>
```

### Fast vs. Slow Magnetosonic Waves

```@raw html
<video src="assets/gallery/magnetosonic_waves.mp4" autoplay loop muted playsinline controls style="max-width:100%; border-radius:8px;"></video>
```

### Non-Uniform SAW Phase Mixing

In a non-uniform slab $v_A(x)$, the perpendicular wavenumber grows linearly in time ($|k_x| \simeq |\omega_A'(x)| t$):

```@raw html
<video src="assets/gallery/phase_mixing.mp4" autoplay loop muted playsinline controls style="max-width:100%; border-radius:8px;"></video>
```
