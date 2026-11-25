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

## Friedrich Diagram Regime Sweep

The animation below sweeps $v_s / v_A$ smoothly between the low-$\beta$ ($v_s < v_A$) and high-$\beta$ ($v_s > v_A$) regimes, illustrating how the fast and slow magnetosonic phase-speed surfaces deform relative to the shear Alfvén cusp:

```@raw html
<video src="assets/gallery/friedrich_sweep.mp4" autoplay loop muted playsinline controls style="max-width:100%; border-radius:8px;"></video>
```

!!! tip "Local Interactive Server"
    For live sliders and rotatable 3D `WGLMakie` scenes locally in your browser, run `julia --project=web web/serve.jl`.

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
