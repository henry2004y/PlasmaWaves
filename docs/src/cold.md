# Cold Plasma, Hall-MHD & Whistler Waves

In the Stix cold-plasma formulation ($T_e = T_i = 0$), the magnetized dielectric tensor $\boldsymbol{\epsilon}_R(\omega)$ depends on $(R, L, P, S, D)$:
$$\tan^2\theta = \frac{-P(n^2 - R)(n^2 - L)}{(S n^2 - RL)(n^2 - P)}$$

In the Hall-MHD regime ($\omega \ll \Omega_e$, $k_\parallel d_i \sim 1$), the ideal MHD Alfvén degeneracy splits into:
$$\frac{\omega}{k_\parallel v_A} = \sqrt{1 + \left(\frac{k_\parallel d_i}{2}\right)^2} \pm \frac{k_\parallel d_i}{2}$$
where `$+$` is the **right-hand Whistler** branch ($\omega > k_\parallel v_A$) and `$-$` is the **left-hand EMIC / ion cyclotron** branch ($\omega < k_\parallel v_A$).

## Parallel Dispersion & CMA Diagram

```@example cold_disp
using CairoMakie, PlasmaWaves
CairoMakie.activate!()
plot_dispersion_parallel()
```

```@example cold_cma
using CairoMakie, PlasmaWaves
CairoMakie.activate!()
plot_cma_diagram()
```

## Wave Propagation Animations

### Whistler vs. EMIC Helical Polarization

```@raw html
<video src="assets/gallery/whistler_emic.mp4" autoplay loop muted playsinline controls style="max-width:100%; border-radius:8px;"></video>
```

### Dispersive Whistler Wave Packet

High-frequency components have larger group velocity $v_g \propto \sqrt{\omega}$ and outrun low-frequency components:

```@raw html
<video src="assets/gallery/whistler_packet.mp4" autoplay loop muted playsinline controls style="max-width:100%; border-radius:8px;"></video>
```

### Faraday Rotation

```@raw html
<video src="assets/gallery/faraday_rotation.mp4" autoplay loop muted playsinline controls style="max-width:100%; border-radius:8px;"></video>
```

### Lower Hybrid Wave

```@raw html
<video src="assets/gallery/hybrid_waves.mp4" autoplay loop muted playsinline controls style="max-width:100%; border-radius:8px;"></video>
```
