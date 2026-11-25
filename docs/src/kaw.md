# Kinetic & Inertial Alfvén Waves

When the perpendicular scale approaches the ion Larmor radius $r_{iL}$ or electron skin depth $d_e$, two distinct dispersive shear Alfvén branches emerge depending on the electron $\beta_e$:

- **Kinetic Alfvén Wave (KAW, $\beta_e \gg m_e / m_i$)**:
  $$\omega^2 = k_\parallel^2 v_A^2 \left[\frac{b_i}{1 - I_0(b_i)e^{-b_i}} + \frac{T_e}{T_i} b_i\right], \qquad b_i = k_\perp^2 r_{iL}^2$$
- **Inertial Alfvén Wave (IAW, $\beta_e \ll m_e / m_i$)**:
  $$\omega^2 = \frac{k_\parallel^2 v_A^2}{1 + k_\perp^2 d_e^2}$$

Unlike ideal MHD Shear Alfvén waves ($E_\parallel = 0$), both KAW and IAW support a finite **parallel electric field** $E_\parallel$ capable of accelerating and trapping electrons along $\mathbf{B}_0$.

## Dispersion & Field Ratios

```@example kaw_disp
using CairoMakie, PlasmaWaves
CairoMakie.activate!()
plot_kaw_iaw_dispersion()
```

## KAW vs. Ideal SAW Electron Trapping

```@raw html
<video src="assets/gallery/kaw_vs_saw.mp4" autoplay loop muted playsinline controls style="max-width:100%; border-radius:8px;"></video>
```
