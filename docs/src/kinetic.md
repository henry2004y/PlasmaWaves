# Kinetic, Landau & Bernstein Waves

## Unmagnetized Vlasov Electrostatic Dispersion

Using the Fried–Conte plasma dispersion function $Z(\zeta) = i\sqrt{\pi}\,\mathrm{erfcx}(-i\zeta)$:
$$\epsilon(\omega, k) = 1 + \sum_s \frac{1}{k^2 \lambda_{Ds}^2}\left[1 + \zeta_s Z(\zeta_s)\right] = 0, \qquad \zeta_s = \frac{\omega}{|k| v_{th,s}}$$

- **Langmuir Waves**: follow Bohm–Gross ($\omega_r^2 \approx \omega_{pe}^2 + 3 k^2 v_{th,e}^2$) at small $k\lambda_{De}$, then suffer strong collisionless **Landau damping** as $k\lambda_{De} \to 0.5$.
- **Ion-Acoustic Waves**: propagate with weak damping when $T_e \gg T_i$ and become heavily ion-Landau-damped as $T_e / T_i \to 1$.

## Electron Bernstein Waves ($k_\parallel = 0$)

Perpendicular electrostatic cyclotron harmonic passbands propagate undamped between consecutive harmonics $(n\Omega_e, (n+1)\Omega_e)$ when $k_\perp r_{Le} \sim \mathcal{O}(1)$:

```@example kin_disp
using CairoMakie, PlasmaWaves
CairoMakie.activate!()
plot_landau_dispersion()
```

```@example kin_ebw
using CairoMakie, PlasmaWaves
CairoMakie.activate!()
plot_bernstein_dispersion()
```

## Kinetic Wave Animations

### Langmuir Wave & $(x, v)$ Phase-Space Landau Damping

```@raw html
<video src="assets/gallery/landau_damping.mp4" autoplay loop muted playsinline controls style="max-width:100%; border-radius:8px;"></video>
```

### Ion-Acoustic Wave Damping ($T_e / T_i = 15$ vs. $1$)

```@raw html
<video src="assets/gallery/ion_acoustic.mp4" autoplay loop muted playsinline controls style="max-width:100%; border-radius:8px;"></video>
```

### Electron Bernstein Wave

```@raw html
<video src="assets/gallery/bernstein_wave.mp4" autoplay loop muted playsinline controls style="max-width:100%; border-radius:8px;"></video>
```
