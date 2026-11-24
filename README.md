# PlasmaWaves.jl

A Julia package combining **analytical plasma wave theory**, **symbolic mathematical derivations** (`Symbolics.jl` + `Latexify.jl`), **kinetic numerical dispersion confirmation** (`PlasmaBO.jl` + `PlasmaDispersionBase.jl` + `PlasmaPopulations.jl`), and **interactive 2D/3D visualizations & web demos** (`Makie.jl` + `Bonito.jl`).

Rather than reinventing another kinetic dispersion solver, `PlasmaWaves.jl` focuses on **understanding wave physics**: deriving and analyzing eigenmodes, polarizations, and transport ratios, confirming where textbook fluid/asymptotic theories hold (and where kinetic damping, Hall splitting, or FLR effects cause them to break down) against `PlasmaBO.jl`, and illustrating wave propagation with publication-grade animations.

## Supported Wave Families

| Regime | Analytical Theory, Confirmation & Animations |
| :--- | :--- |
| **Ideal MHD** | Shear Alfvén (SAW), Circularly Polarized Alfvén (CPAW), Fast & Slow Magnetosonic, Friedrich phase & group diagrams, Non-uniform SAW phase mixing |
| **Cold Plasma & Hall-MHD** | Stix dielectric tensor $(R, L, P, S, D)$, Whistler / Electron Cyclotron, EMIC / Ion Cyclotron, O/X modes, Upper & Lower Hybrid, Faraday rotation, CMA diagram |
| **Dispersive Alfvén** | Kinetic Alfvén Wave (KAW: Bessel & two-fluid Padé), Inertial Alfvén Wave (IAW), parallel $E_\parallel$ & resonant electron trapping |
| **Kinetic / Vlasov** | Plasma dispersion function $Z(\zeta)$, Langmuir / Bohm-Gross + collisionless Landau damping & $(x, v)$ phase space, Ion-Acoustic wave ($T_e/T_i$ damping), Electron & Ion Bernstein Waves (EBW / IBW) |
| **Numerical Confirmation** | Direct branch comparison of any analytical or symbolically derived dispersion against **`PlasmaBO.jl`** (`BOHH`, `BOPBK`, `BOFluid`) |
| **Symbolic Derivations** | Symbolic dispersion relations, group velocities, and LaTeX export for Ideal MHD, Stix CPDR, Hall-MHD, KAW/IAW, and Bohm-Gross |

## Quick Start

### 1. Analytical Wave Theory & `PlasmaBO.jl` Confirmation

```julia
using PlasmaWaves

# Ideal MHD Alfvén vs. full kinetic PlasmaBO (shows Hall splitting at k d_i ~ 1)
comp_alf  = confirm_alfven()
comp_wh   = confirm_whistler()
comp_lang = confirm_langmuir()

# Or pass a PlasmaPopulations.Plasma directly
pp = Plasma(
    Population(Proton, BiMaxwellian; T = 10.0, n = 5e6),
    Population(Electron, BiMaxwellian; T = 10.0, n = 5e6);
    B0 = 1e-8,
)
p = PlasmaParams(pp)
fast = fast_magnetosonic_wave(p, 1e-5, π/3)
```

### 2. Symbolic Math Derivations (`Symbolics.jl` + `Latexify.jl`)

```julia
using PlasmaWaves, Symbolics, Latexify

d_mhd  = derive_mhd_dispersion()       # 3×3 ideal MHD determinant -> Alfvén + Fast/Slow
d_cold = derive_cold_plasma_dispersion() # Stix A n⁴ - B n² + C = 0 + θ=0, π/2 limits
d_hall = derive_hall_mhd_dispersion()  # Whistler/EMIC + symbolic group velocities
d_kaw  = derive_kaw_iaw_dispersion()   # Two-fluid KAW & IAW + E_∥/E_⊥
d_bg   = derive_bohm_gross_dispersion()

println(d_mhd.latex) # Ready to paste into Quarto / KeyNotes!
```

### 3. Animations & Confirmation Figures (`CairoMakie` / `GLMakie`)

```julia
using CairoMakie, PlasmaWaves

save("confirm_alfven.png", plot_confirmation(confirm_alfven()))
save("friedrich.png", plot_friedrich())
animate_alfven_wave("alfven_wave.gif")
animate_magnetosonic_waves("magnetosonic_waves.gif")
animate_kaw_vs_saw("kaw_vs_saw.gif")
animate_landau_damping("landau_damping.gif")
```

To render the entire gallery of figures and GIFs:
```bash
julia --project=examples examples/generate_gallery.jl gallery
```

### 4. Interactive Website (`Bonito.jl` + `WGLMakie`)

```bash
# Live interactive server
julia --project=web web/serve.jl

# Export static multi-page HTML site (including Theory vs. PlasmaBO confirmation page)
julia --project=web web/build_site.jl site
```
