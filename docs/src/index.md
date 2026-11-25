# PlasmaWaves.jl

**PlasmaWaves.jl** combines **analytical plasma wave theory**, **symbolic mathematical derivations** (`Symbolics.jl` + `Latexify.jl`), **kinetic numerical dispersion confirmation** (`PlasmaBO.jl` + `PlasmaDispersionBase.jl` + `PlasmaPopulations.jl`), and **Makie.jl** figures, **MP4** wave animations, and local **Bonito.jl** interactive apps.

!!! tip "Companion Lecture Notes"
    For full step-by-step physical derivations, historical context, and observational discussions, see the companion [KeyNotes Plasma Physics Chapters](https://henry2004y.github.io/KeyNotes/). This site serves as the interactive computational & visual laboratory.

## Wave Regimes Covered

| Regime | Wave Modes & Physics | Page |
| :--- | :--- | :--- |
| **Ideal MHD** | Shear Alfvén (SAW), Circularly Polarized Alfvén (CPAW), Fast & Slow Magnetosonic, Friedrich phase & group diagrams, SAW phase mixing | [MHD & Alfvén Waves](mhd.md) |
| **Cold Plasma & Hall-MHD** | Stix dielectric tensor $(R, L, P, S, D)$, Whistler / EMIC, O/X modes, Upper & Lower Hybrid, Faraday rotation, CMA diagram | [Cold Plasma & Whistlers](cold.md) |
| **Dispersive Alfvén** | Kinetic Alfvén Wave (KAW: Bessel & two-fluid Padé), Inertial Alfvén Wave (IAW), finite $E_\parallel$ & electron trapping | [Kinetic & Inertial Alfvén](kaw.md) |
| **Kinetic / Vlasov** | Plasma dispersion function $Z(\zeta)$, Langmuir / Bohm-Gross + Landau damping, Ion-Acoustic wave ($T_e/T_i$ damping), Bernstein waves | [Kinetic, Landau & Bernstein](kinetic.md) |
| **Kinetic Confirmation** | Tracked branch comparisons of analytical dispersion against `PlasmaBO.jl` (`BOHH`) | [Theory vs. PlasmaBO](confirm.md) |
| **Symbolic Derivations** | Symbolic determinants, group velocities, and LaTeX export | [Symbolic Derivations](symbolic.md) |

## Quick Start

```julia
using PlasmaWaves

# Create a normalized proton-electron plasma (v_A = 1, Ω_i = 1, d_i = 1)
p = normalized_plasma(; β_i = 0.2, β_e = 0.2)

# Ideal MHD & Transport Ratios
fast = fast_magnetosonic_wave(p, 1.0, π / 3)
slow = slow_magnetosonic_wave(p, 1.0, π / 3)
@show magnetic_compressibility(fast) alfven_ratio(slow)

# Confirm textbook theory against full kinetic PlasmaBO.jl
comp = confirm_alfven()
```
