# Theory vs. `PlasmaBO.jl` Kinetic Confirmation

`PlasmaWaves.jl` uses [`PlasmaBO.jl`](https://github.com/JuliaSpacePhysics/PlasmaBO.jl) (`BOHH`) and [`PlasmaDispersionBase.jl`](https://github.com/JuliaPlasma/PlasmaDispersionBase.jl) (`track`) to verify analytical dispersion models against full hot magnetized kinetic solutions, pinpointing exact regime validity limits and kinetic damping rates.

## Ideal MHD Alfvén vs. `PlasmaBO` (Hall Splitting)

Both the left-hand **EMIC** and right-hand **Whistler** branches are tracked backward from $k d_i = 1.2$ into the ideal MHD limit $\omega = k v_A$:

```@example confirm_alfven
using CairoMakie, PlasmaWaves
CairoMakie.activate!()
plot_confirmation(confirm_alfven())
```

## Cold Whistler vs. `PlasmaBO`

```@example confirm_whistler
using CairoMakie, PlasmaWaves
CairoMakie.activate!()
plot_confirmation(confirm_whistler())
```

## Bessel KAW vs. `PlasmaBO`

```@example confirm_kaw
using CairoMakie, PlasmaWaves
CairoMakie.activate!()
plot_confirmation(confirm_kaw())
```

## Bohm–Gross Langmuir vs. `PlasmaBO` (Landau Damping)

```@example confirm_langmuir
using CairoMakie, PlasmaWaves
CairoMakie.activate!()
plot_confirmation(confirm_langmuir())
```
