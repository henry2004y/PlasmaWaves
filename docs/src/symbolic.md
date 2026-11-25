# Symbolic Wave Theory Derivations

When `Symbolics.jl` and `Latexify.jl` are loaded, `PlasmaWavesSymbolicsExt` activates symbolic derivation routines that construct linearized wave matrices, derive closed-form dispersion relations & group velocities, export LaTeX strings for lecture notes, and compile symbolic expressions directly into callables for `confirm_dispersion`.

## Ideal MHD Dispersion Matrix

```@example sym_mhd
using Symbolics, Latexify, PlasmaWaves
d_mhd = derive_mhd_dispersion()
println("LaTeX: ", d_mhd.latex)
d_mhd
```

## Cold-Plasma Stix Dispersion & Principal Limits

```@example sym_cold
using Symbolics, Latexify, PlasmaWaves
d_cold = derive_cold_plasma_dispersion()
println("Parallel (θ=0):      ", d_cold.solutions.parallel_limit)
println("Perpendicular (π/2): ", d_cold.solutions.perpendicular_limit)
d_cold
```

## Hall-MHD & KAW/IAW Symbolic Group Velocities

```@example sym_hall_kaw
using Symbolics, Latexify, PlasmaWaves
d_hall = derive_hall_mhd_dispersion()
d_kaw  = derive_kaw_iaw_dispersion()
println("Whistler v_g = ", d_hall.solutions.vg_whistler)
println("KAW E_z/E_x  = ", d_kaw.solutions.ez_ex_kaw)
```

## Symbolic Bohm–Gross Compiled into `PlasmaBO` Confirmation

```@example sym_confirm
using Symbolics, Latexify, CairoMakie, PlasmaWaves
CairoMakie.activate!()

d_bg = derive_bohm_gross_dispersion()
p = normalized_plasma(; mi_me = 100.0, β_i = 0.05, β_e = 0.5, vA_c = 0.05)
e = get_electron(p)
ωpe = plasma_frequency(e, p)
λDe = debye_length(e, p)

f_bg = symbolic_to_numeric(
    d_bg.solutions.ω_bg,
    d_bg.variables.k,
    Dict(d_bg.variables.ω_pe => ωpe, d_bg.variables.v_th => thermal_speed_1d(e, p)),
)

kλs = range(0.08, 0.45, length = 14)
comp = confirm_dispersion(
    kλ -> f_bg(kλ / λDe),
    p,
    kλs;
    k_vectors = [Wavevector(0.0, Float64(kλ / λDe)) for kλ in kλs],
    k_label = "\$k \\lambda_{De}\$",
    ω_scale = ωpe,
    ω_label = "\$\\omega / \\omega_{pe}\$",
    Ω_ref = ωpe,
    mode = :symbolic_bohm_gross,
    primary_label = "PlasmaBO Kinetic",
)
plot_confirmation(comp)
```
