module PlasmaWaves

using LinearAlgebra
using StaticArrays
using SpecialFunctions
import PlasmaPopulations
import PlasmaDispersionBase
import PlasmaBO
using PlasmaPopulations: Plasma, Population, Proton, Electron, BiMaxwellian, BiKappa, Cold
using PlasmaDispersionBase: Wavevector, AngleSweep, CartesianSweep, DispersionProblem, NormalizedPlasma
using PlasmaBO: BOHH, BOPBK, BOFluid

# Re-export key ecosystem types for seamless interoperability
export Plasma, Population, Proton, Electron, BiMaxwellian, BiKappa, Cold
export Wavevector, AngleSweep, CartesianSweep, DispersionProblem, NormalizedPlasma
export BOHH, BOPBK, BOFluid

# --- Core types & constants ---
export Species, PlasmaParams, WaveSolution
export normalized_plasma, si_plasma, to_bo_plasma, get_electron, get_ion
export plasma_frequency, cyclotron_frequency, mass_density, thermal_pressure
export alfven_speed, sound_speed, thermal_speed, thermal_speed_1d
export inertial_length, larmor_radius, debye_length, plasma_beta

# --- Wave diagnostics & field evaluator ---
export phase_speed, damping_rate, propagation_angle
export polarization_E, polarization_B, magnetic_compressibility
export alfven_ratio, electrostatic_ratio, poynting_flux, evaluate_fields

# --- Ideal MHD waves ---
export friedrich_phase_speeds, friedrich_group_velocities
export shear_alfven_wave, circular_alfven_wave
export fast_magnetosonic_wave, slow_magnetosonic_wave
export saw_phase_mixing

# --- Cold plasma & Hall-MHD waves ---
export stix_parameters, cold_dielectric_tensor, cold_cutoffs, cold_resonances
export cold_plasma_n2, cold_plasma_eigenmode, hall_mhd_dispersion
export whistler_wave, ion_cyclotron_wave
export o_mode_wave, x_mode_wave, upper_hybrid_wave, lower_hybrid_wave
export faraday_rotation_fields

# --- Dispersive Alfvén waves (KAW & IAW) ---
export kinetic_alfven_wave, inertial_alfven_wave, kaw_particle_trajectories

# --- Warm & Analytical Kinetic Electrostatic Waves ---
export plasma_dispersion_Z, plasma_dispersion_dZ
export electrostatic_dielectric, langmuir_wave, ion_acoustic_wave, landau_phase_space
export bernstein_dielectric, bernstein_wave

# --- Numerical Dispersion Confirmation (powered by PlasmaBO.jl) ---
export DispersionComparison, validity_mask
export confirm_dispersion, confirm_alfven, confirm_whistler
export confirm_kaw, confirm_langmuir, confirm_ion_acoustic

# --- Symbolic Derivation stubs (implemented in PlasmaWavesSymbolicsExt) ---
export derive_mhd_dispersion, derive_cold_plasma_dispersion, derive_hall_mhd_dispersion
export derive_kaw_iaw_dispersion, derive_bohm_gross_dispersion, symbolic_to_numeric

# --- Visualization stubs (implemented in PlasmaWavesMakieExt) ---
export plot_friedrich, plot_dispersion_parallel, plot_cma_diagram
export plot_kaw_iaw_dispersion, plot_landau_dispersion, plot_bernstein_dispersion
export plot_transport_ratios, plot_confirmation
export animate_alfven_wave, animate_magnetosonic_waves, animate_whistler_emic
export animate_whistler_packet, animate_faraday_rotation, animate_kaw_vs_saw
export animate_hybrid_waves, animate_landau_damping, animate_ion_acoustic
export animate_bernstein_wave, animate_phase_mixing

# --- Bonito web stubs (implemented in PlasmaWavesBonitoExt) ---
export plasma_waves_routes, build_website, serve_website

include("types.jl")
include("mhd.jl")
include("cold_plasma.jl")
include("dispersive_alfven.jl")
include("kinetic.jl")
include("confirm.jl")

# Stub definitions for Symbolics + Latexify extension
function derive_mhd_dispersion end
function derive_cold_plasma_dispersion end
function derive_hall_mhd_dispersion end
function derive_kaw_iaw_dispersion end
function derive_bohm_gross_dispersion end
function symbolic_to_numeric end

# Stub definitions for Makie extension
function plot_friedrich end
function plot_dispersion_parallel end
function plot_cma_diagram end
function plot_kaw_iaw_dispersion end
function plot_landau_dispersion end
function plot_bernstein_dispersion end
function plot_transport_ratios end
function plot_confirmation end

function animate_alfven_wave end
function animate_magnetosonic_waves end
function animate_whistler_emic end
function animate_whistler_packet end
function animate_faraday_rotation end
function animate_kaw_vs_saw end
function animate_hybrid_waves end
function animate_landau_damping end
function animate_ion_acoustic end
function animate_bernstein_wave end
function animate_phase_mixing end

# Stub definitions for Bonito + WGLMakie extension
function plasma_waves_routes end
function build_website end
function serve_website end

end # module PlasmaWaves
