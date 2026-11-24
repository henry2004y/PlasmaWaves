using Test
using LinearAlgebra
using StaticArrays
using PlasmaWaves

@testset "PlasmaWaves.jl" begin
    @testset "Core Plasma Parameters & PlasmaPopulations Interop" begin
        p = normalized_plasma(; mi_me = 1836.0, β_i = 0.2, β_e = 0.3, vA_c = 0.01)
        @test isapprox(alfven_speed(p), 1.0; rtol = 1e-3)
        @test isapprox(abs(cyclotron_frequency(get_ion(p), p)), 1.0; rtol = 1e-12)
        @test isapprox(cyclotron_frequency(get_electron(p), p), -1836.0; rtol = 1e-12)
        @test isapprox(plasma_beta(p), 0.5; rtol = 1e-12)

        p_si = si_plasma(; B0 = 1e-8, n0 = 5e6, Ti_eV = 10.0, Te_eV = 10.0)
        @test alfven_speed(p_si) > 5e4 && alfven_speed(p_si) < 2e5
        @test plasma_frequency(get_electron(p_si), p_si) > abs(cyclotron_frequency(get_electron(p_si), p_si))

        # Bidirectional adapter with PlasmaPopulations.Plasma & NormalizedPlasma
        norm_p, u = to_bo_plasma(p)
        @test norm_p isa NormalizedPlasma
        @test isapprox(u.omega, 1.0; rtol = 1e-12)

        pp = Plasma(
            Population(Proton, BiMaxwellian; T = 10.0, n = 5e6),
            Population(Electron, BiMaxwellian; T = 10.0, n = 5e6);
            B0 = 1e-8,
        )
        p_from_pp = PlasmaParams(pp)
        @test isapprox(p_from_pp.B0, 1e-8; rtol = 1e-12)
        @test length(p_from_pp.species) == 2
    end

    @testset "Ideal MHD Waves" begin
        p_low = normalized_plasma(; β_i = 0.1, β_e = 0.1)
        p_high = normalized_plasma(; β_i = 1.5, β_e = 1.5)
        vA = alfven_speed(p_low)
        vs = sound_speed(p_low)
        @test vs < vA

        saw = shear_alfven_wave(p_low, 1.5, π / 4)
        @test isapprox(real(saw.ω), 1.5 * cos(π / 4) * vA; rtol = 1e-12)
        @test isapprox(magnetic_compressibility(saw), 0.0; atol = 1e-12)
        @test isapprox(alfven_ratio(saw), 1.0; rtol = 1e-12)
        @test isapprox(abs(saw.δE[1] / saw.δB[2]), vA; rtol = 1e-12)
        @test isapprox(dot(saw.k, saw.δB), 0.0; atol = 1e-12)

        cpaw_r = circular_alfven_wave(p_low, 1.0; hand = :right)
        cpaw_l = circular_alfven_wave(p_low, 1.0; hand = :left)
        @test isapprox(real(polarization_E(cpaw_r)), 1.0; rtol = 1e-12)
        @test isapprox(real(polarization_E(cpaw_l)), -1.0; rtol = 1e-12)
        f0 = evaluate_fields(cpaw_r, SVector(0.0, 0.0, 0.0), 0.0)
        f1 = evaluate_fields(cpaw_r, SVector(0.0, 0.0, 0.7), 1.3)
        @test isapprox(norm(f0.B), norm(f1.B); rtol = 1e-12)

        fast = fast_magnetosonic_wave(p_low, 1.0, π / 3)
        slow = slow_magnetosonic_wave(p_low, 1.0, π / 3)
        @test phase_speed(fast) > phase_speed(slow)
        @test isapprox(dot(fast.k, fast.δB), 0.0; atol = 1e-12)
        @test isapprox(dot(slow.k, slow.δB), 0.0; atol = 1e-12)
        @test real(fast.δp * conj(fast.δB[3])) > 0
        @test real(slow.δp * conj(slow.δB[3])) < 0

        sp_h = friedrich_phase_speeds(alfven_speed(p_high), sound_speed(p_high), 0.0)
        @test isapprox(sp_h.v_fast, sound_speed(p_high); rtol = 1e-12)
        @test isapprox(sp_h.v_slow, alfven_speed(p_high); rtol = 1e-12)

        xs = range(-2.0, 2.0, length = 50)
        by0 = saw_phase_mixing(xs, 0.0)
        @test isapprox(by0[25], exp(-xs[25]^2 / 4); rtol = 1e-12)
    end

    @testset "Cold Plasma & Hall-MHD" begin
        p = normalized_plasma(; mi_me = 100.0, vA_c = 0.1)
        st = stix_parameters(p, 5.0)
        @test isapprox(st.S^2 - st.D^2, st.R * st.L; rtol = 1e-12)

        cut = cold_cutoffs(p)
        res = cold_resonances(p)
        @test abs(stix_parameters(p, cut.ω_pe).P) < 1e-10
        @test abs(stix_parameters(p, cut.ω_R).R) < 1e-10
        @test abs(stix_parameters(p, cut.ω_L).L) < 1e-10
        @test abs(stix_parameters(p, res.ω_UH).S) < 1e-10
        @test abs(stix_parameters(p, res.ω_LH).S) < 1e-10

        wh = whistler_wave(p, 1.5; model = :hall)
        em = ion_cyclotron_wave(p, 1.5; model = :hall)
        @test phase_speed(wh) > alfven_speed(p)
        @test phase_speed(em) < alfven_speed(p)
        @test real(polarization_E(wh)) > 0
        @test real(polarization_E(em)) < 0

        wh_obl = whistler_wave(p, 1.5, π / 6; model = :cold)
        @test real(wh_obl.ω) < abs(cyclotron_frequency(get_electron(p), p)) * cos(π / 6)

        om = o_mode_wave(p, 2.0)
        @test real(om.ω) > cut.ω_pe
        uh = upper_hybrid_wave(p, 0.1)
        @test real(uh.ω) >= res.ω_UH
        lh = lower_hybrid_wave(p, 0.1)
        @test isapprox(real(lh.ω), res.ω_LH; rtol = 1e-12)
        fr = faraday_rotation_fields(p, 1.5 * cut.ω_R, [0.0, 1.0], 0.0)
        @test fr.rotation_rate != 0.0
    end

    @testset "Dispersive Alfvén Waves (KAW & IAW)" begin
        p = normalized_plasma(; β_i = 0.2, β_e = 0.2)
        riL = larmor_radius(get_ion(p), p)
        vA = alfven_speed(p)

        kaw_small = kinetic_alfven_wave(p, 0.1, 1e-4 / riL; model = :bessel)
        @test isapprox(real(kaw_small.ω), 0.1 * vA; rtol = 1e-4)

        kaw_b = kinetic_alfven_wave(p, 0.1, 1.0 / riL; model = :bessel)
        kaw_2 = kinetic_alfven_wave(p, 0.1, 1.0 / riL; model = :twofluid)
        @test real(kaw_b.ω) > 0.1 * vA
        @test isapprox(real(kaw_b.ω), real(kaw_2.ω); rtol = 0.15)
        @test abs(kaw_b.δE[3]) > 0

        p_cold = normalized_plasma(; β_i = 1e-5, β_e = 1e-5)
        de = inertial_length(get_electron(p_cold), p_cold)
        iaw = inertial_alfven_wave(p_cold, 0.1, 1.0 / de)
        @test isapprox(real(iaw.ω), 0.1 * alfven_speed(p_cold) / sqrt(2.0); rtol = 1e-10)
        @test abs(iaw.δE[3]) > 0

        ts = range(0, 1.0, length = 20)
        traj = kaw_particle_trajectories(kaw_b, [0.0, 1.0], [1.0, 1.2], ts)
        @test size(traj.z) == (2, 20)
    end

    @testset "Analytical Kinetic Electrostatic & Bernstein Waves" begin
        @test isapprox(plasma_dispersion_Z(0.0), im * sqrt(π); rtol = 1e-14)
        ζ = 2.5 + 0.3im
        @test isapprox(plasma_dispersion_Z(-conj(ζ)), -conj(plasma_dispersion_Z(ζ)); rtol = 1e-13)
        @test isapprox(real(plasma_dispersion_Z(10.0)), -0.1 * (1 + 1 / 200 + 3 / 40_000 + 15 / 8_000_000); rtol = 1e-7)

        p = normalized_plasma(; β_i = 0.05, β_e = 0.5)
        e = get_electron(p)
        λDe = debye_length(e, p)
        s_lang = langmuir_wave(p, 0.2 / λDe; model = :vlasov)
        s_bg = langmuir_wave(p, 0.2 / λDe; model = :bohm_gross)
        @test isapprox(real(s_lang.ω), real(s_bg.ω); rtol = 0.05)
        @test imag(s_lang.ω) < 0
        @test abs(electrostatic_dielectric(p, s_lang.ω, 0.2 / λDe; species_list = [e])) < 1e-7

        p_hot = normalized_plasma(; β_i = 0.02, β_e = 0.40)
        p_eq  = normalized_plasma(; β_i = 0.20, β_e = 0.20)
        ia_hot = ion_acoustic_wave(p_hot, 1.0; model = :vlasov)
        ia_eq  = ion_acoustic_wave(p_eq, 1.0; model = :vlasov)
        @test imag(ia_hot.ω) < 0 && imag(ia_eq.ω) < 0
        @test abs(imag(ia_eq.ω) / real(ia_eq.ω)) > abs(imag(ia_hot.ω) / real(ia_hot.ω))

        p_ebw = normalized_plasma(; mi_me = 1836.0, β_i = 0.05, β_e = 0.05, vA_c = 0.05)
        Ωe = abs(cyclotron_frequency(get_electron(p_ebw), p_ebw))
        rLe = larmor_radius(get_electron(p_ebw), p_ebw)
        ebw1 = bernstein_wave(p_ebw, 1.5 / rLe; species = :electron, harmonic = 1)
        ebw2 = bernstein_wave(p_ebw, 1.5 / rLe; species = :electron, harmonic = 2)
        @test Ωe < real(ebw1.ω) < 2 * Ωe
        @test 2 * Ωe < real(ebw2.ω) < 3 * Ωe
    end

    @testset "PlasmaBO.jl Numerical Dispersion Confirmation" begin
        # 1. Ideal MHD Alfvén agrees at small k d_i, diverges at k d_i ~ 1 (Hall splitting)
        comp_alf = confirm_alfven(
            normalized_plasma(; mi_me = 100.0, β_i = 0.05, β_e = 0.05, vA_c = 0.05),
            [0.05, 0.1, 0.8];
            alg = BOHH(N = 3, J = 6),
        )
        @test comp_alf.rel_error[1] < 0.05
        @test comp_alf.rel_error[end] > comp_alf.rel_error[1]
        @test validity_mask(comp_alf; tol = 0.05)[1]

        # 2. Cold Whistler agrees with PlasmaBO in the whistler band
        comp_wh = confirm_whistler(
            normalized_plasma(; mi_me = 100.0, β_i = 0.02, β_e = 0.02, vA_c = 0.05),
            [0.5, 1.5];
            alg = BOHH(N = 3, J = 6),
        )
        @test comp_wh.rel_error[1] < 0.05
        @test comp_wh.rel_error[2] < 0.05

        # 3. Bohm-Gross Langmuir agrees at small kλ_D and PlasmaBO captures Landau damping
        comp_lang = confirm_langmuir(
            normalized_plasma(; mi_me = 100.0, β_i = 0.05, β_e = 0.5, vA_c = 0.05),
            [0.12, 0.35];
            alg = BOHH(N = 2, J = 8),
        )
        @test comp_lang.rel_error[1] < 0.05
        @test imag(comp_lang.ω_numerical[2]) < 0 # kinetic Landau damping!
    end

    @testset "Symbolics + Latexify Wave Theory Derivations" begin
        using Symbolics, Latexify

        d_mhd = derive_mhd_dispersion()
        @test d_mhd.name === :ideal_mhd
        @test !isempty(d_mhd.latex)
        # Numerical spot check: symbolic ω_fast matches friedrich_phase_speeds
        f_fast = symbolic_to_numeric(
            d_mhd.solutions.ω_fast,
            d_mhd.variables.k,
            Dict(d_mhd.variables.v_A => 1.0, d_mhd.variables.v_s => 0.6, d_mhd.variables.θ => π / 3),
        )
        @test isapprox(f_fast(2.0), 2.0 * friedrich_phase_speeds(1.0, 0.6, π / 3).v_fast; rtol = 1e-12)

        d_cold = derive_cold_plasma_dispersion()
        @test d_cold.name === :cold_plasma
        @test !isempty(d_cold.latex)

        d_hall = derive_hall_mhd_dispersion()
        @test d_hall.name === :hall_mhd
        f_wh = symbolic_to_numeric(
            d_hall.solutions.ω_whistler,
            d_hall.variables.k_par,
            Dict(d_hall.variables.v_A => 1.0, d_hall.variables.d_i => 1.0),
        )
        p_h = normalized_plasma(; mi_me = 10000.0)
        @test isapprox(f_wh(1.2), hall_mhd_dispersion(p_h, 1.2; branch = :whistler); rtol = 1e-3)

        d_kaw = derive_kaw_iaw_dispersion()
        @test d_kaw.name === :kaw_iaw

        d_bg = derive_bohm_gross_dispersion()
        @test d_bg.name === :bohm_gross
        # Plug symbolic Bohm-Gross directly into PlasmaBO confirmation!
        p_l = normalized_plasma(; mi_me = 100.0, β_i = 0.05, β_e = 0.5, vA_c = 0.05)
        e_l = get_electron(p_l)
        f_bg = symbolic_to_numeric(
            d_bg.solutions.ω_bg,
            d_bg.variables.k,
            Dict(d_bg.variables.ω_pe => plasma_frequency(e_l, p_l), d_bg.variables.v_th => thermal_speed_1d(e_l, p_l)),
        )
        λDe = debye_length(e_l, p_l)
        comp_sym = confirm_dispersion(
            k -> f_bg(k),
            p_l,
            [0.12 / λDe];
            alg = BOHH(N = 2, J = 8),
            Ω_ref = plasma_frequency(e_l, p_l),
        )
        @test comp_sym.rel_error[1] < 0.05
    end

    @testset "Makie Visualization Extension Smoke Tests" begin
        using CairoMakie
        mktempdir() do dir
            fig = plot_friedrich()
            out_png = joinpath(dir, "friedrich.png")
            save(out_png, fig)
            @test isfile(out_png) && filesize(out_png) > 1000

            comp = confirm_alfven(
                normalized_plasma(; mi_me = 100.0, β_i = 0.05, β_e = 0.05, vA_c = 0.05),
                [0.05, 0.2, 0.6];
                alg = BOHH(N = 3, J = 6),
            )
            out_conf = joinpath(dir, "confirm.png")
            save(out_conf, plot_confirmation(comp))
            @test isfile(out_conf) && filesize(out_conf) > 1000

            out_gif = joinpath(dir, "alfven.gif")
            animate_alfven_wave(out_gif; nframes = 6, fps = 6)
            @test isfile(out_gif) && filesize(out_gif) > 1000
        end
    end
end
