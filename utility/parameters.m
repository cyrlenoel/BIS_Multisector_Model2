parameters 
    @#for se in 1:nsector
        omega_c_@{se}, rho_om_c_@{se}, sig_om_c_@{se}, omega_g_@{se}, omega_ns_@{se}, 
        omega_i_@{se}, omega_nd_@{se}, omega_y_@{se}, tau_c_@{se}, 
        rho_a_@{se}, rho_mu_@{se}, sig_a_@{se}, sig_mu_@{se}, thet_@{se}, chi_p_@{se},  
        
        ss_gam_@{se},ss_gam_f_@{se}, ss_gam_x_@{se},ss_w_@{se}, ss_rk_@{se}, ss_nva_@{se}, ss_n_@{se}, ss_y_@{se}, ss_x_@{se},
        ss_c_@{se},ss_i_@{se},ss_g_@{se},ss_a_@{se}, ss_mu_@{se}, ss_z_@{se},  
        ss_y_va_@{se}, ss_f_@{se}, ss_k_@{se}, ss_ks_@{se}, c_mu_mm_@{se}, c_tfp_mm_@{se},

        @#for se_2 in 1 : nsector
            ss_xx_@{se_2}_@{se},
              ii_weights_@{se_2}_@{se},
        @#endfor

    @#endfor

    eta, eta_g, zeta, vphi, psi, xi, bet, h, capup, eps_w, nu, delt, Spp, chi_w, rho_r, phi_pie, phi_pie_NYFed, phi_pie_core, phi_pie_bar, phi_pie_forecast, phi_pie_headmin, phi_y, phi_gap, sig_r, rho_g, sig_g, sig_cshift,
    A_N, thet_w, rk_wedge, omega_r, rho_trans, sig_trans, rho_tfp, sig_tfp, sig_mm, 

    ss_gam_I, ss_gam_G, ss_w, ss_n, ss_trans, ss_c_r, ss_c, ss_c_nr, ss_g, 
    ss_lam_r, ss_lam, ss_lam_nr, ss_lab,  ss_invest, ss_y_va, ss_hrs, ss_r, ss_rr, n_ait
    // , ss_cshift, rho_cshift, sig_cshift
    ;