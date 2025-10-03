var             

// Aggregate non-expectational variables

       invest, invest_flex, y_va, y_va_flex, g, g_flex, r, rr, rr_flex, y_va_nonservice, y_va_service,
       gam_I, gam_I_flex,

             gam_G,             gam_G_flex,
        pie_I,
        
        w, w_flex, 
        pie_W, pie_ye,
        Lpie, L2pie, 
        lab, lab_flex,

        hrs,  hrs_flex,  
        
        c, c_flex, eps_m,
        
        lam,lam_flex,
        trans, gap, 
    lam_nr, lam_nr_flex, pie_core, pie_average, pie_average_NYFed, pie_forecast, pie_headmin, tfp, INFL, FFR, FFRDELTA, YGROWTH, YGAP, VA_MINING, DEFL_MINING, VA_MANUFACTURING, DEFL_MANUFACTURING, INFL_Annual, YGROWTH_Annual, Pie_wr, Lpie_va_2, L2pie_va_2, Ly_va_2, L2y_va_2, L3y_va_2, Lpie_va_5, L2pie_va_5, Ly_va_5, L2y_va_5, L3y_va_5, Ly_va, L2y_va, L3y_va

// Lags for AIT rule        
@#for nait in 0 : 16
        lpie_@{nait}
                        
        @#endfor

// Leads for forecast rule        
@#for nfore in 1 : 4
        fpie_@{nfore}
                        
        @#endfor


@#for se in 1 : nsector

            // Sectoral non-expectational variables
            y_@{se}, y_flex_@{se}, 
            y_va_@{se}, y_va_flex_@{se},
            nva_@{se}, nva_flex_@{se},
            f_@{se}, f_flex_@{se},
            x_@{se}, x_flex_@{se}, 
            mc_@{se},
            c_@{se}, c_flex_@{se},
            i_@{se}, i_flex_@{se},
            g_@{se}, g_flex_@{se},
            k_@{se}, k_flex_@{se},
            ks_@{se}, ks_flex_@{se}, 
            capu_@{se}, capu_flex_@{se},
            n_@{se}, n_flex_@{se},
            w_@{se}, w_flex_@{se}, 
            gam_@{se}, gam_flex_@{se}, 
            gam_f_@{se}, gam_f_flex_@{se}, 
            gam_x_@{se}, gam_x_flex_@{se},
            gam_va_@{se},      gam_va_flex_@{se} ,
            a_@{se}, om_c_@{se}, pie_va_@{se}, mu_@{se},
                
               
                       

          @#for se_2 in 1 : nsector
               xx_@{se_2}_@{se},
                xx_flex_@{se_2}_@{se},
            @#endfor
        @#endfor
        

        // Aggregate expectational variables
        c_r, c_r_flex,
        c_nr, c_nr_flex, 
        pie,

        lam_r,lam_r_flex,
        cshift 
         
                 
        // Sectoral expectational variables

@#for se in 1 : nsector
       ,     q_@{se},  q_flex_@{se}, 
            rk_@{se},  rk_flex_@{se},
            z_@{se}, z_flex_@{se},
           

        pie_W_@{se},       pie_@{se}
                        
        @#endfor

        

        ;



