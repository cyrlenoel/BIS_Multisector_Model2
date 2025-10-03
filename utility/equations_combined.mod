model(linear);

/// Equations 95-104

/// Capital accumulation (95)
// CHECK OK
@#for se in 1:nsector
[name = 'capital accumulation - industry', type = 'sticky']
k_@{se}      - delt*z_@{se}      = (1-delt) * k_@{se}(-1);
@#endfor

@#for se in 1:nsector
[name = 'capital accumulation - industry - flex', type = 'flex']
k_flex_@{se} - delt*z_flex_@{se} = (1-delt) * k_flex_@{se}(-1);
@#endfor

/// Consumer price index (96)
// CHECK OK
[name = 'consumer price index', type = 'sticky']
0 = 
@#for se in 1:nsector
     + (omega_c_@{se}*ss_gam_@{se}^(1-eta)) * gam_@{se} 
@#endfor 
;

[name = 'consumer price index - flex', type = 'flex']
0 = 
@#for se in 1:nsector
     + (omega_c_@{se}*ss_gam_@{se}^(1-eta)) * gam_flex_@{se} 
@#endfor 
;

/// Investment price index (97)
// CHECK OK
[name = 'investment price index', type = 'sticky']
(ss_gam_I^(1-eta))*gam_I = 
@#for se in 1:nsector
     + omega_i_@{se}*ss_gam_@{se}^(1-eta) * gam_@{se} 
@#endfor 
;

[name = 'investment price index', type = 'flex']
(ss_gam_I^(1-eta))*gam_I_flex = (
@#for se in 1:nsector
     + omega_i_@{se}*ss_gam_@{se}^(1-eta) * gam_flex_@{se} 
@#endfor 
);

/// Government expenditure price index (98)
// CHECK OK
[name = 'government price index', type = 'sticky']
ss_gam_G^(1-eta_g) * gam_G = 
@#for se in 1:nsector 
    + omega_g_@{se}*ss_gam_@{se}^(1-eta_g) * gam_@{se} 
@#endfor
;

[name = 'government price index', type = 'flex']
ss_gam_G^(1-eta_g) * gam_G_flex = 
@#for se in 1:nsector 
    + omega_g_@{se}*ss_gam_@{se}^(1-eta_g) * gam_flex_@{se} 
@#endfor
;

/// Aggregate labour supply index (102)
// CHECK OK
[name = 'labour supply index', type = 'sticky']
ss_lab^((xi+1)/xi) * lab - (
@#for se in 1:nsector 
     + omega_ns_@{se}^(-1/xi) * ss_n_@{se}^((xi+1)/xi) * n_@{se} 
@#endfor
) = 0;

[name = 'labour supply index - flex', type = 'flex']
ss_lab^((xi+1)/xi) * lab_flex - (
@#for se in 1:nsector 
     + omega_ns_@{se}^(-1/xi) * ss_n_@{se}^((xi+1)/xi) * n_flex_@{se} 
@#endfor
) = 0;

/// Investment price inflation (103)
// CHECK OK
[name = 'investment price inflation', type = 'sticky']
pie_I - pie = gam_I - gam_I(-1);

/// Wage inflation (104)
// CHECK OK 
#WAGE_BILL = (
@#for se in 1 : nsector
+ ss_w_@{se} * ss_n_@{se}
@#endfor
);

[name = 'wage inflation', type = 'sticky']
pie_W - (
@#for se in 1 : nsector
+ ss_w_@{se} * ss_n_@{se} / WAGE_BILL * pie_W_@{se}
@#endfor
) = 0 ;

/// Aggregate wage index
// CHECK OK
[name = 'wage index', type = 'sticky']
w = 
(
@#for se in 1:nsector                                       
	+ omega_ns_@{se} * (ss_w_@{se} / ss_w)^(1 + xi) * w_@{se}    
@#endfor
) ;

[name = 'wage index - flex', type = 'flex']
w_flex = 
(
@#for se in 1:nsector                                       
	+ omega_ns_@{se} * (ss_w_@{se} / ss_w)^(1 + xi) * w_flex_@{se}    
@#endfor
) ;

/// Investment market clearing
// CHECK OK
[name = 'investment market clearing', type = 'sticky']
ss_invest * invest = 
(
@#for se in 1:nsector
        +  ss_z_@{se} * z_@{se}
@#endfor 
) ;

[name = 'investment market clearing - flex', type = 'flex']
ss_invest * invest_flex = 
(
@#for se in 1:nsector
        +  ss_z_@{se} * z_flex_@{se}
@#endfor 
) ;

/// Consumption choice for Ricardian consumers
// NOTE: DR - HAVE EXCLUDED SHOCKS 
[name = 'consumption ricardian', type = 'sticky']
(1 + h^2 * bet) * c_r = bet * h * c_r(+1) - 
                        (1 - h) * (1 - bet * h)*lam_r + 
                        h * c_r(-1) + (1 - h) * cshift - (1 - h) * bet * cshift(+1) ;

[name = 'consumption ricardian - flex', type = 'flex']
(1 + h^2 * bet) * c_r_flex = bet * h * c_r_flex(+1) - 
                             (1 - h) * (1 - bet * h)*lam_r_flex + 
                             h * c_r_flex(-1)  + (1 - h) * cshift - (1 - h) * bet * cshift(+1) ;


/// Euler equation for Ricardian consumers
// CHECK OK
[name = 'euler ricardian', type = 'sticky']
r + lam_r(+1) - pie(+1) = lam_r;

[name = 'euler ricardian - flex', type = 'flex']
 rr_flex + lam_r_flex(+1) = lam_r_flex  ;

/// Capital stock choice
// CHECK OK
@#for se in 1:nsector
[name = 'capital choice', type = 'sticky']
 -lam_r + lam_r(+1) + bet * (1-delt)*q_@{se}(+1) + 
            ((bet * ss_rk_@{se}) / rk_wedge) * rk_@{se}(+1) = q_@{se};   
@#endfor

@#for se in 1:nsector
[name = 'capital choice - flex', type = 'flex']
 -lam_r_flex + lam_r_flex(+1) + bet * (1-delt)*q_flex_@{se}(+1) + 
                 ((bet * ss_rk_@{se}) / rk_wedge) * rk_flex_@{se}(+1) = q_flex_@{se}  ;   
@#endfor

/// Investment choice
// CHECK NOT OK: CHECK COEFFICIENTS
@#for se in 1:nsector
[name = 'investment choice', type = 'sticky']
	(1 + bet) * z_@{se} = q_@{se}/Spp + bet * z_@{se}(+1) + z_@{se}(-1) ;  
@#endfor

@#for se in 1:nsector
[name = 'investment choice - flex', type = 'flex']
	(1 + bet) * z_flex_@{se} = q_flex_@{se}/Spp + bet * z_flex_@{se}(+1) + z_flex_@{se}(-1) ;  
@#endfor

/// Consumption choice for non-Ricardian consumers
// CHECK OK
[name = 'consumption non-ricardian', type = 'sticky']
c_nr = (1/ss_c_nr) * (ss_w*ss_lab*(w + lab)) ;
//+ (ss_trans*trans_hat));

[name = 'consumption non-ricardian - flex', type = 'flex']
c_nr_flex = (1/ss_c_nr) * (ss_w*ss_lab*(w_flex + lab_flex)) ;

/// Marginal utility of consumption for non-Ricardian consumers
// CHECK OK
[name = 'muc non-ricardian', type = 'sticky']
(1 - h) * (1 - bet * h) * lam_nr = - (1 + bet * h^2)*c_nr + 
                                    bet * h * c_nr(+1) +
                                    h * c_nr(-1) + 
                                    (1 - h) * cshift - (1 - h) * bet * cshift(+1);

[name = 'muc non-ricardian - flex', type = 'flex']
(1 - h) * (1 - bet * h) * lam_nr_flex = - (1 + bet * h^2)*c_nr_flex + 
                                        bet * h * c_nr_flex(+1) +
                                        h * c_nr_flex(-1) + 
                                       (1 - h) * cshift - (1 - h) * bet * cshift(+1);

/// Aggregate marginal utility (116)
// CHECK OK
[name = 'muc aggregate', type = 'sticky']
ss_lam * lam = omega_r * ss_lam_r * lam_r + 
       (1 - omega_r) * ss_lam_nr * lam_nr;

[name = 'muc aggregate - flex', type = 'flex']
ss_lam * lam_flex = omega_r * ss_lam_r  * lam_r_flex + 
           (1 - omega_r) * ss_lam_nr * lam_nr_flex;

/// Aggregate consumption (115)
// CHECK OK
[name = 'consumption aggregate', type = 'sticky']
ss_c * c = omega_r * ss_c_r * c_r + 
            (1 - omega_r) * ss_c_nr * c_nr ;

[name = 'consumption aggregate - flex', type = 'flex']
ss_c * c_flex = omega_r * ss_c_r * c_r_flex + 
                (1 - omega_r) * ss_c_nr * c_nr_flex ;

/// Consumption variety choice (99)
// CHECK OK
@#for se in 1:nsector
    [name = 'consumption variety', type = 'sticky']
	c_@{se} = c - eta * gam_@{se};
@#endfor

@#for se in 1:nsector
    [name = 'consumption variety - flex', type = 'flex']
	c_flex_@{se} = c_flex - eta * gam_flex_@{se};
@#endfor

/// Investment variety choice (100)
// CHECK OK
@#for se in 1:nsector
    [name = 'investment variety', type = 'sticky']
	i_@{se} = invest - eta * (gam_@{se} - gam_I);
@#endfor

@#for se in 1:nsector
    [name = 'investment variety - flex', type = 'flex']
	i_flex_@{se} = invest_flex - eta * (gam_flex_@{se} - gam_I_flex);
@#endfor

/// Government expenditure variety choice (101)
// CHECK OK
@#for se in 1:nsector
    [name = 'government variety', type = 'sticky']
	g_@{se} = g - eta * (gam_@{se} - gam_G);
@#endfor

@#for se in 1:nsector
    [name = 'government variety - flex', type = 'flex']
	g_flex_@{se} = g_flex - eta * (gam_flex_@{se} - gam_G_flex);
@#endfor

/// Wage growth (117)
// CHECK OK
#TEMP_W = 1 + bet * chi_w ;
#KAPPA  = (1 - bet * thet_w) * (1 - thet_w) / thet_w * xi / (xi + eps_w) ;
@#for se in 1:nsector
        [name = 'wage growth', type = 'sticky']
        pie_W_@{se} = bet / TEMP_W * pie_W_@{se}(+1) - 
                        KAPPA * lam + 
                        KAPPA * (nu - 1 / xi) * lab + 
                        KAPPA / xi * n_@{se} - 
                        KAPPA * w_@{se} + 
                        chi_w / TEMP_W * pie_W_@{se}(-1) ;
@#endfor

/// Wage choice - flexible price
// CHECK OK
@#for se in 1:nsector
    [name = 'wage choice - flex', type = 'flex']
    lam_flex = (nu - 1 / xi) * lab_flex + 1 / xi * n_flex_@{se} - 
                w_flex_@{se} ;
@#endfor

/// Gross output in sector se (118)
// CHECK OK
#TEMP_GO = (vphi - 1) / vphi ;
@#for se in 1:nsector
    [name = 'gross output', type = 'sticky']
    (ss_y_@{se} / ss_a_@{se})^TEMP_GO * (y_@{se} - (a_@{se} + tfp )) = omega_y_@{se}^(1/vphi) * ss_f_@{se}^TEMP_GO * f_@{se} + 
                                                                (1 - omega_y_@{se})^(1/vphi) * ss_x_@{se}^TEMP_GO * x_@{se} ;
@#endfor

@#for se in 1:nsector
    [name = 'gross output - flex', type = 'flex']
    (ss_y_@{se} / ss_a_@{se})^TEMP_GO * (y_flex_@{se} - (a_@{se} + tfp )) = omega_y_@{se}^(1/vphi) * ss_f_@{se}^TEMP_GO * f_flex_@{se} + 
                                                                   (1 - omega_y_@{se})^(1/vphi) * ss_x_@{se}^TEMP_GO * x_flex_@{se} ;
@#endfor

/// Marginal cost in sector se (119)
// CHECK OK
@#for se in 1:nsector
    [name = 'marginal cost', type = 'sticky']
	mc_@{se} + (a_@{se} + tfp ) = omega_y_@{se} * (ss_gam_f_@{se} * (1 + tau_c_@{se}) / (ss_gam_@{se} * ss_a_@{se}))^(1 - vphi) * gam_f_@{se} + 
                          (1 - omega_y_@{se}) * (ss_gam_x_@{se} * (1 + tau_c_@{se}) / (ss_gam_@{se} * ss_a_@{se}))^(1 - vphi) * gam_x_@{se} - 
                          gam_@{se} ;
@#endfor

@#for se in 1:nsector
    [name = 'marginal cost - flex', type = 'flex']
	(a_@{se} + tfp ) = omega_y_@{se} * (ss_gam_f_@{se} * (1 + tau_c_@{se}) / (ss_gam_@{se} * ss_a_@{se}))^(1 - vphi) * gam_f_flex_@{se} + 
               (1 - omega_y_@{se}) * (ss_gam_x_@{se} * (1 + tau_c_@{se}) / (ss_gam_@{se} * ss_a_@{se}))^(1 - vphi) * gam_x_flex_@{se} - 
               gam_flex_@{se} ;
@#endfor

/// Factor demand in sector se (120)
// CHECK OK
@#for se in 1:nsector
    [name = 'factor demand go', type = 'sticky']
	vphi * gam_f_@{se} + f_@{se} - vphi * gam_x_@{se} - x_@{se} = 0;
@#endfor

@#for se in 1:nsector
    [name = 'factor demand go - flex', type = 'flex']
	vphi * gam_f_flex_@{se} + f_flex_@{se} - vphi * gam_x_flex_@{se} - x_flex_@{se} = 0;
@#endfor

/// Gam f price index in sector se (121)
// CHECK OK
@#for se in 1:nsector
    [name = 'gam f index', type = 'sticky']
	(ss_gam_f_@{se}^(1-zeta)) * gam_f_@{se} = 
        (omega_nd_@{se} * (ss_w_@{se}^(1 - zeta)) * w_@{se}) + 
        (1-omega_nd_@{se}) * (ss_rk_@{se}^(1 - zeta))*rk_@{se} ;
@#endfor

@#for se in 1:nsector
    [name = 'gam f index - flex', type = 'flex']
	(ss_gam_f_@{se}^(1-zeta)) * gam_f_flex_@{se} = 
        (omega_nd_@{se} * (ss_w_@{se}^(1 - zeta)) * w_flex_@{se}) + 
        (1-omega_nd_@{se}) * (ss_rk_@{se}^(1 - zeta))*rk_flex_@{se} ;
@#endfor

/// Intermediate good price index in sector se (122)
// CHECK OK
@#for se in 1:nsector
    [name = 'gam x index', type = 'sticky']
    (ss_gam_x_@{se}^(1 - psi))*gam_x_@{se} = 
    (
    @#for se_2 in 1:nsector 
        + ii_weights_@{se_2}_@{se} * (ss_gam_@{se_2}^(1-psi)) * gam_@{se_2}
	@#endfor
    ) ;
@#endfor

@#for se in 1:nsector
    [name = 'gam x index - flex', type = 'flex']
    (ss_gam_x_@{se}^(1 - psi))*gam_x_flex_@{se} = 
    (
    @#for se_2 in 1:nsector 
        + ii_weights_@{se_2}_@{se} * (ss_gam_@{se_2}^(1-psi)) * gam_flex_@{se_2}
	@#endfor
    ) ;
@#endfor

/// Labour demand in sector se (123)
// CHECK OK
@#for se in 1:nsector
    [name = 'labour demand', type = 'sticky']
	f_@{se} - zeta * w_@{se} + zeta * gam_f_@{se} - n_@{se} = 0;
@#endfor

@#for se in 1:nsector
    [name = 'labour demand - flex', type = 'flex']
	f_flex_@{se} - zeta * w_flex_@{se} + zeta * gam_f_flex_@{se} - n_flex_@{se} = 0;
@#endfor

/// Capital demand in sector se (124)
// CHECK OK
@#for se in 1:nsector
    [name = 'capital demand', type = 'sticky']
	f_@{se} - zeta * rk_@{se} + zeta * gam_f_@{se} - ks_@{se} = 0; 
@#endfor

@#for se in 1:nsector
    [name = 'capital demand - flex', type = 'flex']
	f_flex_@{se} - zeta * rk_flex_@{se} + zeta * gam_f_flex_@{se} - ks_flex_@{se} = 0; 
@#endfor

/// Capital services
// CHECK OK
@#for se in 1:nsector
    [name = 'capital services', type = 'sticky']
	ks_@{se} = capu_@{se} + k_@{se}(-1);
@#endfor

@#for se in 1:nsector
    [name = 'capital services - flex', type = 'flex']
	ks_flex_@{se} = capu_flex_@{se} + k_flex_@{se}(-1);
@#endfor

/// Capital utilisation in sector 
// CHECK OK
@#for se in 1:nsector
    [name = 'capital utilisation', type = 'sticky']
	capu_@{se} = capup * rk_@{se};   
@#endfor

@#for se in 1:nsector
    [name = 'capital utilisation - flex', type = 'flex']
	capu_flex_@{se} = capup * rk_flex_@{se};   
@#endfor

/// (125) Intermediate good k demand in sector se
// CHECK OK
@#for se in 1:nsector
	@#for se2 in 1:nsector
    [name = 'int goods demand', type = 'sticky']
    x_@{se} - xx_@{se2}_@{se} - psi * gam_@{se2} + psi * gam_x_@{se} ;
	@#endfor
@#endfor

@#for se in 1:nsector
	@#for se2 in 1:nsector
    [name = 'int goods demand - flex', type = 'flex']
    x_flex_@{se} - xx_flex_@{se2}_@{se} - psi * gam_flex_@{se2} + psi * gam_x_flex_@{se} ;
	@#endfor
@#endfor

/// (126) Definition of relative price in sector se
// CHECK OK
@#for se in 1:nsector
    [name = 'industry relative price', type = 'sticky']
	gam_@{se} + pie - gam_@{se}(-1) = pie_@{se} ;
@#endfor


/// (127) Phillips curve in sector se
// CHECK OK
@#for se in 1:nsector
    #TEMP_PC_@{se} = 1 + bet * chi_p_@{se} ;
    #KAPPA_PC_@{se} = (1 - thet_@{se}) * (1 - thet_@{se} * bet) / thet_@{se} ;
    [name = 'phillips curve', type = 'sticky']
	pie_@{se} = bet / TEMP_PC_@{se} * pie_@{se}(+1) + KAPPA_PC_@{se} / TEMP_PC_@{se} * mc_@{se} + chi_p_@{se} / TEMP_PC_@{se} * pie_@{se}(-1) + mu_@{se};
@#endfor

/// (128) Link between wage inflation and real wages in sector se
//CHECK OK
@#for se in 1:nsector
    [name = 'wage inflation link', type = 'sticky']
	pie_W_@{se} = w_@{se} + pie - w_@{se}(-1);
@#endfor

/// (129) Market clearing in sector se
//CHECK OK
@#for se in 1:nsector
    [name = 'market clearing', type = 'sticky']
	ss_y_@{se} * y_@{se} = ss_c_@{se} * c_@{se} + 
                           ss_i_@{se} * i_@{se} + 
                           ss_g_@{se} * g_@{se} + (
	@#for k in 1:nsector
		+ ss_xx_@{se}_@{k} * xx_@{se}_@{k} 
	@#endfor 
    );
@#endfor

@#for se in 1:nsector
    [name = 'market clearing - flex', type = 'flex']
	ss_y_@{se} * y_flex_@{se} = ss_c_@{se} * c_flex_@{se} + 
                                ss_i_@{se} * i_flex_@{se} + 
                                ss_g_@{se} * g_flex_@{se} + (
	@#for k in 1:nsector
		+ ss_xx_@{se}_@{k} * xx_flex_@{se}_@{k} 
	@#endfor 
    );
@#endfor

/// (130) Taylor rule
[name = 'Inflation lag 0', type = 'flex']
lpie_0 = pie;

[name = 'Inflation lag 1', type = 'flex']
lpie_1 = pie(-1);

[name = 'Inflation lag 2', type = 'flex']
lpie_2 = lpie_1(-1);

[name = 'Inflation lag 3', type = 'flex']
lpie_3 = lpie_2(-1);

[name = 'Inflation lag 4', type = 'flex']
lpie_4 = lpie_3(-1);

[name = 'Inflation lag 5', type = 'flex']
lpie_5 = lpie_4(-1);

[name = 'Inflation lag 6', type = 'flex']
lpie_6 = lpie_5(-1);

[name = 'Inflation lag 7', type = 'flex']
lpie_7 = lpie_6(-1);

[name = 'Inflation lag 8', type = 'flex']
lpie_8 = lpie_7(-1);

[name = 'Inflation lag 9', type = 'flex']
lpie_9 = lpie_8(-1);

[name = 'Inflation lag 10', type = 'flex']
lpie_10 = lpie_9(-1);

[name = 'Inflation lag 11', type = 'flex']
lpie_11 = lpie_10(-1);

[name = 'Inflation lag 12', type = 'flex']
lpie_12 = lpie_11(-1);

[name = 'Inflation lag 13', type = 'flex']
lpie_13 = lpie_12(-1);

[name = 'Inflation lag 14', type = 'flex']
lpie_14 = lpie_13(-1);

[name = 'Inflation lag 15', type = 'flex']
lpie_15 = lpie_14(-1);

[name = 'Inflation lag 16', type = 'flex']
lpie_16 = lpie_15(-1);

//[name = 'arithmetic mean of inflation', type = 'flex']
pie_average = (
@#for se in 0:n_ait_loop-1
        +lpie_@{se} 
@#endfor 
)/n_ait ;


[name = 'Inflation lead 1', type = 'flex']
fpie_1 = pie(+1);

[name = 'Inflation lead 2', type = 'flex']
fpie_2 = fpie_1(+1);

[name = 'Inflation lead 3', type = 'flex']
fpie_3 = fpie_2(+1);

[name = 'Inflation lead 4', type = 'flex']
fpie_4 = fpie_3(+1);

//[name = 'forecast of inflation', type = 'flex']
pie_forecast = (
@#for se in 1:4
        +fpie_@{se} 
@#endfor 
)/4 ;


//CHECK OK
//[name = 'taylor rule', type = 'sticky']
//r = rho_r * r(-1) + (1 - rho_r) * (phi_pie * pie_average + phi_gap * gap) + sig_r * eps_r;

//r = rho_r * r(-1) + (1 - rho_r) * (phi_pie * pie + phi_pie_core * pie_core + (phi_gap/4) * gap) + sig_r * eps_r;

//[name = 'AIT rule (equations (4) and (5) in BIS Working Paper No 1156)', type = 'sticky']
//r = rho_r * r(-1) + (1 - rho_r) * (phi_pie * pie_average + phi_gap * gap) + sig_r * eps_r;
[name = 'discounted sum of past deviations of inflation from target', type = 'flex']
pie_average_NYFed = pie - 0 + 0.93*pie_average_NYFed(-1);

//generic rule
[name = 'taylor rule', type = 'sticky']
r = rho_r * r(-1) + (1 - rho_r) * (phi_pie * pie + phi_pie_bar * pie_average +  phi_pie_core * pie_core + phi_pie_NYFed * pie_average_NYFed + phi_pie_headmin * pie_headmin + phi_pie_forecast * pie_forecast + phi_gap * gap) + sig_r * eps_r; 

/// (131) output gap
//CHECK OK
[name = 'output gap', type = 'common']
gap = y_va - y_va_flex;

/// (133) Aggregate value added
//CHECK OK
[name = 'aggregate value added', type = 'sticky']
ss_y_va * y_va = (
@#for se in 1:nsector
	+ ss_nva_@{se} * y_va_@{se}
@#endfor
);

[name = 'aggregate value added - flex', type = 'sticky - flex']
ss_y_va * y_va_flex = (
@#for se in 1:nsector
	+ ss_nva_@{se} * y_va_flex_@{se}
@#endfor
);

/// Sectoral value added
// CHECK OK
@#for se in 1 : nsector
[name = 'sectoral value added', type = 'sticky']
 y_va_@{se} = ss_gam_@{se} * ss_y_@{se} / ss_nva_@{se} * y_@{se} - 
                                ss_gam_x_@{se} * ss_x_@{se} / ss_nva_@{se} * x_@{se} ;

@#endfor

@#for se in 1 : nsector
[name = 'sectoral value added - flex', type = 'flex']
 y_va_flex_@{se} = ss_gam_@{se} * ss_y_@{se} / ss_nva_@{se} * y_flex_@{se} - 
                                    ss_gam_x_@{se} * ss_x_@{se} / ss_nva_@{se} * x_flex_@{se} ;

@#endfor

/// Sectoral value added inflation
//CHECK OK
@#for se in 1 : nsector
    [name = 'sectoral value added inflation', type = 'sticky']
    pie_va_@{se} = pie + gam_va_@{se} - gam_va_@{se}(-1) ;
@#endfor

/// Sectoral "nominal value added"
// CHECK OK
@#for se in 1 : nsector
[name = 'sectoral nominal value added', type = 'sticky']
    ss_nva_@{se} * nva_@{se} = ss_gam_@{se} * ss_y_@{se} * (gam_@{se} + y_@{se}) - 
                                ss_gam_x_@{se} * ss_x_@{se} * (gam_x_@{se} + x_@{se}) ;
@#endfor

@#for se in 1 : nsector
[name = 'sectoral nominal value added - flex', type = 'flex']
    ss_nva_@{se} * nva_flex_@{se} = ss_gam_@{se} * ss_y_@{se} * (gam_flex_@{se} + y_flex_@{se}) - 
                                    ss_gam_x_@{se} * ss_x_@{se} * (gam_x_flex_@{se} + x_flex_@{se}) ;
@#endfor


/// Sectoral value added price
// CHECK OK
@#for se in 1 : nsector
[name = 'sectoral value added price', type = 'sticky']
    y_va_@{se} = nva_@{se} - gam_va_@{se} ;
@#endfor

@#for se in 1 : nsector
[name = 'sectoral value added price - flex', type = 'flex']
    y_va_flex_@{se} = nva_flex_@{se} - gam_va_flex_@{se} ;
@#endfor

/// (132) Year-on-year inflation
// CHECK OK
[name = 'inflation - yoy', type = 'sticky']
pie_ye = pie + Lpie + L2pie + L2pie(-1);

// Lagged inflation
Lpie = pie(-1) ;

// 2 period lagged inflation
L2pie = Lpie(-1) ;

/// Hours worked
//CHECK OK
[name = 'hours worked', type = 'sticky']
ss_hrs * hrs = +
@#for se in 1 : nsector
    + ss_n_@{se} * n_@{se}
@#endfor
;

[name = 'hours worked - flex', type = 'flex']
ss_hrs * hrs_flex = +
@#for se in 1 : nsector
    + ss_n_@{se} * n_flex_@{se}
@#endfor
;

/// Real interest rate - sticky price
//CHECK OK
[name = 'real interest rate', type = 'sticky']
rr = r - pie(+1) ;

/// (134) Productivity shock process
//CHECK OK
@#for se in 1:nsector
[name = 'industry productivity', type = 'common']
	a_@{se} = rho_a_@{se} * a_@{se}(-1) + sig_a_@{se} * eps_a_@{se} + c_tfp_mm_@{se} * sig_mm * eps_m;
@#endfor

@#for se in 1 : nsector
[name = 'industry consumption shifter', type = 'common']
    om_c_@{se} = 0.9* om_c_@{se}(-1);

@#endfor

/// (135) Aggregate government expenditure shock process
//CHECK OK
[name = 'government expenditure', type = 'common']
g = rho_g * g(-1) + sig_g * eps_g;

g_flex = rho_g * g_flex(-1) + sig_g * eps_g;

// Consumption shifter
cshift = 0.9*cshift(-1) + sig_cshift * eps_cshift ;

/// (136) Transfers
// DR - DELETED THIS FOR NOW
//trans = rho_trans * trans(-1) + eps_trans;
trans = 0.9*trans(-1) + eps_trans ;

/// (137) Markup shock process
//CHECK OK
@#for se in 1:nsector
[name = 'industry markup', type = 'common']
	mu_@{se} = rho_mu_@{se} * mu_@{se}(-1) + sig_mu_@{se} * eps_mu_@{se} + c_mu_mm_@{se} * sig_mm * eps_m;//
@#endfor


/// (138) Aggreagte TFP shock process
tfp = rho_tfp * tfp(-1) + sig_tfp * eps_tfp;

// ========================================================================
// Additional equations

// Core inflation
//#TOTAL_CORE_WEIGHT = ss_c - ss_gam_1 * ss_c_1 - ss_gam_2 * ss_c_2 ;
//[name = 'core inflation', type = 'common']
//pie_core = 
//@#for se in 3 : nsector
//    + ss_gam_@{se} * ss_c_@{se} * pie_@{se} / TOTAL_CORE_WEIGHT
//@#endfor
//;

// Core inflation proxied by service sector inflation
#TOTAL_CORE_WEIGHT = ss_c - ss_c_1 - ss_c_2 - ss_c_3 - ss_c_4 - ss_c_5 ;
[name = 'core inflation', type = 'common']
pie_core = 
@#for se in 6 : nsector
    + ss_c_@{se} * pie_@{se} / TOTAL_CORE_WEIGHT
@#endfor
;


#TOTAL_NONSERVICE_WEIGHT = ss_nva_1 + ss_nva_2 + ss_nva_3 + ss_nva_4 + ss_nva_5;
#TOTAL_SERVICE_WEIGHT = 1 - TOTAL_NONSERVICE_WEIGHT ;

y_va_nonservice = 
@#for se in 1 : 5
    + ss_nva_@{se} * y_va_@{se} / TOTAL_NONSERVICE_WEIGHT
@#endfor
;
y_va_service = 
@#for se in 6 : nsector
    + ss_nva_@{se} * y_va_@{se} / TOTAL_SERVICE_WEIGHT
@#endfor
;

// Some aggreagte variables
Lpie_va_2 = pie_va_2(-1) ;  //pie_va_2(-1)
L2pie_va_2 = Lpie_va_2(-1) ; // pie_va_2(-2)
Ly_va_2 = y_va_2(-1) ; // y_va_2(-1)
L2y_va_2 = Ly_va_2(-1) ; // y_va_2(-2)
L3y_va_2 = L2y_va_2(-1) ;// y_va_2(-3)
Lpie_va_5 = pie_va_5(-1) ;  //pie_va_2(-1)
L2pie_va_5 = Lpie_va_5(-1) ; // pie_va_2(-2)
Ly_va_5 = y_va_5(-1) ; // y_va_2(-1)
L2y_va_5 = Ly_va_5(-1) ; // y_va_2(-2)
L3y_va_5 = L2y_va_5(-1) ;// y_va_2(-3)
Ly_va = y_va(-1) ; // y_va_2(-1)
L2y_va = Ly_va(-1) ; // y_va_2(-2)
L3y_va = L2y_va(-1) ;// y_va_2(-3)

INFL         = 400*pie;
FFR          = 400*r;
FFRDELTA     = 400*(r-r(-1));
YGROWTH      = 400*(y_va-y_va(-1));
YGAP         = 100*gap;
VA_MINING    = 100*(y_va_2-L3y_va_2(-1)); //100*(y_va_2+y_va_2(-1)+y_va_2(-2)+y_va_2(-3)-(y_va_2(-4)+y_va_2(-5)+y_va_2(-6)+y_va_2(-7)));
DEFL_MINING  = 100*(pie_va_2 + Lpie_va_2 + L2pie_va_2 + L2pie_va_2(-1));
VA_MANUFACTURING = 100*(y_va_5-L3y_va_5(-1)); //100*(y_va_5+y_va_5(-1)+y_va_5(-2)+y_va_5(-3)-(y_va_5(-4)+y_va_5(-5)+y_va_5(-6)+y_va_5(-7)));
DEFL_MANUFACTURING = 100*(pie_va_5 + Lpie_va_5 + L2pie_va_5 + L2pie_va_5(-1));
INFL_Annual  = 100*(pie + lpie_1 + lpie_2 + lpie_3);
YGROWTH_Annual = 100*(y_va-L3y_va(-1));
Pie_wr = w-w(-1);
eps_m = eps_mm; 
pie_headmin = 0.5*pie + 0.5*pie_2;
end;