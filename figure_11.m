clc
clear
close all

% Parameters
eta = 0.9;
a_1 = 0.3;
a_2 = 0.7;
epsilon = 0.01;

L_1 = 60;           % Link length from S to IRS
L_2 = 60;           % Link length from IRS to M_1
L_M1 = L_1 + L_2;   % NLOS through IRS total path
L_M2 = 40;          % LOS path from M_1 to M_2
L_M11 = 60;         % LOS path length
L_B_M2 = 100;        % LOS from Buoy to M_2

c_lambda = 0.056;   % extinction coefficient for pure sea water
W_0 = 0.01;
lambda = 470e-9;
D_R = 0.05;         % receiver aperture diameter
R_s = 1;
zeta = 1;
N = 64;

gamma_th_1 = 2^R_s - 1;
gamma_th_2 = gamma_th_1;

% Parameters for B.L = 16.5L/min in Salt water
w_s = 0.4951;
lamb_s = 0.1368;
a_s = 0.0161;
b_s = 3.2033;
c_s = 82.1030;

% Parameters for B.L = 4.7L/min
w_f = 0.2109;
lamb_f = 0.4603;
a_f = 1.2526;
b_f = 1.1501;
c_f = 41.3258;

% Path loss
h_l_M1  = exp(-c_lambda * L_M1);    % NLOS via IRS
h_l_M11 = exp(-c_lambda * L_M11);   % LOS
h_l_M2 = exp(-c_lambda * L_M2 );    % LOS from M_1 to M_2
h_l_B_M2 = exp(-c_lambda * L_B_M2 );    % LOS from M_1 to M_2

% Aperture / beam spread for NLOS via IRS
W_Z_M1 = W_0 * sqrt(1 + ((lambda * L_M1) / (pi * W_0^2))^2);
v_1_M1 = sqrt(pi/2) * (D_R / (2 * W_Z_M1));
A_0_M1 = (erf(v_1_M1))^2;

% Aperture / beam spread for LOS
W_Z_M11 = W_0 * sqrt(1 + ((lambda * L_M11) / (pi * W_0^2))^2);
v_1_M11 = sqrt(pi/2) * (D_R / (2 * W_Z_M11));
A_0_M11 = (erf(v_1_M11))^2;

% Aperture / beam spread for LOS from M_1 to M_2
W_Z_2 = W_0 * sqrt(1 + ((lambda * L_M2) / (pi * W_0^2))^2);
v_2 = sqrt(pi/2) * (D_R / (2 * W_Z_2));
A_0_2 = (erf(v_2))^2;

% Aperture / beam spread for LOS from buoy to M_2
W_Z_B_2 = W_0 * sqrt(1 + ((lambda * L_B_M2) / (pi * W_0^2))^2);
v_B_2 = sqrt(pi/2) * (D_R / (2 * W_Z_B_2));
A_0_B_2 = (erf(v_B_2))^2;

% SNR range
gamma_dB = 0:10:100;
gamma_values = 10.^(gamma_dB/10);

% Output arrays
P_out_ana_LOS_s  = zeros(size(gamma_dB));
P_out_ana_LOS_f  = zeros(size(gamma_dB));

P_out_ana_NLOS_s = zeros(size(gamma_dB));
P_out_ana_NLOS_f = zeros(size(gamma_dB));

P_outage_M2_s = zeros(size(gamma_dB));
P_outage_M2_f = zeros(size(gamma_dB));

P_outage_B_M2_s = zeros(size(gamma_dB));
P_outage_B_M2_f = zeros(size(gamma_dB));

P_out_ana_auto_s = zeros(size(gamma_dB));
P_out_ana_auto_f = zeros(size(gamma_dB));

P_out_B_M2_s = zeros(size(gamma_dB));
P_out_B_M2_f = zeros(size(gamma_dB));

P_out_mc_LOS_s   = zeros(size(gamma_dB));
P_out_mc_LOS_f   = zeros(size(gamma_dB));

P_out_mc_NLOS_s  = zeros(size(gamma_dB));
P_out_mc_NLOS_f  = zeros(size(gamma_dB));

P_out_mc_M2_s = zeros(size(gamma_dB));
P_out_mc_M2_f = zeros(size(gamma_dB));

P_out_mc_auto_s  = zeros(size(gamma_dB));
P_out_mc_auto_f  = zeros(size(gamma_dB));

P_out_mc_B_M2_s = zeros(size(gamma_dB));
P_out_mc_B_M2_f = zeros(size(gamma_dB));

P_outage_mc_B_M2_s = zeros(size(gamma_dB));
P_outage_mc_B_M2_f = zeros(size(gamma_dB));

N_sim = 1e5;

for i = 1:length(gamma_values)
    rho = gamma_values(i);

    % ================= Analytical LOS in salt water =================
    alpha_LOS = gamma_th_2 / ((a_2 - (a_1 *  gamma_th_2)) * rho * eta);
    beta_LOS  = gamma_th_1 / ((a_1 - (a_2 * epsilon * gamma_th_1)) * rho * eta);
    V_LOS = min(alpha_LOS, beta_LOS);

    arg_1_LOS_s = V_LOS / (A_0_M11 * h_l_M11 * lamb_s);
    product_term_1_LOS_s = (w_s * zeta^2) * arg_1_LOS_s^(zeta^2);
    G_Meijer_1_LOS_s = meijerG(1-zeta^2,1,[1-zeta^2,0],-zeta^2,arg_1_LOS_s);
    term_1_LOS_s = product_term_1_LOS_s * G_Meijer_1_LOS_s;

    arg_2_LOS_s = V_LOS / (A_0_M11 * h_l_M11 * b_s);
    product_term_2_LOS_s = ((1-w_s) * zeta^2) / (gamma(a_s) * c_s) * arg_2_LOS_s^(zeta^2);

    a_par_s = 1 - zeta^2/c_s;
    b_par_1_s = a_s - (zeta^2/c_s);
    b_par_2_s = -zeta^2/c_s;

    G_Meijer_2_LOS_s = meijerG(a_par_s,1,[b_par_1_s,0],b_par_2_s,arg_2_LOS_s^c_s);
    term_2_LOS_s = product_term_2_LOS_s * G_Meijer_2_LOS_s;

    P_out_ana_LOS_s(i) = term_1_LOS_s + term_2_LOS_s;

    % ================= Analytical LOS in salt water =================
    arg_1_LOS_f = V_LOS / (A_0_M11 * h_l_M11 * lamb_f);
    product_term_1_LOS_f = (w_f * zeta^2) * arg_1_LOS_f^(zeta^2);
    G_Meijer_1_LOS_f = meijerG(1-zeta^2,1,[1-zeta^2,0],-zeta^2,arg_1_LOS_f);
    term_1_LOS_f = product_term_1_LOS_f * G_Meijer_1_LOS_f;

    arg_2_LOS_f = V_LOS / (A_0_M11 * h_l_M11 * b_f);
    product_term_2_LOS_f = ((1-w_f) * zeta^2) / (gamma(a_f) * c_f) * arg_2_LOS_f^(zeta^2);

    a_par_f = 1 - zeta^2/c_f;
    b_par_1_f = a_f - (zeta^2/c_f);
    b_par_2_f = -zeta^2/c_f;

    G_Meijer_2_LOS_f = meijerG(a_par_f,1,[b_par_1_f,0],b_par_2_f,arg_2_LOS_f^c_f);
    term_2_LOS_f = product_term_2_LOS_f * G_Meijer_2_LOS_f;

    P_out_ana_LOS_f(i) = term_1_LOS_f + term_2_LOS_f;

    % ================= Analytical NLOS via IRS =================
    alpha_NLOS = gamma_th_2 / ((a_2 - (a_1 *  gamma_th_2)) * rho * N^2 * eta);
    beta_NLOS  = gamma_th_1 / ((a_1 - (a_2 * epsilon * gamma_th_1)) * rho * N^2 * eta);
    V_NLOS = min(alpha_NLOS, beta_NLOS);

    arg_1_NLOS_s = V_NLOS / (A_0_M1 * h_l_M1 * lamb_s);
    product_term_1_NLOS_s = (w_s * zeta^2) * arg_1_NLOS_s^(zeta^2);
    G_Meijer_1_NLOS_s = meijerG(1-zeta^2,1,[1-zeta^2,0],-zeta^2,arg_1_NLOS_s);
    term_1_NLOS_s = product_term_1_NLOS_s * G_Meijer_1_NLOS_s;

    arg_2_NLOS_s = V_NLOS / (A_0_M1 * h_l_M1 * b_s);
    product_term_2_NLOS_s = ((1-w_s) * zeta^2) / (gamma(a_s) * c_s) * arg_2_NLOS_s^(zeta^2);

    G_Meijer_2_NLOS_s = meijerG(a_par_s,1,[b_par_1_s,0],b_par_2_s,arg_2_NLOS_s^c_s);
    term_2_NLOS_s = product_term_2_NLOS_s * G_Meijer_2_NLOS_s;

    P_out_ana_NLOS_s(i) = term_1_NLOS_s + term_2_NLOS_s;
    
    % ================= Analytical NLOS via IRS =================
    arg_1_NLOS_f = V_NLOS / (A_0_M1 * h_l_M1 * lamb_f);
    product_term_1_NLOS_f = (w_f * zeta^2) * arg_1_NLOS_f^(zeta^2);
    G_Meijer_1_NLOS_f = meijerG(1-zeta^2,1,[1-zeta^2,0],-zeta^2,arg_1_NLOS_f);
    term_1_NLOS_f = product_term_1_NLOS_f * G_Meijer_1_NLOS_f;

    arg_2_NLOS_f = V_NLOS / (A_0_M1 * h_l_M1 * b_f);
    product_term_2_NLOS_f = ((1-w_f) * zeta^2) / (gamma(a_f) * c_f) * arg_2_NLOS_f^(zeta^2);

    G_Meijer_2_NLOS_f = meijerG(a_par_f,1,[b_par_1_f,0],b_par_2_f,arg_2_NLOS_f^c_f);
    term_2_NLOS_f = product_term_2_NLOS_f * G_Meijer_2_NLOS_f;

    P_out_ana_NLOS_f(i) = term_1_NLOS_f + term_2_NLOS_f;

    % Final outage at U_1:
    P_out_ana_auto_s(i) = P_out_ana_LOS_s(i) *  P_out_ana_NLOS_s(i);
    P_out_ana_auto_f(i) = P_out_ana_LOS_f(i) *  P_out_ana_NLOS_f(i);
    
    % ================= Analytical LOS from M_1 to M_2 =================

    V_M2 = gamma_th_1 / (a_2 * rho  * eta);

    arg_1_M2_s = V_M2 / (A_0_2 * h_l_M2 *lamb_s);
    product_term_1_M2_s = (w_s * zeta^2 ) * arg_1_M2_s^(zeta^2);
    G_Meijer_1_M2_s = meijerG(1-zeta^2,1,[1-zeta^2,0],-zeta^2,arg_1_M2_s);
    term_1_M2_s = product_term_1_M2_s * G_Meijer_1_M2_s ;
    
    arg_2_M2_s = V_M2 / (A_0_2 * h_l_M2 * b_s);
    product_term_2_M2_s = ((1-w_s) * zeta^2) / (gamma(a_s) *c_s) * arg_2_M2_s^(zeta^2);
    G_Meijer_2_M2_s = meijerG(a_par_s,1,[b_par_1_s,0],b_par_2_s,arg_2_M2_s^c_s);
    term_2_M2_s = product_term_2_M2_s * G_Meijer_2_M2_s ;

    P_outage_M2_s(i) = term_1_M2_s + term_2_M2_s;

    arg_1_M2_f = V_M2 / (A_0_2 * h_l_M2 *lamb_f);
    product_term_1_M2_f = (w_f * zeta^2 ) * arg_1_M2_f^(zeta^2);
    G_Meijer_1_M2_f = meijerG(1-zeta^2,1,[1-zeta^2,0],-zeta^2,arg_1_M2_f);
    term_1_M2_f = product_term_1_M2_f * G_Meijer_1_M2_f ;
    
    arg_2_M2_f = V_M2 / (A_0_2 * h_l_M2 * b_f);
    product_term_2_M2_f = ((1-w_f) * zeta^2) / (gamma(a_f) *c_f) * arg_2_M2_f^(zeta^2);
    G_Meijer_2_M2_f = meijerG(a_par_f,1,[b_par_1_f,0],b_par_2_f,arg_2_M2_f^c_f);
    term_2_M2_f = product_term_2_M2_f * G_Meijer_2_M2_f ;

    P_outage_M2_f(i) = term_1_M2_f + term_2_M2_f;

    % ================= Analytical LOS from buoy to M_2 =================

    V_B_M2 = gamma_th_2 / ((a_2 - (a_1 *  gamma_th_2)) * rho * eta);

    arg_1_B_M2_s = V_B_M2 / (A_0_B_2 * h_l_B_M2 *lamb_s);
    product_term_1_B_M2_s = (w_s * zeta^2 ) * arg_1_B_M2_s^(zeta^2);
    G_Meijer_1_B_M2_s = meijerG(1-zeta^2,1,[1-zeta^2,0],-zeta^2,arg_1_B_M2_s);
    term_1_B_M2_s = product_term_1_B_M2_s * G_Meijer_1_B_M2_s ;
    
    arg_2_B_M2_s = V_B_M2 / (A_0_B_2 * h_l_B_M2 * b_s);
    product_term_2_B_M2_s = ((1-w_s) * zeta^2) / (gamma(a_s) *c_s) * arg_2_B_M2_s^(zeta^2);
    G_Meijer_2_B_M2_s = meijerG(a_par_s,1,[b_par_1_s,0],b_par_2_s,arg_2_B_M2_s^c_s);
    term_2_B_M2_s = product_term_2_B_M2_s * G_Meijer_2_B_M2_s ;

    P_outage_B_M2_s(i) = term_1_B_M2_s + term_2_B_M2_s;

    P_out_B_M2_s(i) = P_outage_M2_s(i) *  P_outage_B_M2_s(i);

    arg_1_B_M2_f = V_B_M2 / (A_0_B_2 * h_l_B_M2 *lamb_f);
    product_term_1_B_M2_f = (w_f * zeta^2 ) * arg_1_B_M2_f^(zeta^2);
    G_Meijer_1_B_M2_f = meijerG(1-zeta^2,1,[1-zeta^2,0],-zeta^2,arg_1_B_M2_f);
    term_1_B_M2_f = product_term_1_B_M2_f * G_Meijer_1_B_M2_f ;
    
    arg_2_B_M2_f = V_B_M2 / (A_0_B_2 * h_l_B_M2 * b_f);
    product_term_2_B_M2_f = ((1-w_f) * zeta^2) / (gamma(a_f) *c_f) * arg_2_B_M2_f^(zeta^2);
    G_Meijer_2_B_M2_f = meijerG(a_par_f,1,[b_par_1_f,0],b_par_2_f,arg_2_B_M2_f^c_f);
    term_2_B_M2_f = product_term_2_B_M2_f * G_Meijer_2_B_M2_f ;

    P_outage_B_M2_f(i) = term_1_B_M2_f + term_2_B_M2_f;

    P_out_B_M2_f(i) = P_outage_M2_f(i) *  P_outage_B_M2_f(i);

    % ================= Monte Carlo channels =================
    h_t_s = eggrand_salt(N_sim, w_s, lamb_s, a_s, b_s, c_s);
    h_t_f = eggrand_fresh(N_sim, w_f, lamb_f, a_f, b_f, c_f);

    h_p_M1  = pointing_error_model(A_0_M1,  zeta, N_sim);   % NLOS
    h_p_M11 = pointing_error_model(A_0_M11, zeta, N_sim);   % LOS
    h_p_M2 = pointing_error_M2(A_0_2, zeta, N_sim);
    h_p_B_M2 = pointing_error_B_M2(A_0_B_2, zeta, N_sim);

    h_M1_s  = h_l_M1  * h_p_M1  .* h_t_s;     % NLOS via IRS
    h_M11_s = h_l_M11 * h_p_M11 .* h_t_s;     % LOS
    h_2_s = h_l_M2 * h_p_M2 .* h_t_s;
    h_M12_s = h_l_B_M2 * h_p_B_M2 .* h_t_s;

    h_M1_f  = h_l_M1  * h_p_M1  .* h_t_f;     % NLOS via IRS
    h_M11_f = h_l_M11 * h_p_M11 .* h_t_f;     % LOS
    h_2_f = h_l_M2 * h_p_M2 .* h_t_f;
    h_M12_f = h_l_B_M2 * h_p_B_M2 .* h_t_f;

    % LOS SINRs
    gamma_X1_LOS_s = (eta * a_1 * rho * h_M11_s) ./ ((eta * epsilon * a_2 * rho * h_M11_s) + 1);
    gamma_X2_LOS_s = (eta * a_2 * rho * h_M11_s) ./ ((eta * a_1 * rho * h_M11_s) + 1);

    outage_LOS_s = (gamma_X2_LOS_s < gamma_th_2) & (gamma_X1_LOS_s < gamma_th_1);
    P_out_mc_LOS_s(i) = mean(outage_LOS_s);

    gamma_X1_LOS_f = (eta * a_1 * rho * h_M11_f) ./ ((eta * epsilon * a_2 * rho * h_M11_f) + 1);
    gamma_X2_LOS_f = (eta * a_2 * rho * h_M11_f) ./ ((eta * a_1 * rho * h_M11_f) + 1);

    outage_LOS_f = (gamma_X2_LOS_f < gamma_th_2) & (gamma_X1_LOS_f < gamma_th_1);
    P_out_mc_LOS_f(i) = mean(outage_LOS_f);

    % NLOS via IRS SINRs
    gamma_X1_NLOS_s = (eta * a_1 * rho * N^2 * h_M1_s) ./ ((eta * epsilon * a_2 * rho * N^2 * h_M1_s) + 1);
    gamma_X2_NLOS_s = (eta * a_2 * rho * N^2 * h_M1_s) ./ ((eta *  a_1 * rho * N^2 * h_M1_s) + 1);

    outage_NLOS_s = (gamma_X2_NLOS_s < gamma_th_2) & (gamma_X1_NLOS_s < gamma_th_1);
    P_out_mc_NLOS_s(i) = mean(outage_NLOS_s);

    gamma_X1_NLOS_f = (eta * a_1 * rho * N^2 * h_M1_f) ./ ((eta * epsilon * a_2 * rho * N^2 * h_M1_f) + 1);
    gamma_X2_NLOS_f = (eta * a_2 * rho * N^2 * h_M1_f) ./ ((eta *  a_1 * rho * N^2 * h_M1_f) + 1);

    outage_NLOS_f = (gamma_X2_NLOS_f < gamma_th_2) & (gamma_X1_NLOS_f < gamma_th_1);
    P_out_mc_NLOS_f(i) = mean(outage_NLOS_f);

    P_out_mc_auto_s(i) =P_out_mc_NLOS_s(i) * P_out_mc_LOS_s(i);
    P_out_mc_auto_f(i) =P_out_mc_NLOS_f(i) * P_out_mc_LOS_f(i);
    
    % LOS SINRs of M_2
    gamma_X2_L2_s = eta * a_2 * rho * h_2_s ;   % LOS from M_1 to M_2
    gamma_B_X2_LOS_s = (eta * a_2 * rho * h_M12_s) ./ ((eta * a_1 * rho * h_M12_s) + 1);  % LOS from buoy to M_2

    outage_M2_s = gamma_X2_L2_s < gamma_th_2;
    P_out_mc_M2_s(i) = mean(outage_M2_s);

    outage_B_M2_s = gamma_B_X2_LOS_s < gamma_th_2;
    P_out_mc_B_M2_s(i) = mean(outage_B_M2_s);

    P_outage_mc_B_M2_s(i) =P_out_mc_M2_s(i) * P_out_mc_B_M2_s(i);

    gamma_X2_L2_f = eta * a_2 * rho * h_2_f ;   % LOS from M_1 to M_2
    gamma_B_X2_LOS_f = (eta * a_2 * rho * h_M12_f) ./ ((eta * a_1 * rho * h_M12_f) + 1);  % LOS from buoy to M_2

    outage_M2_f = gamma_X2_L2_f < gamma_th_2;
    P_out_mc_M2_f(i) = mean(outage_M2_f);

    outage_B_M2_f = gamma_B_X2_LOS_f < gamma_th_2;
    P_out_mc_B_M2_f(i) = mean(outage_B_M2_f);

    P_outage_mc_B_M2_f(i) =P_out_mc_M2_f(i) * P_out_mc_B_M2_f(i);

end

% ================= Plot =================
figure;

semilogy(gamma_dB, P_out_ana_auto_s, 'b-',  'LineWidth', 2);hold on;
semilogy(gamma_dB, P_out_mc_auto_s,  'bo',  'LineWidth', 1.5);
semilogy(gamma_dB, P_out_ana_auto_f, 'r-',  'LineWidth', 2);hold on;
semilogy(gamma_dB, P_out_mc_auto_f,  'ro',  'LineWidth', 1.5);
semilogy(gamma_dB, P_out_B_M2_s, 'b--', 'LineWidth', 1.5);
semilogy(gamma_dB, P_outage_mc_B_M2_s, 'bo', 'LineWidth', 1.5);
semilogy(gamma_dB, P_out_B_M2_f, 'r--', 'LineWidth', 1.5);
semilogy(gamma_dB, P_outage_mc_B_M2_f, 'ro', 'LineWidth', 1.5);


xlabel('SNR (dB)');
ylabel('Outage Probability');
legend('Analytical for U_1 in Salt water B.L=16.5L/min', 'Simulation for U_1 in Salt water B.L=16.5L/min', ...
    'Analytical for U_1 in fresh water B.L=4.7L/min', 'Simulation for U_1 in fresh water B.L=4.7L/min' , ...
    'Analytical for U_2 in Salt water B.L=16.5L/min', 'Simulation for U_2 in Salt water B.L=16.5L/min', ...
    'Analytical for U_2 in fresh water B.L=4.7L/min', 'Simulation for U_2 in fresh water B.L=4.7L/min','Location', 'southwest');
grid on;
xlim([0 70]);


% ================= Functions =================
function h_t = eggrand_salt(N_sim, w_s, lamb_s, a_s, b_s, c_s)
    U = rand(N_sim,1);
    h_t = zeros(N_sim,1);

    % Exponential branch
    idx1 = (U <= w_s);
    h_t(idx1) = exprnd(lamb_s, sum(idx1),1);

    % Generalized Gamma branch
    idx2 = ~idx1;
    Y = gamrnd(a_s,1,sum(idx2),1);   % Gamma(a,1)
    h_t(idx2) = b_s * (Y.^(1/c_s));  % GG transformation
end
function h_t = eggrand_fresh(N_sim, w_f, lamb_f, a_f, b_f, c_f)
    U = rand(N_sim,1);
    h_t = zeros(N_sim,1);

    % Exponential branch
    idx1 = (U <= w_f);
    h_t(idx1) = exprnd(lamb_f, sum(idx1),1);

    % Generalized Gamma branch
    idx2 = ~idx1;
    Y = gamrnd(a_f,1,sum(idx2),1);   % Gamma(a,1)
    h_t(idx2) = b_f * (Y.^(1/c_f));  % GG transformation
end


function h_p = pointing_error_model(A_0, zeta, N_sim)
    U = rand(N_sim,1);
    R = sqrt(-2 * log(U));
    h_p = A_0 * exp(-(R.^2) ./ (2 * zeta^2));
end

% Pointing error random generation
function h_p_2 = pointing_error_M2(A_0_2, zeta, N_sim)
    U = rand(N_sim,1);        % Uniform(0,1)
    h_p_2 = A_0_2 * (U.^(1/zeta^2));
end

% Pointing error random generation
function h_p_B_M2 = pointing_error_B_M2(A_0_B_2, zeta, N_sim)
    U = rand(N_sim,1);        % Uniform(0,1)
    h_p_B_M2 = A_0_B_2 * (U.^(1/zeta^2));
end