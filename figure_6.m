clc
clear
close all

% Parameters
eta = 0.9;          % responsivity  
a_1 = 0.3;          % Power allocation coefficient for near user
a_2 = 0.7;          % Power allocation coefficient for far user
epsilon = 0.01;     % residual  interference


L_1 = 60;           % Link length from S to IRS
L_2 = 60;           % Link length from IRS to U_1
L_M1 = L_1 + L_2;   % NLOS path of U_1 through IRS
L_M2 = 40;          % LOS path from U_1 to U_2
L_M11 = 70;         % LOS path length of U_1
L_B_M2 = 80;        % LOS from Buoy to U_2

c_lambda = 0.056;   % extinction coefficient for pure sea water
W_0 = 0.01;         % spot size of gaussian beam
lambda = 470e-9;    % wavelength of light
D_R = 0.05;         % receiver aperture diameter
R_s = 1;            % data rate
zeta = 1;           % pointing error coefficient
N = 64;             % number of reflecting elements in IRS

% Parameters for B.L = 4.7L/min
w = 0.4589;
lamb = 0.3449;
a = 1.0421;
b = 1.5768;
c = 35.9424;

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

gamma_th_dB = [1 3];
gamma_th_val = 10.^(gamma_th_dB / 10);

% SNR range
gamma_dB = 0:5:80;
gamma_values = 10.^(gamma_dB/10);

% Output arrays
P_out_ana_LOS  = zeros(size(gamma_dB));
P_out_ana_NLOS = zeros(size(gamma_dB));
P_outage_M2 = zeros(size(gamma_dB));
P_outage_B_M2 = zeros(size(gamma_dB));
P_out_ana_auto = zeros(size(gamma_dB));
P_out_B_M2 = zeros(size(gamma_dB));

P_out_mc_LOS   = zeros(size(gamma_dB));
P_out_mc_NLOS  = zeros(size(gamma_dB));
P_out_mc_M2 = zeros(size(gamma_dB));
P_out_mc_auto  = zeros(size(gamma_dB));
P_out_mc_B_M2 = zeros(size(gamma_dB));
P_outage_mc_B_M2 = zeros(size(gamma_dB));

N_sim = 1e5;

for j = 1:length(gamma_th_val)
    gamma_th = gamma_th_val(j);

    for i = 1:length(gamma_values)
        rho = gamma_values(i);
    
        % ================= Analytical LOS =================
        alpha_LOS = gamma_th / ((a_2 - (a_1 *  gamma_th)) * rho * eta);
        beta_LOS  = gamma_th / ((a_1 - (a_2 * epsilon * gamma_th)) * rho * eta);
        V_LOS = min(alpha_LOS, beta_LOS);
    
        arg_1_LOS = V_LOS / (A_0_M11 * h_l_M11 * lamb);
        product_term_1_LOS = (w * zeta^2) * arg_1_LOS^(zeta^2);
        G_Meijer_1_LOS = meijerG(1-zeta^2,1,[1-zeta^2,0],-zeta^2,arg_1_LOS);
        term_1_LOS = product_term_1_LOS * G_Meijer_1_LOS;
    
        arg_2_LOS = V_LOS / (A_0_M11 * h_l_M11 * b);
        product_term_2_LOS = ((1-w) * zeta^2) / (gamma(a) * c) * arg_2_LOS^(zeta^2);
    
        a_par = 1 - zeta^2/c;
        b_par_1 = a - (zeta^2/c);
        b_par_2 = -zeta^2/c;
    
        G_Meijer_2_LOS = meijerG(a_par,1,[b_par_1,0],b_par_2,arg_2_LOS^c);
        term_2_LOS = product_term_2_LOS * G_Meijer_2_LOS;
    
        P_out_ana_LOS(j,i) = term_1_LOS + term_2_LOS;
    
        % ================= Analytical NLOS via IRS =================
        alpha_NLOS = gamma_th / ((a_2 - (a_1 *  gamma_th)) * rho * N^2 * eta);
        beta_NLOS  = gamma_th / ((a_1 - (a_2 * epsilon * gamma_th)) * rho * N^2 * eta);
        V_NLOS = min(alpha_NLOS, beta_NLOS);
    
        arg_1_NLOS = V_NLOS / (A_0_M1 * h_l_M1 * lamb);
        product_term_1_NLOS = (w * zeta^2) * arg_1_NLOS^(zeta^2);
        G_Meijer_1_NLOS = meijerG(1-zeta^2,1,[1-zeta^2,0],-zeta^2,arg_1_NLOS);
        term_1_NLOS = product_term_1_NLOS * G_Meijer_1_NLOS;
    
        arg_2_NLOS = V_NLOS / (A_0_M1 * h_l_M1 * b);
        product_term_2_NLOS = ((1-w) * zeta^2) / (gamma(a) * c) * arg_2_NLOS^(zeta^2);
    
        G_Meijer_2_NLOS = meijerG(a_par,1,[b_par_1,0],b_par_2,arg_2_NLOS^c);
        term_2_NLOS = product_term_2_NLOS * G_Meijer_2_NLOS;
    
        P_out_ana_NLOS(j,i) = term_1_NLOS + term_2_NLOS;
    
        % Final outage at U_1:
        P_out_ana_auto(j,i) = P_out_ana_LOS(j,i) *  P_out_ana_NLOS(j,i);
        
        % ================= Analytical LOS from M_1 to M_2 =================
    
        V_M2 = gamma_th / (a_2 * rho  * eta);
    
        arg_1_M2 = V_M2 / (A_0_2 * h_l_M2 *lamb);
        product_term_1_M2 = (w * zeta^2 ) * arg_1_M2^(zeta^2);
        G_Meijer_1_M2 = meijerG(1-zeta^2,1,[1-zeta^2,0],-zeta^2,arg_1_M2);
        term_1_M2 = product_term_1_M2 * G_Meijer_1_M2 ;
        
        arg_2_M2 = V_M2 / (A_0_2 * h_l_M2 * b);
        product_term_2_M2 = ((1-w) * zeta^2) / (gamma(a) *c) * arg_2_M2^(zeta^2);
        a_par = 1-zeta^2/c;
        b_par_1 = a-(zeta^2 /c);
        b_par_2 = -zeta^2/c ;
        G_Meijer_2_M2 = meijerG(a_par,1,[b_par_1,0],b_par_2,arg_2_M2^c);
        term_2_M2 = product_term_2_M2 * G_Meijer_2_M2 ;
    
        P_outage_M2(j,i) = term_1_M2 + term_2_M2;
    
        % ================= Analytical LOS from buoy to M_2 =================
    
        V_B_M2 = gamma_th / ((a_2 - (a_1 *  gamma_th)) * rho * eta);
    
        arg_1_B_M2 = V_B_M2 / (A_0_B_2 * h_l_B_M2 *lamb);
        product_term_1_B_M2 = (w * zeta^2 ) * arg_1_B_M2^(zeta^2);
        G_Meijer_1_B_M2 = meijerG(1-zeta^2,1,[1-zeta^2,0],-zeta^2,arg_1_B_M2);
        term_1_B_M2 = product_term_1_B_M2 * G_Meijer_1_B_M2 ;
        
        arg_2_B_M2 = V_B_M2 / (A_0_B_2 * h_l_B_M2 * b);
        product_term_2_B_M2 = ((1-w) * zeta^2) / (gamma(a) *c) * arg_2_B_M2^(zeta^2);
        a_par = 1-zeta^2/c;
        b_par_1 = a-(zeta^2 /c);
        b_par_2 = -zeta^2/c ;
        G_Meijer_2_B_M2 = meijerG(a_par,1,[b_par_1,0],b_par_2,arg_2_B_M2^c);
        term_2_B_M2 = product_term_2_B_M2 * G_Meijer_2_B_M2 ;
    
        P_outage_B_M2(j,i) = term_1_B_M2 + term_2_B_M2;
    
        P_out_B_M2(j,i) = P_outage_M2(j,i) *  P_outage_B_M2(j,i);
    
        % ================= Monte Carlo channels =================
        h_t = eggrand(N_sim, w, lamb, a, b, c);
    
        h_p_M1  = pointing_error_model(A_0_M1,  zeta, N_sim);   % NLOS
        h_p_M11 = pointing_error_model(A_0_M11, zeta, N_sim);   % LOS
        h_p_M2 = pointing_error_M2(A_0_2, zeta, N_sim);
        h_p_B_M2 = pointing_error_B_M2(A_0_B_2, zeta, N_sim);
    
        h_M1  = h_l_M1  * h_p_M1  .* h_t;     % NLOS via IRS
        h_M11 = h_l_M11 * h_p_M11 .* h_t;     % LOS
        h_2 = h_l_M2 * h_p_M2 .* h_t;
        h_M12 = h_l_B_M2 * h_p_B_M2 .* h_t;
    
        % LOS SINRs
        gamma_X1_LOS = (eta * a_1 * rho * h_M11) ./ ((eta * epsilon * a_2 * rho * h_M11) + 1);
        gamma_X2_LOS = (eta * a_2 * rho * h_M11) ./ ((eta * a_1 * rho * h_M11) + 1);
    
        outage_LOS = (gamma_X2_LOS < gamma_th) & (gamma_X1_LOS < gamma_th);
        P_out_mc_LOS(j,i) = mean(outage_LOS);
    
        % NLOS via IRS SINRs
        gamma_X1_NLOS = (eta * a_1 * rho * N^2 * h_M1) ./ ((eta * epsilon * a_2 * rho * N^2 * h_M1) + 1);
        gamma_X2_NLOS = (eta * a_2 * rho * N^2 * h_M1) ./ ((eta *  a_1 * rho * N^2 * h_M1) + 1);
    
        outage_NLOS = (gamma_X2_NLOS < gamma_th) & (gamma_X1_NLOS < gamma_th);
        P_out_mc_NLOS(j,i) = mean(outage_NLOS);
        
        % LOS SINRs of M_2
        gamma_X2_L2 = eta * a_2 * rho * h_2 ;   % LOS from M_1 to M_2
        gamma_B_X2_LOS = (eta * a_2 * rho * h_M12) ./ ((eta * a_1 * rho * h_M12) + 1);  % LOS from buoy to M_2
    
    
        P_out_mc_auto(j,i) =P_out_mc_NLOS(j,i) * P_out_mc_LOS(j,i);
    
        outage_M2 = gamma_X2_L2 < gamma_th;
        P_out_mc_M2(j,i) = mean(outage_M2);
    
        outage_B_M2 = gamma_B_X2_LOS < gamma_th;
        P_out_mc_B_M2(j,i) = mean(outage_B_M2);
    
        P_outage_mc_B_M2(j,i) =P_out_mc_M2(j,i) * P_out_mc_B_M2(j,i);
    
    end
end

% ================= Plot =================
figure;

semilogy(gamma_dB, P_out_ana_auto(1,:), 'b-',  'LineWidth', 2);hold on;
semilogy(gamma_dB, P_out_mc_auto(1,:),  'ro',  'LineWidth', 1.5);
semilogy(gamma_dB, P_out_B_M2(1,:), 'b--', 'LineWidth', 1.5);
semilogy(gamma_dB, P_outage_mc_B_M2(1,:), 'ro', 'LineWidth', 1.5);

semilogy(gamma_dB, P_out_ana_auto(2,:), 'b-',  'LineWidth', 2);
semilogy(gamma_dB, P_out_mc_auto(2,:),  'r*',  'LineWidth', 1.5);
semilogy(gamma_dB, P_out_B_M2(2,:), 'b--', 'LineWidth', 1.5);
semilogy(gamma_dB, P_outage_mc_B_M2(2,:), 'r*', 'LineWidth', 1.5);



xlabel('SNR (dB)');
ylabel('Outage Probability');
legend('Analytical of U_1 at \gamma_{th}=1dB','Simulation of U_1 at \gamma_{th}=1dB', ... 
    'Analytical of U_2 at \gamma_{th}=1dB','Simulation of U_2 at \gamma_{th}=1dB', ... 
    'Analytical of U_1 at \gamma_{th}=3dB','Simulation of U_1 at \gamma_{th}=3dB', ... 
    'Analytical of U_2 at \gamma_{th}=3dB','Simulation of U_2at \gamma_{th}=3dB', ... 
    'Location', 'southwest');
grid on;
xlim([20 60]);


% ================= Functions =================
function h_t = eggrand(N_sim, w, lamb, a, b, c)
    U = rand(N_sim,1);
    h_t = zeros(N_sim,1);

    % Exponential branch
    idx1 = (U <= w);
    h_t(idx1) = exprnd(lamb, sum(idx1), 1);

    % Generalized Gamma branch
    idx2 = ~idx1;
    Y = gamrnd(a, 1, sum(idx2), 1);
    h_t(idx2) = b * (Y.^(1/c));
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