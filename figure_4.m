clc
clear
close all

% Parameters
eta = 0.9;
epsilon = 0.01;

L_1 = 60;           % Link length from S to IRS
L_2 = 60;           % Link length from IRS to M_1
L_M1 = L_1 + L_2;   % NLOS through IRS total path
L_M2 = 40;          % LOS path from M_1 to M_2
L_M11 = 70;         % LOS path length
L_B_M2 = 80;        % LOS from Buoy to M_2

c_lambda = 0.056;   % extinction coefficient for pure sea water
W_0 = 0.01;
lambda = 470e-9;
D_R = 0.05;         % receiver aperture diameter
R_s = 1;
zeta = 1;
N = 16;

gamma_th_1 = 2^R_s - 1;
gamma_th_2 = gamma_th_1;

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

% SNR range
gamma_dB = 0:10:80;
gamma_values = 10.^(gamma_dB/10);

a_2_val = [0.6 0.8];
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
for j = 1:length(a_2_val)
    a_2 = a_2_val(j);
    a_1 = 1- a_2;

    for i = 1:length(gamma_values)
        rho = gamma_values(i);
    
        % ================= Analytical LOS =================
        alpha_LOS = gamma_th_2 / ((a_2 - (a_1 *  gamma_th_2)) * rho * eta);
        beta_LOS  = gamma_th_1 / ((a_1 - (a_2 * epsilon * gamma_th_1)) * rho * eta);
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
        alpha_NLOS = gamma_th_2 / ((a_2 - (a_1 *  gamma_th_2)) * rho * N^2 * eta);
        beta_NLOS  = gamma_th_1 / ((a_1 - (a_2 * epsilon * gamma_th_1)) * rho * N^2 * eta);
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
    
        V_M2 = gamma_th_1 / (a_2 * rho  * eta);
    
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
    
        P_outage_M2(i) = term_1_M2 + term_2_M2;
    
        % ================= Analytical LOS from buoy to M_2 =================
    
        V_B_M2 = gamma_th_2 / ((a_2 - (a_1 *  gamma_th_2)) * rho * eta);
    
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
    
        P_outage_B_M2(i) = term_1_B_M2 + term_2_B_M2;
    
        P_out_B_M2(j,i) = P_outage_M2(i) *  P_outage_B_M2(i);
    
        %P_out_B_M2(i) = (P_out_ana_LOS(i) * P_out_ana_NLOS(i) * (1+P_outage_M2(i)) + P_outage_M2(i)) *  P_outage_B_M2(i);
    
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
    
        outage_LOS = (gamma_X2_LOS < gamma_th_2) & (gamma_X1_LOS < gamma_th_1);
        P_out_mc_LOS(j,i) = mean(outage_LOS);
    
        % NLOS via IRS SINRs
        gamma_X1_NLOS = (eta * a_1 * rho * N^2 * h_M1) ./ ((eta * epsilon * a_2 * rho * N^2 * h_M1) + 1);
        gamma_X2_NLOS = (eta * a_2 * rho * N^2 * h_M1) ./ ((eta *  a_1 * rho * N^2 * h_M1) + 1);
    
        outage_NLOS = (gamma_X2_NLOS < gamma_th_2) & (gamma_X1_NLOS < gamma_th_1);
        P_out_mc_NLOS(j,i) = mean(outage_NLOS);
        
        % LOS SINRs of M_2
        gamma_X2_L2 = eta * a_2 * rho * h_2 ;   % LOS from M_1 to M_2
        gamma_B_X2_LOS = (eta * a_2 * rho * h_M12) ./ ((eta * a_1 * rho * h_M12) + 1);  % LOS from buoy to M_2
    
        gamma_X2_U1 = min(gamma_X2_NLOS, gamma_X2_LOS);
    
        P_out_mc_auto(j,i) =P_out_mc_NLOS(j,i) * P_out_mc_LOS(j,i);
    
        outage_M2 = gamma_X2_L2 < gamma_th_2;
        P_out_mc_M2(i) = mean(outage_M2);
    
        outage_B_M2 = gamma_B_X2_LOS < gamma_th_2;
        P_out_mc_B_M2(i) = mean(outage_B_M2);
    
        P_outage_mc_B_M2(j,i) =P_out_mc_M2(i) * P_out_mc_B_M2(i);
    end
end

% ================= Plot =================
figure;
semilogy(gamma_dB, P_out_ana_LOS(1,:), 'g-', 'LineWidth', 1.5); hold on;
semilogy(gamma_dB, P_out_mc_LOS(1,:), 'gp', 'LineWidth', 1.5);
semilogy(gamma_dB, P_out_ana_LOS(2,:), 'r--', 'LineWidth', 1.5);
semilogy(gamma_dB, P_out_mc_LOS(2,:), 'ro', 'LineWidth', 1.5);
semilogy(gamma_dB, P_out_B_M2(1,:), 'b-', 'LineWidth', 1.5); 
semilogy(gamma_dB, P_outage_mc_B_M2(1,:), 'bp', 'LineWidth', 1.5);
semilogy(gamma_dB, P_out_B_M2(2,:), 'b-', 'LineWidth', 1.5);
semilogy(gamma_dB, P_outage_mc_B_M2(2,:), 'bo', 'LineWidth', 1.5);

%%% for N= 32;
P_out_ana_auto_val_1 = [0.999927196046063	0.423787203237010	0.099285158592416	0.005780399001969	0.000179048885521	0.000003657731868	0.000000061411895	0.000000000923464	0.000000000012937];
P_out_mc_auto_val_1 = [0.999950000000000	0.424220000000000	0.098187574100000	0.005813408400000	0.000172184500000	0.000004182500000	0.000000032400000	0	0];

P_out_ana_auto_val_2 = [0.977865207486988	0.343130890194672	0.059532935493474	0.003217744561853	0.000089893892604	0.000001753172150	0.000000028721835	0.000000000425410	0.000000000005899];
P_out_mc_auto_val_2 = [0.978200000000000	0.342280000000000	0.059439087400000	0.003221015100000	0.000087220800000	0.000001533300000	0.000000033600000	0.000000001200000	0];

%%%%% for N=64
P_out_ana_auto_val_1_64 = [0.630443545802800	0.194685287668375	0.034668563777668	0.001795577756647	0.000052512307971	0.000001037414107	0.000000017040785	0.000000000252327	0.000000000003495];
P_out_mc_auto_val_1_64 = [0.631500000000000	0.193950000000000	0.034589813500000	0.001731071100000	0.000043895600000	0.000000821800000	0.000000016700000	0	0];

P_out_ana_auto_val_2_64 = [0.529375372144143	0.147459235419814	0.020190865332634	0.000986286542880	0.000026170565776	0.000000494928185	0.000000007944357	0.000000000115968	0.000000000001591];
P_out_mc_auto_val_2_64 = [0.531020000000000	0.146970000000000	0.020137228000000	0.001000010700000	0.000028857600000	0.000000249600000	0	0	0];

plot(gamma_dB,P_out_ana_auto_val_1,'k-','LineWidth',2);
plot(gamma_dB, P_out_mc_auto_val_1,'ko');
plot(gamma_dB,P_out_ana_auto_val_2,'k-','LineWidth',2 );
plot(gamma_dB, P_out_mc_auto_val_2,'k*');

plot(gamma_dB,P_out_ana_auto_val_1_64 ,'r-','LineWidth',2);
plot(gamma_dB,P_out_mc_auto_val_1_64,'ro' );
plot(gamma_dB,P_out_ana_auto_val_2_64,'r-','LineWidth',2  );
plot(gamma_dB,P_out_mc_auto_val_2_64 ,'r*');


xlabel('SNR (dB)');
ylabel('Outage Probability of U_1');
legend('U_1(No IRS) Ana. ,a_2=0.6', ...
    'U_1(No IRS) sim. ,a_2=0.6', ...
    'U_1(No IRS) Ana. ,a_2=0.8', ...
    'U_1(No IRS) sim. ,a_2=0.8', ...
    'U_2, Ana., a_2=0.6', ...
       'U_2, Sim., a_2=0.6', ...
       'U_2, Ana., a_2=0.8', ...
       'U_2, Sim., a_2=0.8', ...
       'U_1 (IRS), Ana., a_2=0.6 and N=32', ...
       'U_1 (IRS), Sim., a_2=0.6 and N=32', ...
       'U_1 (IRS), Ana., a_2=0.8 and N=32' , ...
       'U_1 (IRS), Sim., a_2=0.8 and N=32', ...
       'U_1 (IRS), Ana., a_2=0.6 and N=64', ...
       'U_1 (IRS), Sim., a_2=0.6 and N=64', ...
       'U_1 (IRS), Ana., a_2=0.8 and N=64' , ...
       'U_1 (IRS), Sim., a_2=0.8 and N=64', ...
       'Location', 'southwest');
grid on;
xlim([0 70]);

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