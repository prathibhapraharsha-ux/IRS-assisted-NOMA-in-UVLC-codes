clc
clear
close all

% Parameters
eta = 0.9;          % responsivity  
a_1 = 0.3;          % Power allocation coefficient for near user
a_2 = 0.7;          % Power allocation coefficient for far user
epsilon = 0.01;     % residual  interference

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

gamma_th_1 = 2^R_s - 1 ;
gamma_th_2 = gamma_th_1;
gamma_th = 10^(0.3);

gamma_th_dB = [0.1 1.5 3.5];
gamma_th_val = 10.^(gamma_th_dB / 10);

% SNR range
gamma_dB = 0:5:80;
gamma_values = 10.^(gamma_dB/10);

length_user = 0:20:180;

% Output arrays
P_out_ana_NLOS = zeros(size(length_user));
P_out_mc_NLOS  = zeros(size(length_user));

N_sim = 1e5;

for j = 1:length(gamma_th_val)
    rho = gamma_th_val(j);

    for i = 1:length(length_user)
        L_M1 = length_user(i);
        
        h_l_M1  = exp(-c_lambda * L_M1);    % NLOS via IRS

        % Aperture / beam spread for NLOS via IRS
        W_Z_M1 = W_0 * sqrt(1 + ((lambda * L_M1) / (pi * W_0^2))^2);
        v_1_M1 = sqrt(pi/2) * (D_R / (2 * W_Z_M1));
        A_0_M1 = (erf(v_1_M1))^2;

        
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
        a_par = 1-zeta^2/c;
        b_par_1 = a-(zeta^2 /c);
        b_par_2 = -zeta^2/c ;
        G_Meijer_2_NLOS = meijerG(a_par,1,[b_par_1,0],b_par_2,arg_2_NLOS^c);
        term_2_NLOS = product_term_2_NLOS * G_Meijer_2_NLOS;
    
        P_out_ana_NLOS(j,i) = term_1_NLOS + term_2_NLOS;
    
        % ================= Monte Carlo channels =================
        h_t = eggrand(N_sim, w, lamb, a, b, c);
        h_p_M1  = pointing_error_model(A_0_M1,  zeta, N_sim);   % NLOS
        h_M1  = h_l_M1  * h_p_M1  .* h_t;     % NLOS via IRS
     
 
        % NLOS via IRS SINRs
        gamma_X1_NLOS = (eta * a_1 * rho * N^2 * h_M1) ./ ((eta * epsilon * a_2 * rho * N^2 * h_M1) + 1);
        gamma_X2_NLOS = (eta * a_2 * rho * N^2 * h_M1) ./ ((eta *  a_1 * rho * N^2 * h_M1) + 1);
    
        outage_NLOS = (gamma_X2_NLOS < gamma_th_1) & (gamma_X1_NLOS < gamma_th_1);
        P_out_mc_NLOS(j,i) = mean(outage_NLOS);
       
    end
end

% ================= Plot =================
figure;

semilogy(length_user, P_out_ana_NLOS(1,:), 'b-',  'LineWidth', 2);hold on;
semilogy(length_user, P_out_mc_NLOS(1,:),  'ro',  'LineWidth', 1.5);
semilogy(length_user, P_out_ana_NLOS(2,:), 'b-',  'LineWidth', 2);
semilogy(length_user, P_out_mc_NLOS(2,:),  'r*',  'LineWidth', 1.5);
semilogy(length_user, P_out_ana_NLOS(3,:), 'b-',  'LineWidth', 2);
semilogy(length_user, P_out_mc_NLOS(3,:),  'r^',  'LineWidth', 1.5);

xlabel('User U_1 NLOS distance(d_1 + d_2)');
ylabel('Outage Probability');
legend('Analytical of U_1 at \gamma_{th}=1dB','Simulation of U_1 at \gamma_{th}=1dB', ...
    'Analytical of U_1 at \gamma_{th}=2dB','Simulation of U_1 at \gamma_{th}=2dB', ...
    'Analytical of U_1 at \gamma_{th}=3dB','Simulation of U_1 at \gamma_{th}=3dB', ...
    'Location', 'southwest');
grid on;
%xlim([20 60]);


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