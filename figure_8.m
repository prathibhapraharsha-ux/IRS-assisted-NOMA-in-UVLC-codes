clc
clear
close all

% Parameters given
eta = 0.9;
N = 64;
a_1 = 0.3;
a_2 = 0.7;
epsilon = 0.01;
L_1 = 40;     % Link length from S to IRS
L_2 = 40;     % Link length from IRS to M_1
L_IRS = L_1 + L_2 ;
L_M1 =60;
L_M2N = 40;
L_M2 = 70;
c_lambda = 0.056;    % extrinsic coefficient for pure sea water
W_0 = 0.01;
lambda = 470 .* 10.^(-9);
D_R = 0.05;     % receiver apeature diamater
R = 0.01;       % random radial displacement
R_s = 0.5;
zeta =1;
gamma_th_1 = 2^R_s - 1 ;
gamma_th_2 = gamma_th_1;


% Parameters for B.L = 4.7L/min
w = 0.4589;
lamb = 0.3449;
a = 1.0421;
b = 1.5768;
c = 35.9424;

h_l_IRS = exp(-c_lambda * L_IRS );
h_l_M1 = exp(-c_lambda * L_M1 );
h_l_M2 = exp(-c_lambda * L_M2 );
h_l_M2N = exp(-c_lambda * L_M2N );

W_Z_IRS = W_0 * (1 + ((lambda * L_IRS) / (pi * W_0^2))^2)^(1/2) ;
v_IRS = sqrt(pi/2) * (D_R / (2 * W_Z_IRS));
A_0_IRS = (erf(v_IRS))^2;

W_Z_1 = W_0 * (1 + ((lambda * L_M1) / (pi * W_0^2))^2)^(1/2) ;
v_1 = sqrt(pi/2) * (D_R / (2 * W_Z_1));
A_0 = (erf(v_1))^2;

W_Z_2 = W_0 * sqrt(1 + ((lambda * L_M2) / (pi * W_0^2))^2);
v_2 = sqrt(pi/2) * (D_R / (2 * W_Z_2));
A_0_2 = (erf(v_2))^2;

W_Z_2N = W_0 * sqrt(1 + ((lambda * L_M2N) / (pi * W_0^2))^2);
v_2N = sqrt(pi/2) * (D_R / (2 * W_Z_2N));
A_0_2N = (erf(v_2N))^2;

gamma_dB = 0:5:60; % Example range, adjust as needed
gamma_values = 10.^(gamma_dB / 10); % SNR in linear scale

P_outage_M1 = zeros(size(gamma_dB));
P_outage_M2 = zeros(size(gamma_dB));
P_outage_M2N = zeros(size(gamma_dB));
P_outage_IRS = zeros(size(gamma_dB));

% Monte Carlo settings
N_sim = 1e5;   % number of realizations
P_out_MC = zeros(size(gamma_values));
P_out_MC_IRS = zeros(size(gamma_values));
P_out_MC_M2 = zeros(size(gamma_values));
P_out_MC_M2N = zeros(size(gamma_values));

for i = 1:length(gamma_values)
    rho = gamma_values(i);
    
    %%%%  Analytical Expression for L_1 IRS %%%%%%%%%%%

    alpha_IRS = gamma_th_2 / ((a_2 - (a_1 * epsilon * gamma_th_2) ) * rho * N^2 * eta);
    beta_IRS = gamma_th_1 / ((a_1 - (a_2 * epsilon * gamma_th_1) ) * rho * N^2 * eta);
    VM1_IRS = min(alpha_IRS,beta_IRS);
   
    arg_1_IRS = VM1_IRS / (A_0_IRS * h_l_IRS *lamb);
    product_term_1_IRS = (w * zeta^2 ) * arg_1_IRS^(zeta^2);
    G_Meijer_1_IRS = meijerG(1-zeta^2,1,[1-zeta^2,0],-zeta^2,arg_1_IRS);
    term_1_IRS = product_term_1_IRS * G_Meijer_1_IRS ;
    
    arg_2_IRS = VM1_IRS / (A_0_IRS * h_l_IRS * b);
    product_term_2_IRS = ((1-w) * zeta^2) / (gamma(a) *c) * arg_2_IRS^(zeta^2);
    a_par = 1-zeta^2/c;
    b_par_1 = a-(zeta^2 /c);
    b_par_2 = -zeta^2/c ;
    G_Meijer_2_IRS = meijerG(a_par,1,[b_par_1,0],b_par_2,arg_2_IRS^c);
    term_2_IRS = product_term_2_IRS * G_Meijer_2_IRS ;

    P_outage_IRS(i) = term_1_IRS + term_2_IRS;
    
    %%%%  Analytical Expression for L_1 without IRS %%%%%%%%%%%

    alpha = gamma_th_2 / ((a_2 - (a_1 * epsilon * gamma_th_2) ) * rho  * eta);
    beta = gamma_th_1 / ((a_1 - (a_2 * epsilon * gamma_th_1) ) * rho  * eta);
    V_M1 = min (alpha,beta);

    arg_1_M1 = V_M1 / (A_0 * h_l_M1 *lamb);
    product_term_1_M1 = (w * zeta^2 ) * arg_1_M1^(zeta^2);
    G_Meijer_1_M1 = meijerG(1-zeta^2,1,[1-zeta^2,0],-zeta^2,arg_1_M1);
    term_1_M1 = product_term_1_M1 * G_Meijer_1_M1 ;
    
    arg_2_M1 = V_M1 / (A_0 * h_l_M1 * b);
    product_term_2_M1 = ((1-w) * zeta^2) / (gamma(a) *c) * arg_2_M1^(zeta^2);
    G_Meijer_2_M1 = meijerG(a_par,1,[b_par_1,0],b_par_2,arg_2_M1^c);
    term_2_M1 = product_term_2_M1 * G_Meijer_2_M1 ;

    P_outage_M1(i) = term_1_M1 + term_2_M1;

    %%%%  Analytical Expression for L_2 with NOMA %%%%%%%%%%%

    V_M2N = gamma_th_1 / (a_2 * rho  * eta);

    arg_1_M2N = V_M2N / (A_0_2N * h_l_M2N *lamb);
    product_term_1_M2N = (w * zeta^2 ) * arg_1_M2N^(zeta^2);
    G_Meijer_1_M2N = meijerG(1-zeta^2,1,[1-zeta^2,0],-zeta^2,arg_1_M2N);
    term_1_M2N = product_term_1_M2N * G_Meijer_1_M2N ;
    
    arg_2_M2N = V_M2N / (A_0_2N * h_l_M2N * b);
    product_term_2_M2N = ((1-w) * zeta^2) / (gamma(a) *c) * arg_2_M2N^(zeta^2);
    G_Meijer_2_M2N = meijerG(a_par,1,[b_par_1,0],b_par_2,arg_2_M2N^c);
    term_2_M2N = product_term_2_M2N * G_Meijer_2_M2N ;

    P_outage_M2N(i) = term_1_M2N + term_2_M2N;

    %%%%  Analytical Expression for L_2 without NOMA %%%%%%%%%%%

    alpha_M2 = gamma_th_2 / ((a_2 - (a_1 *  epsilon * gamma_th_2) ) * rho * eta);
    
    arg_1_M2 = alpha_M2 / (A_0_2 * h_l_M2 *lamb);
    product_term_1_M2 = (w * zeta^2 ) * arg_1_M2^(zeta^2);
    G_Meijer_1_M2 = meijerG(1-zeta^2,1,[1-zeta^2,0],-zeta^2,arg_1_M1);
    term_1_M2 = product_term_1_M2 * G_Meijer_1_M2 ;

    arg_2_M2 = alpha_M2 / (A_0_2 * h_l_M2 * b);
    product_term_2_M2 = ((1-w) * zeta^2) / (gamma(a) *c) * arg_2_M2^(zeta^2);
    G_Meijer_2_M2 = meijerG(a_par,1,[b_par_1,0],b_par_2,arg_2_M2^c);
    term_2_M2 = product_term_2_M2 * G_Meijer_2_M2 ;

    P_outage_M2(i) = term_1_M2 + term_2_M2;
    

    %%%%%%%%%%%%%   SImulation COde %%%%%%%%%%%%%%%%

    % --- Channel gain simulation ---
    % Turbulence model: mixture distribution (Exponential + Generalized Gamma)
    h_p_M1 = pointing_error(A_0, zeta, N_sim);
    h_p_M2 = pointing_error_M2(A_0_2, zeta, N_sim);
    h_p_M2N = pointing_error_M2N(A_0_2N, zeta, N_sim);
    h_t = eggrand(N_sim, w, lamb, a, b, c);
    h_p_IRS = pointing_error_IRS(A_0_IRS, zeta, N_sim);

    % Total channel = path loss * pointing error * turbulence
    h_M1 = h_l_M1 * h_p_M1 .* h_t;
    h_IRS = h_l_IRS * h_p_IRS .* h_t;
    h_M2 = h_l_M2 * h_p_M2 .* h_t;
    h_M2N = h_l_M2N * h_p_M2N .* h_t;
    
    % Instantaneous SNR for M1 user
    gamma_X1_L1_IRS = (eta * a_1 * rho * N^2 * h_IRS) ./ ((eta * epsilon * a_2 * rho * N^2 * h_IRS)+1) ;
    gamma_X2_L1_IRS = (eta * a_2 * rho * N^2 * h_IRS) ./ ((eta * epsilon * a_1 * rho * N^2 * h_IRS)+1) ;

    gamma_X1_L1 = (eta * a_1 * rho * h_M1) ./ ((eta * epsilon * a_2 * rho  * h_M1)+1) ;
    gamma_X2_L1 = (eta * a_2 * rho * h_M1) ./ ((eta * epsilon * a_1 * rho  * h_M1)+1) ;
    
    gamma_X2_L2N = eta * a_2 * rho * h_M2N ;

    gamma_X2_L2 = (eta * a_2 * rho * h_M2) ./ ((eta  * epsilon * a_1 * rho  * h_M2)+1) ;
   
    % Outage check: log2(1+SNR) < R_s
    outage_IRS = (gamma_X2_L1_IRS < gamma_th_2) & (gamma_X1_L1_IRS < gamma_th_1);
    P_out_MC_IRS(i) = mean(outage_IRS);

    outage = (gamma_X2_L1 < gamma_th_2) & (gamma_X1_L1 < gamma_th_1);
    P_out_MC(i) = mean(outage);

    outage_M2 = gamma_X2_L2 < gamma_th_2;
    P_out_MC_M2(i) = mean(outage_M2);

    outage_M2N = gamma_X2_L2N < gamma_th_2;
    P_out_MC_M2N(i) = mean(outage_M2N);
    
end
% **Plot the results**
figure;
semilogy(gamma_dB, P_outage_IRS, 'b--', 'LineWidth', 1.5);
hold on;
semilogy(gamma_dB, P_out_MC_IRS, 'rs', 'LineWidth', 1.5);
semilogy(gamma_dB, P_outage_M1, 'b--', 'LineWidth', 1.5);
semilogy(gamma_dB, P_out_MC, 'r^', 'LineWidth', 1.5);
semilogy(gamma_dB, P_outage_M2, 'b-', 'LineWidth', 1.5);
semilogy(gamma_dB, P_out_MC_M2, 'r*', 'LineWidth', 1.5);
semilogy(gamma_dB, P_outage_M2N, 'b-', 'LineWidth', 1.5);
semilogy(gamma_dB, P_out_MC_M2N, 'rp', 'LineWidth', 1.5);
xlabel('SNR (dB)');
ylabel('Outage Probability');
legend('Analytical of L_1 user with IRS','Simulation of L_1 user with IRS', ...
   'Analytical of L_1 user without IRS','Simulation of L_1 user without IRS', ...
   'Analytical of L_2 user without NOMA','Simulation of L_2 user without NOMA', ...
   'Analytical of L_2 user with NOMA','Simulation of L_2 user with NOMA','Location','south west');
grid on;
xlim([20 60]);
%ylim([1e6 1]);

function h_t = eggrand(N_sim, w, lamb, a, b, c)
    U = rand(N_sim,1);
    h_t = zeros(N_sim,1);

    % Exponential branch
    idx1 = (U <= w);
    h_t(idx1) = exprnd(lamb, sum(idx1),1);

    % Generalized Gamma branch
    idx2 = ~idx1;
    Y = gamrnd(a,1,sum(idx2),1);   % Gamma(a,1)
    h_t(idx2) = b * (Y.^(1/c));  % GG transformation
end

% Pointing error random generation
function h_p_M1 = pointing_error(A_0, zeta, N_sim)
    U = rand(N_sim,1);        % Uniform(0,1)
    R = sqrt(-2*log(U));
    h_p_M1 = A_0 * exp(- (R.^2) ./ (2*zeta^2));
end

% Pointing error random generation
function h_p_IRS = pointing_error_IRS(A_0_IRS, zeta, N_sim)
    U = rand(N_sim,1);        % Uniform(0,1)
    R = sqrt(-2*log(U));
    h_p_IRS = A_0_IRS * exp(- (R.^2) ./ (2*zeta^2));
end

% Pointing error random generation
function h_p_M2 = pointing_error_M2(A_0_2, zeta, N_sim)
    U = rand(N_sim,1);        % Uniform(0,1)
    h_p_M2 = A_0_2 * (U.^(1/zeta^2));
end

% Pointing error random generation
function h_p_M2N = pointing_error_M2N(A_0_2N, zeta, N_sim)
    U = rand(N_sim,1);        % Uniform(0,1)
    h_p_M2N = A_0_2N * (U.^(1/zeta^2));
end
