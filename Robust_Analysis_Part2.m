%% Robust Control SC42145 - Part 2: Robust Analysis and Controller Design
% Group 49
% Date: December 31, 2025

clear; close all; clc;

%% 1. Load Data
% Load assignment data
if exist('Assignment_Data_SC42145_2025.mat', 'file')
    load('Assignment_Data_SC42145_2025.mat');
else
    warning('Assignment_Data_SC42145_2025.mat not found.');
end

% Import transfer functions
if exist('FWT', 'var')
    beta_to_omega_r = FWT(1,1);
    beta_to_z = FWT(2,1);
    tau_to_omega_r = FWT(1,2);
    tau_to_z = FWT(2,2);
    % Construct plant
    G = [tf(beta_to_omega_r), tf(tau_to_omega_r);
         tf(beta_to_z), tf(tau_to_z)];
    G_nom = G; % Use G as G_nom for compatibility
else
    warning('FWT variable not found in loaded data.');
end

%% 2. Define Weights
% Construct performance weight Wp
f_cutoff = 0.3;
w_c = 2*pi*f_cutoff;
attenuation = 1e-4;
H_inf_norm = 2.4;
Wp11 = tf([1/H_inf_norm, w_c], [1, w_c*attenuation]);
Wp = [Wp11, 0; 0, 0.1];

% Construct control weight Wu
Wu22 = tf([5e-3, 7e-4, 5e-5], [1, 14e-4, 1e-6]);
Wu = [0.01, 0; 0, Wu22];

% construct weight matrix for input and output deviation
s = tf('s');
W_i1 = ((1/(16*pi))*s + 0.2)/((1/(64*pi))*s + 1);
W_i2 = ((1/(16*pi))*s + 0.2)/((1/(64*pi))*s + 1);
W_o1 = (0.05*s + 0.25)/(0.01*s + 1);
W_o2 = (0.05*s + 0.25)/(0.01*s + 1);

W_i = [W_i1, 0; 0, W_i2];
W_o = [W_o1, 0; 0, W_o2];

disp('Weights defined successfully.');

%% 4. Define Generalized Plant for Robust Analysis (Part 2.1)
% Construct generalized plant P_gen using user-specified structure

% Plant
G_d = G;
G_d.u = {'u_in(1)', 'u_in(2)'};
G_d.y = {'y(1)', 'y(2)'};

% Performance weight
Wp_blk = Wp; 
Wp_blk.u = {'e(1)', 'e(2)'};
Wp_blk.y = {'z1(1)', 'z1(2)'};

% Control weight
Wu_blk = Wu;
Wu_blk.u = {'u(1)', 'u(2)'};
Wu_blk.y = {'z2(1)', 'z2(2)'};

% Input uncertainty weight
W_in = W_i;
W_in.u = {'u(1)', 'u(2)'};
W_in.y = {'y_di(1)', 'y_di(2)'};

% Output uncertainty weight
W_out = W_o;
W_out.u = {'y(1)', 'y(2)'};
W_out.y = {'y_do(1)', 'y_do(2)'};

% Summation blocks
Sum_u1 = sumblk('u_in(1) = u(1) + u_di(1)');
Sum_u2 = sumblk('u_in(2) = u(2) + u_di(2)');
Sum_o1 = sumblk('y_o(1) = y(1) + u_do(1)');
Sum_o2 = sumblk('y_o(2) = y(2) + u_do(2)');
Sum_e1 = sumblk('e(1) = y_o(1) - w(1)');
Sum_e2 = sumblk('e(2) = y_o(2) - w(2)');

% Build interconnected system P_gen (Open Loop for Uncertainty)
P_gen = connect(G_d, W_out, W_in, Wp_blk, Wu_blk, ...
    Sum_u1, Sum_o1, Sum_e1, Sum_u2, Sum_o2, Sum_e2, ...
    {'w(1)', 'w(2)', 'u(1)', 'u(2)', 'u_di(1)', 'u_di(2)', 'u_do(1)', 'u_do(2)'}, ...
    {'z1(1)', 'z1(2)', 'z2(1)', 'z2(2)', 'y_di(1)', 'y_di(2)', 'y_do(1)', 'y_do(2)', 'e(1)', 'e(2)'});

% Define Uncertainty Blocks for Synthesis/Analysis
delta_i1 = ultidyn('delta_i1', [1 1]);
delta_i2 = ultidyn('delta_i2', [1 1]);
delta_o1 = ultidyn('delta_o1', [1 1]);
delta_o2 = ultidyn('delta_o2', [1 1]);

% Create Uncertain System P_rob by closing the loop with Delta
Delta_i_blk = blkdiag(delta_i1, delta_i2);
Delta_i_blk.u = {'y_di(1)', 'y_di(2)'};
Delta_i_blk.y = {'u_di(1)', 'u_di(2)'};

Delta_o_blk = blkdiag(delta_o1, delta_o2);
Delta_o_blk.u = {'y_do(1)', 'y_do(2)'};
Delta_o_blk.y = {'u_do(1)', 'u_do(2)'};

P_rob = connect(G_d, W_out, W_in, Wp_blk, Wu_blk, ...
    Delta_i_blk, Delta_o_blk, ...
    Sum_u1, Sum_o1, Sum_e1, Sum_u2, Sum_o2, Sum_e2, ...
    {'w(1)', 'w(2)', 'u(1)', 'u(2)'}, ...
    {'z1(1)', 'z1(2)', 'z2(1)', 'z2(2)', 'e(1)', 'e(2)'});

%% 3. Synthesize Nominal Mixed-Sensitivity Controller (Part 1)
% Re-synthesize K_nom using the new P_gen structure (ignoring uncertainty)
% P_nom_equiv inputs: w(1,2), u(1,2)
% P_nom_equiv outputs: z1(1,2), z2(1,2), e(1,2)
% This corresponds to P_rob with Delta=0.
P_nom = P_rob.NominalValue; 

nmeas = 2;
ncont = 2;

disp('Synthesizing Nominal H-infinity Controller...');
[K_nom, CL_nom, gamma_nom] = hinfsyn(P_nom, nmeas, ncont);
disp(['Nominal Gamma: ', num2str(gamma_nom)]);

% Check stability
if isstable(CL_nom)
    disp('Nominal closed-loop system is stable.');
else
    warning('Nominal closed-loop system is UNSTABLE.');
end

%% 5. Robustness Analysis (NS, NP, RS, RP)
% Form Closed Loop Uncertain System
CL_rob = lft(P_rob, K_nom);

% Analyze using robuststab and robustperf (or mussv)
disp('Performing Robustness Analysis...');

% Nominal Stability (NS)
% Already checked, but can check if CL_rob.NominalValue is stable
if isstable(CL_rob.NominalValue)
    disp('NS: OK');
else
    disp('NS: FAILED');
end

% Nominal Performance (NP)
% H-infinity norm of nominal closed loop
np_norm = norm(CL_rob.NominalValue, inf);
disp(['Nominal Performance (Gamma): ', num2str(np_norm)]);
if np_norm < 1
    disp('NP: OK (<1)');
else
    disp('NP: Not met (>1)');
end

% Robust Stability (RS)
% Calculate structured singular value for stability
% robuststab computes the margins
[stabmarg, destabunc, report, info] = robuststab(CL_rob);
disp('Robust Stability Report:');
disp(report);

% Robust Performance (RP)
% robustperf computes the margins
[perfmarg, perfunc, report_perf, info_perf] = robustperf(CL_rob);
disp('Robust Performance Report:');
disp(report_perf);

% Plotting Singular Values of Uncertain Plant (Task 2.1.3)
figure;
% sigma(G_unc); % G_unc not explicitly defined in new structure
if exist('G_unc', 'var')
    sigma(G_unc);
    title('Singular Values of Uncertain Plant G_{unc}');
else
    sigma(P_rob.NominalValue);
    title('Singular Values of Generalized Plant P_{gen}');
end
grid on;

% Plotting Uncertainty Weights (Task 2.1.3)
figure;
bodemag(W_i1, 'b', W_i2, 'r', W_o1, 'g', W_o2, 'k');
legend('W_{i1}', 'W_{i2}', 'W_{o1}', 'W_{o2}');
title('Uncertainty Weights Frequency Response');
grid on;

% Plotting Mu for NP, RS, RP (Task 2.1.5)
% We can use mussv explicitly to get the mu plots
% Extract the M matrix for the interconnection
% M = lft(P_rob, K_nom) is the uncertain closed loop.
% To plot mu, we need the frequency response.

% Plotting Mu for NP, RS, RP (Task 2.1.5)
% Use robustperf info to get bounds
if exist('info_perf', 'var') && isfield(info_perf, 'MussvBnds')
    mubounds = info_perf.MussvBnds;
    % Plot Mu
    figure;
    semilogx(mubounds.Frequency, mubounds(1,1).ResponseData(:), 'b-');
    hold on;
    semilogx(mubounds.Frequency, mubounds(1,2).ResponseData(:), 'r--');
    title('Structured Singular Value \mu (Robust Performance)');
    xlabel('Frequency (rad/s)');
    ylabel('\mu');
    legend('Lower Bound', 'Upper Bound');
    grid on;
else
    disp('Could not extract Mu bounds from robustperf info.');
end

disp('Analysis Complete.');

%% 6. Robust Controller Design (D-K Iteration) - Section 2.2

disp('-------------------------------------------------');
disp('Starting Robust Controller Design (D-K Iteration)');

% Controller Synthesis

% Options for musyn
nmeas = 2;
ncont = 2;
opts = musynOptions('Display', 'short', 'MaxIter', 10);

disp('Running D-K Iteration (musyn)...');
[K_rob, CL_rob_opt, bnd] = musyn(P_rob, nmeas, ncont, opts);

disp('D-K Iteration Complete.');
% disp(['Achieved Robust Performance (Mu): ', num2str(bnd)]);

% Form closed loop with robust controller
CL_rob_final = lft(P_rob, K_rob);

% Nominal Stability (NS)
if isstable(CL_rob_final.NominalValue)
    disp('NS (Robust Controller): OK');
else
    disp('NS (Robust Controller): FAILED');
end

% Nominal Performance (NP)
np_norm_rob = norm(CL_rob_final.NominalValue, inf);
disp(['Nominal Performance (Gamma) with Robust Controller: ', num2str(np_norm_rob)]);

% Robust Stability (RS)
[stabmarg_rob, destabunc_rob, report_rob, info_rob] = robuststab(CL_rob_final);
disp('Robust Stability:');
disp(report_rob);

% Robust Performance (RP)
[perfmarg_rob, perfunc_rob, report_perf_rob, info_perf_rob] = robustperf(CL_rob_final);
disp('Robust Performance:');
disp(report_perf_rob);

% Compare Mu Plots (Nominal vs Robust Controller)
if exist('info_perf', 'var') && isfield(info_perf, 'MussvBnds') && ...
   exist('info_perf_rob', 'var') && isfield(info_perf_rob, 'MussvBnds')

    mubounds_nom = info_perf.MussvBnds;
    mubounds_rob = info_perf_rob.MussvBnds;

    figure;
    semilogx(mubounds_nom.Frequency, mubounds_nom(1,2).ResponseData(:), 'r--', 'LineWidth', 1.5);
    hold on;
    semilogx(mubounds_rob.Frequency, mubounds_rob(1,2).ResponseData(:), 'b-', 'LineWidth', 1.5);
    title('Comparison of Robust Performance (\mu)');
    xlabel('Frequency (rad/s)');
    ylabel('\mu (Upper Bound)');
    legend('Mixed-Sensitivity (K_{nom})', 'Robust Control (K_{rob})');
    grid on;
else
    disp('Could not extract Mu bounds for comparison.');
end

%% 7. Time-Domain Simulations (Nominal vs Uncertain)

% Define Simulation Scenarios
t = 0:0.1:100;
% Step input on reference w (acting as disturbance proxy)
% w has 2 channels. Let's step w(1) (omega_r reference).
u_sim = zeros(length(t), 2);
u_sim(:,1) = 1; % Step on w(1)

% Label Controller I/O for connect
K_nom.u = {'e(1)', 'e(2)'};
K_nom.y = {'u(1)', 'u(2)'};

K_rob.u = {'e(1)', 'e(2)'};
K_rob.y = {'u(1)', 'u(2)'};

% Nominal Simulation
disp('Running Nominal Simulations...');

% Create simulation closed loops that output 'e'
% P_nom has inputs: w, u. Outputs: z1, z2, e.
% We close the loop on u, e. We want to observe e.
CL_sim_nom_Knom = connect(P_nom, K_nom, {'w(1)', 'w(2)'}, {'e(1)', 'e(2)'});
CL_sim_nom_Krob = connect(P_nom, K_rob, {'w(1)', 'w(2)'}, {'e(1)', 'e(2)'});

e_nom_Knom = lsim(CL_sim_nom_Knom, u_sim, t);
y_nom_Knom = e_nom_Knom + u_sim; 

e_nom_Krob = lsim(CL_sim_nom_Krob, u_sim, t);
y_nom_Krob = e_nom_Krob + u_sim;

% Plot Nominal Comparison
figure;
subplot(2,1,1);
plot(t, y_nom_Knom(:,1), 'r', 'LineWidth', 1); hold on;
plot(t, y_nom_Krob(:,1), 'b', 'LineWidth', 1);
title('Nominal Response to Reference Step w(1): Generator Speed \omega_r');
legend('K_{nom}', 'K_{rob}');
grid on;

subplot(2,1,2);
plot(t, y_nom_Knom(:,2), 'r', 'LineWidth', 1); hold on;
plot(t, y_nom_Krob(:,2), 'b', 'LineWidth', 1);
title('Nominal Response to Reference Step w(1): Tower Displacement z');
legend('K_{nom}', 'K_{rob}');
grid on;

% Uncertain Simulation (Worst-Case or Random Samples)
disp('Running Uncertain Simulations...');
num_samples = 10;
P_rand = usample(P_rob, num_samples);

figure;
sgtitle('Uncertain Responses (Random Samples)');
subplot(2,1,1); hold on; title('Generator Speed \omega_r');
subplot(2,1,2); hold on; title('Tower Displacement z');

for i = 1:num_samples
    P_sample = P_rand(:,:,i);
    
    % Create simulation closed loops for this sample
    CL_sim_unc_Knom = connect(P_sample, K_nom, {'w(1)', 'w(2)'}, {'e(1)', 'e(2)'});
    CL_sim_unc_Krob = connect(P_sample, K_rob, {'w(1)', 'w(2)'}, {'e(1)', 'e(2)'});
    
    e_unc_Knom = lsim(CL_sim_unc_Knom, u_sim, t);
    y_unc_Knom = e_unc_Knom + u_sim;
    
    e_unc_Krob = lsim(CL_sim_unc_Krob, u_sim, t);
    y_unc_Krob = e_unc_Krob + u_sim;
    
    subplot(2,1,1);
    plot(t, y_unc_Knom(:,1), 'r', 'LineWidth', 0.5);
    plot(t, y_unc_Krob(:,1), 'b', 'LineWidth', 0.5);
    
    subplot(2,1,2);
    plot(t, y_unc_Knom(:,2), 'r', 'LineWidth', 0.5);
    plot(t, y_unc_Krob(:,2), 'b', 'LineWidth', 0.5);
end

% Add legend to first plot only (dummy handles)
subplot(2,1,1);
h1 = plot(NaN,NaN,'r');
h2 = plot(NaN,NaN,'b');
legend([h1, h2], 'K_{nom} (Uncertain)', 'K_{rob} (Uncertain)');

disp('All tasks complete.');
