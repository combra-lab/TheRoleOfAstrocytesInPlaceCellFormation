%% Single CA1 cell experiment
% A long familiarization traversal through a location results in the
% emergence of a place field in a CA1 cell under the astrocytic modulation.
% Multiple revisits to the location demonstrate the transient nature of the
% place field, which fades when the astrocytic modulation is canceled by
% the return of feedforward inhibition to its high resting levels.
% This file recreates Figure 3, A-E.

clear all
clc

%Constants 
N_SN = 100; % Number of CA3 neurons
N_PM = 2; % Number of CA1 neurons

dt = 0.001; % Integration Time Step ms
dt_r = 100*dt; % Record Time Step ms
T = 15000; % Simulation time

tau_dec = 30; % Voltage decay time constant
g_bar = 0.225/(0.8*N_SN); % Synaptic conductance

time = [dt:dt:T]; % Simulation time Array
time_r = [dt:dt_r:T]; % Recording time Array

%% Input Specification
I_ex_i = zeros(N_SN, T/dt);
V_ex_i = zeros(N_SN, T/dt);

% Setting baseline stimulation (off-location)
base_fr = 2; % Baseline firing rate
for i=1:1:(N_SN) 
    V_ex_i(i, :) = PoissonTrain(dt, T, base_fr);
end

% Setting place-dependent sensory stimulation (on-location)
% Load the 40 neurons to recreate the experiment in the paper or set a 
% different group of neurons to be activated
active_neurons = dlmread('40_active_neurons.txt'); 
exc_fr = 40;
pulses = 7;
stim_starts = [1000, 3000, 4000, 6000, 8000, 10000, 13000];
stim_ends = [2500, 3150, 4150, 6150, 8150, 10150, 13150];

stim_start = zeros(pulses, length(active_neurons));
stim_end = zeros(pulses, length(active_neurons));
for pp = 1:pulses
    stim_start(pp, :) = stim_starts(pp)*ones(1, length(active_neurons));
    stim_end(pp, :) = stim_ends(pp)*ones(1, length(active_neurons));
end

fr = exc_fr*ones(1, length(active_neurons));

for j =1:pulses
    for i =1:length(active_neurons)
        V_ex_i(active_neurons(i), stim_start(j, i)/dt:stim_end(j, i)/dt-1) = PoissonTrain(dt, stim_end(j, i)-stim_start(j, i), fr(i));
    end
end
for i=2:1:(T/dt)
    I_ex_i(:, i) = I_ex_i(:, i-1) - dt*I_ex_i(:, i-1)/tau_dec;
    for j=1:1:N_SN    
        if V_ex_i(j, i) > 0
            I_ex_i(j, i) = I_ex_i(j, i) + 0.9;
        end
    end
end
clear V_ex_i;
%% Initialize astrocytic dynamics
Ca = 1.2*ones(1, T/dt); % Extracellular calcium conentration

% Using the values of the astrocytic membrane voltage that are given by the
% biologically-realistic model presented in "The neuroglial potassium cycle 
% during neurotransmission: role of kir4. 1 channels", Sibille et al., PLoS 
% Computational Biology, 2015. The code for the model was provided after 
% personal communication with the authors.
V_A = dlmread('V_A.txt'); % Astrocytic voltage
%% Initialize CA3 neurons

V_SN = zeros(N_SN, 3);
V_SN_r = zeros(N_SN, T/dt_r);
V_SN(:, 3) = -70; 
V_SN_r(:, 1) = V_SN(:, 3);
AP_SN = sparse(N_SN, T/dt);
I_out_SN = zeros(N_SN,1);
Tr_pre_SN = zeros(N_SN, 1);

h_Na_SN = h_inf_Na(V_SN(:, 3));
n_K_SN = n_inf_K(V_SN(:, 3));
b_A_SN = b_inf_A(V_SN(:, 3));
z_M_SN = z_inf_M(V_SN(:, 3));

modes_SN = zeros(N_SN, 1);
%% Initialize CA1 neurons
I_stim_PM = zeros(N_PM, 1);
V_PM = zeros(N_PM, 3);
V_PM_r = zeros(N_PM, T/dt_r);
V_PM(:, 3) = -70; 
V_PM_r(:, 1) = V_PM(:, 3);
AP_PM = sparse(N_PM, T/dt);
Tr_post_PM = zeros(N_PM, 1);

h_Na_PM = h_inf_Na(V_PM(:, 3));
n_K_PM = n_inf_K(V_PM(:, 3));
b_A_PM = b_inf_A(V_PM(:, 3));
z_M_PM = z_inf_M(V_PM(:, 3));

modes_PM = [0 1];
%% Initialize CA3-CA1 connectivity
w_r = zeros(N_SN, N_PM, T/dt_r);
w(:, 1, 1) = normrnd(8, 1, [N_SN, 1]);
w(:, 2, 1) = w(:, 1, 1);

w_max = 30;
w_min = 0;
%% Run Simulation
tic
idx = 1;
for ii = 1:1:(T/dt)-1 % Update network state in each timestep
    if mod(ii, (T/dt)/100) == 0
        percent = (ii/(T/dt))*100 % Print simulation progress (%) in the command line
    end
    % Replace neuron voltage values to create space for the new timestep
    V_SN(:, 1) = V_SN(:, 2);
    V_SN(:, 2) = V_SN(:, 3);

    V_PM(:, 1) = V_PM(:, 2);
    V_PM(:, 2) = V_PM(:, 3);
    
    % Calculate the new values of the neuron voltage
    [V_SN(:, 3), h_Na_SN, n_K_SN, b_A_SN, z_M_SN] = HH(dt, V_SN(:, 2), I_ex_i(:, ii), h_Na_SN, n_K_SN, b_A_SN, z_M_SN, modes_SN, Ca(ii));
    [V_PM(:, 3), h_Na_PM, n_K_PM, b_A_PM, z_M_PM] = HH(dt, V_PM(:, 2), I_stim_PM, h_Na_PM, n_K_PM, b_A_PM, z_M_PM, modes_PM, Ca(ii));
    
    % Calculate the new values of the astrocyte-modulated extracellular Ca
    [Ca(ii+1), x] = astro(dt, Ca(ii), V_A(ii), I_stim_PM(2));
    if ii>1
        % Detect spikes in the neurons
        AP_SN(:, ii) = spike_detect(V_SN(:, 1), V_SN(:, 2), V_SN(:, 3));
        AP_PM(:, ii) = spike_detect(V_PM(:, 1), V_PM(:, 2), V_PM(:, 3));

        spiking = find(AP_SN(:, ii)==1);
        non_spiking = find(AP_SN(:, ii)==0);
        
        % Update output currents of the spiking neurons
        I_out_SN = I_out_SN - dt*I_out_SN/tau_dec;
        if numel(spiking) > 0
            I_out_SN(spiking) = I_out_SN(spiking) + g_bar*V_SN(spiking, 2);
        end

        I_out_SN_weighted = bsxfun(@times, I_out_SN, w/w_max);
        I_stim_PM = sum(I_out_SN_weighted)';
    end
    
    % Update connection weights between CA3 and CA1 neurons
    [w, Tr_pre_SN, Tr_post_PM] = STDP(dt, w, AP_SN(:, ii), AP_PM(:, ii), Tr_pre_SN, Tr_post_PM); % STDP updates

    % Sample and record values every dt_r time steps to reduce memory usage
    if mod(ii*dt, dt_r) == 0 % record data
        idx = idx + 1;
        V_SN_r(:, idx) = V_SN(:, 3);
        V_PM_r(:, idx) = V_PM(:, 3);
        w_r(:, :, idx) = w;
    end        
end
toc
%% Plots
% Raster Plot 
AP = [AP_SN; AP_PM];
figure();
[i, j] = find(AP == 1);
scatter(j*dt, i, 15, 'k', '.')
axis([-5 (T+5) 0.5 (N_SN + N_PM + 0.5)])
yticks([1:(N_SN + N_PM)]);

subplots = 3;
figure()
subplot(subplots, 1, 1)
plot(time_r, V_PM_r(2,:))
xlabel('Time (ms)');
ylabel('CA1 neuron Voltage (mV)');
axis([0 T -100 100])
subplot(subplots, 1, 2)
plot(time, V_A)
xlabel('Time (ms)');
ylabel('Astrocytic Voltage (mV)');
subplot(subplots, 1, 3)
plot(time, Ca)
xlabel('Time (ms)');
ylabel('Calcium Concentration(mM)');