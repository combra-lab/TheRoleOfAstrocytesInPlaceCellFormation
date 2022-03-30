%% Experiment with a population of CA1 cells
% A traversal through a location results in the emergence of place fields 
% in a population of CA1 cells under the astrocytic modulation. This file 
% generates the spatial encodings of the track shown in Fig. 3, G and 
% Fig. 6, A-B

clear all
clc
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Track specification 
scale = 1; % parameter to rescale the experiment in space
track_l = 1.87/scale; % track length in meters
lap_time = round(4100/(scale)); % lap time in ms
v = track_l/lap_time; % speed in m/ms
track_shape = 1; % Track shape: 0 for linear and 1 for circular/treadmill
number_laps = 1;
rest_laps = []; % Specification of resting laps (without sensory stimulation)

%% %%%%%%%%%%%%% CA3 PC firing rate specification as a function of position 
N_CA3 = 200;
N_CA1 = 400;

T = lap_time; % Simulation time set to lap time
dt = 0.001; % Simulation timestep
dt_r = 100*dt; % Sampling timestep
time = dt:dt:T;
time_r = dt:dt_r:T;
tau_dec = 30; % Current deacy time constant
x = v*time; % position of the mouse

CA3_PC_peak_dist = track_l/N_CA3;
y_CA3_PC = CA3_PC_peak_dist:CA3_PC_peak_dist:track_l;

d = zeros(N_CA3, length(x));
if track_shape == 0
    for ii = 1:N_CA3 % Distance for a Linear track
       d(ii, :) = -abs((x - y_CA3_PC(ii)));  
    end
else
    for ii = 1:N_CA3 % Distance for a Circular track
        ind1 = find(x >= y_CA3_PC(ii) + track_l/2); 
        d(ii, ind1) = (-track_l + x(ind1) - y_CA3_PC(ii));
        ind2 = find((x > y_CA3_PC(ii) - track_l/2) & (x < y_CA3_PC(ii) + track_l/2));
        d(ii, ind2) = (x(ind2) - y_CA3_PC(ii));  
        ind3 = find(x <= y_CA3_PC(ii) - track_l/2);
        d(ii, ind3) = (x(ind3) - y_CA3_PC(ii)) + track_l;
    end
end


sigma = 0.21/scale;
R = exp(-(d/sigma).^2); 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Create Neuron to Neuron Connectivity
g_bar = 5.25/N_CA3; 
g_bar_feed = 0.04/N_CA3;

w_max = 25; % Max-Min weights
w_min = 0;
w_CA3_CA1 = normrnd(12, 6, [N_CA3, N_CA1]);
conn_prob = 1;

for ii = 1:N_CA3
    ind = randperm(N_CA1, round((1-conn_prob)*N_CA1));
    w_CA3_CA1(ii, ind) = 0;
end

w_CA3_CA1(w_CA3_CA1 > w_max) = w_max;
w_CA3_CA1(w_CA3_CA1 < w_min) = 0.1;
init_tot_wgt = sum(w_CA3_CA1);

init_w_CA1_CA1_feed_inh = -1000*(ones(N_CA1, N_CA1) - eye(N_CA1, N_CA1));
w_CA1_CA1_feed_inh = init_w_CA1_CA1_feed_inh;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Create Neuron to Astro Connectivity
g_bar_A = 0.05/300000;

w_N_A = zeros(N_CA1, N_CA1/4);
for j=1:1:N_CA1/4
    w_N_A((j-1)*4 + 1:(j-1)*4 + 4, j) = 1;
end
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LAP ITERATION
exp_folder = sprintf(datestr(datetime('now'))); % Create folder to save data
mkdir(exp_folder)
tic
for l = 1:1:number_laps
    l
    %% %%%%%%%%%%%%%%%%%%%%%%%% CA1 PC Input specification as CA3 PC firing 
    V_ex_i = zeros(N_CA3, length(time));
    I_ex_i = zeros(N_CA3, length(time));
    
    if ismember(l, rest_laps) == 0
        fr = 40;
        for i=1:1:(N_CA3)
            V_ex_i(i, :) = PoissonTrain_mod(dt, T, fr, R(i, :));
        end

        for i=2:1:size(V_ex_i, 2)
                I_ex_i(:, i) = I_ex_i(:, i-1) - dt*I_ex_i(:, i-1)/tau_dec;
            for j=1:1:N_CA3    
                if V_ex_i(j, i) == 1
                    I_ex_i(j, i) = I_ex_i(j, i) + 0.9;  
                end
            end
        end
        I_stim = I_ex_i;
        clear I_ex_i;
    else
        I_stim = I_ex_i;
    end

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialize CA1 PCs
    if l == 1 % Initialize before the first lap
    	V_CA1 = zeros(N_CA1, 3);
        V_CA1(:, 3) = -70;
        AP_CA1 = sparse(N_CA1, T/dt);
        I_ex_CA1 = zeros(N_CA3, T/dt);
        I_feed_inh = zeros(N_CA1, 1);

        V_CA1_r = zeros(N_CA1, T/dt_r);
        V_CA1_r(:, 1) = V_CA1(:, 3);

        % Tr_pre = zeros(N, 1);
        Tr_CA1 = zeros(N_CA1, 1);
        Tr_CA3 = zeros(N_CA3, 1);

        h_Na_CA1 = h_inf_Na(V_CA1(:, 3));
        n_K_CA1 = n_inf_K(V_CA1(:, 3));
        b_A_CA1 = b_inf_A(V_CA1(:, 3));
        z_M_CA1 = z_inf_M(V_CA1(:, 3));
        
        modes_CA1 = 1*ones(N_CA1, 1);
        Ca_CA1 = 1.2*ones(N_CA1, 1);
    else % Or initialize with the last values of the previous lap after the first lap
        temp = AP_CA1(:, end);
        AP_CA1 = sparse(N_CA1, T/dt);
        AP_CA1(:, 1) = temp;
        
        temp = V_CA1_r(:, end);
        V_CA1_r = zeros(N_CA1, T/dt_r);
        V_CA1_r(:, 1) = temp;
    end    

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialize Astros
    if l == 1
        V_A = -80*ones(N_CA1/4, 1);
        Ca = 1.2*ones(N_CA1/4, 1);
        I_A = zeros(N_CA1/4, 1);
        V_A_r = zeros(N_CA1/4, T/dt_r);
        Ca_r = zeros(N_CA1/4, T/dt_r);
        I_A_r = zeros(N_CA1/4, T/dt_r);
    else
        temp = V_A_r(:, end);
        V_A_r = zeros(N_CA1/4, T/dt_r);
        V_A_r(:,1) = temp;
        
        temp = Ca_r(:, end);
        Ca_r = zeros(N_CA1/4, T/dt_r);
        Ca_r(:,1) = temp;
        
        temp = I_A_r(:, end);
        I_A_r = zeros(N_CA1/4, T/dt_r);
        I_A_r(:,1) = temp;
    end
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Run Network
    idx = 0;
    for ii = 1:1:(T/dt)%-1
            if mod(ii, 10*(T/dt)/100) == 0
                percent = (ii/(T/dt))*100 % Print experiment evolution (%)
            end
        I_ex_CA1 = I_stim(:, ii);
        I_ex_CA1_weighted = g_bar*sum(bsxfun(@times, I_ex_CA1, w_CA3_CA1/w_max))';

        I_A = I_A - dt*I_A/tau_dec;
        I_A = I_A + g_bar_A*sum(bsxfun(@times, I_ex_CA1_weighted, w_N_A))';

        I_feed_inh_weighted = g_bar*sum(bsxfun(@times, I_feed_inh, w_CA1_CA1_feed_inh))';
        I_ex_CA1_weighted = I_ex_CA1_weighted + I_feed_inh_weighted;

        V_CA1(:, 1) = V_CA1(:, 2);
        V_CA1(:, 2) = V_CA1(:, 3);
        Ca_CA1 = sum(bsxfun(@times, Ca, w_N_A'));


        [V_CA1(:, 3), h_Na_CA1, n_K_CA1, b_A_CA1, z_M_CA1] = HH(dt, V_CA1(:, 2), I_ex_CA1_weighted, h_Na_CA1, n_K_CA1, b_A_CA1, z_M_CA1, modes_CA1, Ca_CA1);
        [Ca, V_A] = astro(dt, Ca, V_A, I_A);
        AP_CA1(:, ii) = spike_detect(V_CA1(:, 1), V_CA1(:, 2), V_CA1(:, 3));

        spiking = find(AP_CA1(:, ii)==1);

        I_feed_inh = I_feed_inh - dt*I_feed_inh/tau_dec;
        if numel(spiking) > 0
            I_feed_inh(spiking) = I_feed_inh(spiking) + g_bar_feed*V_CA1(spiking, 2);
        end
        
        [w_CA3_CA1, Tr_CA3, Tr_CA1] = STDP(dt, w_CA3_CA1, V_ex_i(:, ii), AP_CA1(:, ii), Tr_CA3, Tr_CA1);
        new_tot_wgt = sum(w_CA3_CA1);
        mult = (1 - ((abs(init_tot_wgt-new_tot_wgt))./init_tot_wgt)).^2;
        w_CA1_CA1_feed_inh = bsxfun(@times, mult', init_w_CA1_CA1_feed_inh');
    end
    
    cd(exp_folder);
    lap_folder = sprintf('Lap%d', l); 
    mkdir(lap_folder)
    cd(lap_folder);
    save('AP_CA1.mat', 'AP_CA1');
    cd('../..')
    
end
toc
%% Plot CA1 spiking activity 
figure();
[i, j] = find(AP_CA1 == 1);
scatter(j*dt, i, 15, 'k', '.')
axis([-5 T+5 0.5 (N_CA1 + 0.5)])
%% Plot CA1 spiking activity sorted over time
fAP = zeros(N_CA1, 1);
for kk = 1:N_CA1
    if numel(find(AP_CA1(kk, :) == 1)) > 0
        fAP(kk) = min(find(AP_CA1(kk, :) == 1))*dt;
    end
end
fAP(:, 2) = 1:N_CA1; 
fAP_sorted = sortrows(fAP);
AP_CA1_sorted = AP_CA1(fAP_sorted(:, 2), :);

figure();
[i, j] = find(AP_CA1_sorted == 1);
scatter(j*dt, i, 15, 'k', '.')
axis([-5 T+5 0.5 (N_CA1 + 0.5)])

%% Bin spiking activity of each CA1 neuron and plot it sorted w.r.t. its peak
total_AP = AP_CA1;
bin_size = 5; % bin size in cm
bins = floor(track_l*(100/bin_size));
activity = zeros(N_CA1, bins);
time_bin = floor(length(total_AP)/bins);
for kk = 1:bins
    activity(:, kk) = sum(total_AP(:, round((kk-1)*time_bin+1):min(length(total_AP), round(kk*time_bin))), 2);
end

bin_of_max_act = zeros(N_CA1, 1);
for ll = 1:N_CA1
    if numel(find(activity(ll, :) > 0)) > 0
        bin_of_max_act(ll) = min(find(activity(ll, :) == max(activity(ll, :))));
    end
end
bin_of_max_act(:, 2) = 1:N_CA1; 
bin_of_max_act_sorted = sortrows(bin_of_max_act);
activity_sorted = activity(bin_of_max_act_sorted(:, 2), :);

figure();
imagesc(activity_sorted)
