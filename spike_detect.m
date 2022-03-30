function spikes = spike_detect(v_p, v_c, v_n)
%   spike_detect : Detect spikes from the voltage of a HH neuron by
%   comparing three consecutive values of the neuron's voltage (v_p, v_c, v_n)
% param(scalar) : v_p (1st voltage value)
% param(scalar) : v_c (2nd voltage value)
% param(scalar) : v_n (3rd voltage value)

% Detect spikes in the continuous value of the membrane potential of a HH
% neuron. The function detects a spike when the difference between two
% consecutive pairs of voltage values switch sign.

    back_diff = (v_c - v_p) > 0; 
    forw_diff = (v_n - v_c) < 0;
    pos_volt = v_c > 0;
    
    spikes = and(and(back_diff,forw_diff), pos_volt);
end

