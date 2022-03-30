function [weights, TrPre, TrPost] = STDP(dt, weights, firing_pre, firing_post, TrPre, TrPost)
%   STDP : Spike-timing-dependent plasticity (STDP) learning dynamics 
% param(scalar) : dt (timestep)
% param(array) : weights (connection weights modified by the STDP dynamics)
% param(vector) : firing_pre (binary vector with 1s indicating spikes of the respective presynaptic neuron at the current timestep)
% param(vector) : firing_post (binary vector with 1s indicating spikes of the respective postsynaptic neuron at the current timestep)
% param(vector) : Tr_pre (spike trace vector of the presynaptic neurons)
% param(vector) : Tr_pre (spike trace vector of the postsynaptic neurons)

% Updating the synaptic weights of a set of connections using the STDP
% learning rule
    a = 0.175;
    b = 0.5*a;
    
    tau_dec_pre_tr = 60;
    tau_dec_post_tr = 60;
    
    w_max = 25;
    w_min = 0;
    
    TrPre = TrPre - dt*TrPre/tau_dec_pre_tr;
    TrPost = TrPost - dt*TrPost/tau_dec_post_tr;
    
    TrPre(firing_pre == 1) = 1;
    TrPost(firing_post == 1) = 1;
    
    
    if sum(firing_pre) > 0
        weights(find(firing_pre == 1), :) = weights(find(firing_pre == 1), :) - b*TrPost';
    end
    if sum(firing_post) > 0
        weights(:, find(firing_post == 1)) = weights(:, find(firing_post == 1)) + a*TrPre;
    end
    
    weights(weights > w_max) = w_max;
    weights(weights < w_min) = w_min;  
end

