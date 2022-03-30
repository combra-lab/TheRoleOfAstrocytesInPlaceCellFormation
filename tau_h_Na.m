function tau = tau_h_Na(V)
%   tau_h_Na : Dynamics for the time constant of the h parameter of the Na channel
% param : V (membrane voltage) 

%   Calculation of the time constant of the h parameter for the Na channel of the conditionally-bursting pyramidal neuron
%   The formulation and the parameters follow the model in Golomb, D., Yue, C. and Yaari, Y., 2006. Contribution of persistent Na+ current and M-type K+ current to somatic bursting in CA1 pyramidal cells: combined experimental and modeling study. Journal of neurophysiology, 96(4), pp.1912-1926.
    theta_ht = -40.5;
    sigma_ht = -6;
    
    tau = 0.1 + 0.75*(1./(1 + exp(-(V - theta_ht)/sigma_ht)));
end

