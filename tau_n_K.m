function tau = tau_n_K(V)
%   tau_n_K : Dynamics for the time constant of the n parameter of the K channel
% param : V (membrane voltage) 

%   Calculation of the time constant of the n parameter for the Na channel of the conditionally-bursting pyramidal neuron
%   The formulation and the parameters follow the model in Golomb, D., Yue, C. and Yaari, Y., 2006. Contribution of persistent Na+ current and M-type K+ current to somatic bursting in CA1 pyramidal cells: combined experimental and modeling study. Journal of neurophysiology, 96(4), pp.1912-1926.
    theta_nt = -27;
    sigma_nt = -15;
    
    tau = 0.1 + 0.5*(1./(1 + exp(-(V - theta_nt)/sigma_nt)));
end
