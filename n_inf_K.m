function n = n_inf_K(V)
%   n_inf_K : Dynamics for the n parameter of the K channel
% param : V (membrane voltage) 

%   Calculation of the n parameter for the K channel of the conditionally-bursting pyramidal neuron
%   The formulation and the parameters follow the model in Golomb, D., Yue, C. and Yaari, Y., 2006. Contribution of persistent Na+ current and M-type K+ current to somatic bursting in CA1 pyramidal cells: combined experimental and modeling study. Journal of neurophysiology, 96(4), pp.1912-1926.
    theta_n = -35;
    sigma_n = 10;
    
    n = 1./(1 + exp((-(V - theta_n))/sigma_n));
end
