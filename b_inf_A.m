function b = b_inf_A(V)
%   b_inf_A : Dynamics for the b parameter of the A-type K channel
% param : V (membrane voltage) 

%   Calculation of the b parameter for the A-type K channel of the conditionally-bursting pyramidal neuron
%   The formulation and the parameters follow the model in Golomb, D., Yue, C. and Yaari, Y., 2006. Contribution of persistent Na+ current and M-type K+ current to somatic bursting in CA1 pyramidal cells: combined experimental and modeling study. Journal of neurophysiology, 96(4), pp.1912-1926.
    theta_b = -80;
    sigma_b = -6;
    
    b = 1./(1 + exp((-(V - theta_b))/sigma_b));
end
