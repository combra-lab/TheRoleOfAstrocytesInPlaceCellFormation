function a = a_inf_A(V)
%   a_inf_A : Dynamics for the a parameter of the delayed rectifier K channel
% param : V (membrane voltage) 

%   Calculation of the a parameter for the delayed rectifier K channel of the conditionally-bursting pyramidal neuron
%   The formulation and the parameters follow the model in Golomb, D., Yue, C. and Yaari, Y., 2006. Contribution of persistent Na+ current and M-type K+ current to somatic bursting in CA1 pyramidal cells: combined experimental and modeling study. Journal of neurophysiology, 96(4), pp.1912-1926.
    theta_a = -50;
    sigma_a = 20;
    
    a = 1./(1 + exp((-(V - theta_a))/sigma_a));
end
