function p = p_inf_NaP(V, theta)
%   p_inf_NaP : Dynamics for the p parameter of the persistent Na channel
% param : V (membrane voltage)
% param : theta (activation threshold dependent on the extracellular calcium conentration) 

%   Calculation of the h parameter for the Na channel of the conditionally-bursting pyramidal neuron
%   The formulation and the parameters follow the model in Golomb, D., Yue, C. and Yaari, Y., 2006. Contribution of persistent Na+ current and M-type K+ current to somatic bursting in CA1 pyramidal cells: combined experimental and modeling study. Journal of neurophysiology, 96(4), pp.1912-1926.

    sigma_p = 3;
    
    p = 1./(1 + exp((-(V - theta))/sigma_p));
end
