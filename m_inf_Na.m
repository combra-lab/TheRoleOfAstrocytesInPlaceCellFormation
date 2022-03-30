function m = m_inf_Na(V)
%   m_inf_Na : Dynamics for the h parameter of the Na channel
% param : V (membrane voltage) 

%   Calculation of the m parameter for the Na channel of the conditionally-bursting pyramidal neuron
%   The formulation and the parameters follow the model in Golomb, D., Yue, C. and Yaari, Y., 2006. Contribution of persistent Na+ current and M-type K+ current to somatic bursting in CA1 pyramidal cells: combined experimental and modeling study. Journal of neurophysiology, 96(4), pp.1912-1926.   

    theta_m = -30;
    sigma_m = 9.5;
    
    m = 1./(1 + exp((-(V - theta_m))/sigma_m));
end

