function z = z_inf_M(V)
%   m_inf_Na : Dynamics for the z parameter of the muscarinic-sensitive K channel
% param : V (membrane voltage) 

%   Calculation of the z parameter for the muscarinic-sensitive K channel of the conditionally-bursting pyramidal neuron
%   The formulation and the parameters follow the model in Golomb, D., Yue, C. and Yaari, Y., 2006. Contribution of persistent Na+ current and M-type K+ current to somatic bursting in CA1 pyramidal cells: combined experimental and modeling study. Journal of neurophysiology, 96(4), pp.1912-1926.
    theta_z = -39;
    sigma_z = 5;
    
    z = 1./(1 + exp((-(V - theta_z))/sigma_z));
end
