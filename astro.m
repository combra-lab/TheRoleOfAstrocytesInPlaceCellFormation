function [Ca_c, V_A_c] = astro(dt, Ca_p, V_A_p, I_ex)
%   astro : Astrocytic dynamics 
% param(scalar) : dt (timestep)
% param(vector) : Ca_p (value of the extracellular calcium concentration in the vicinity of an astrocyte in the previous timestep)
% param(vector) : V_A_p (value of the astrocytic membrane voltage in the previous timestep)

% Astrocytic dynamics for the membrane potentail and the extracellular
% calcium in the vicinity

    V_rest = -80;
    tau_Ca_dec = 2000;
    tau_Ca_rise = 1000;
    tau_volt_dec = 1000;    
    Ca_base = 1.2;

    V_A_c = V_A_p + dt*I_ex - dt*(V_A_p - V_rest)/tau_volt_dec;
    
    Ca_c = Ca_p - dt*Ca_p.*(V_A_p - V_rest)/tau_Ca_dec + dt*(Ca_base - Ca_p)/tau_Ca_rise;

end

