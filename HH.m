function [v_c, h_Na_c, n_K_c, b_A_c, z_M_c] = HH(dt, v_p, I_p, h_Na_p, n_K_p, b_A_p, z_M_p, modes, Ca_e)
%   HH : Hodgkin-Huxley model of the conditionally-bursting pyramidal neuron
% param(scalar) : dt (timestep)
% param(vector) : v_p (value of the membrane voltage in the previous timestep)
% param(vector) : I_p (value of the input current in the previous timestep)
% param(vector) : h_Na_p (value of the h parameter of the Na current in the previous timestep)
% param(vector) : n_K_p (value of the n parameter of the K current in the previous timestep)
% param(vector) : b_A_p (value of the b parameter of the A-type K current in the previous timestep)
% param(vector) : z_M_p (value of the z parameter of the muscarinic-sensitive K current in the previous timestep)
% param(vector) : modes (type of neuron, 0 for unmodulated and 1 for conditionally-bursting neurons)
% param(vector) : Ca_e (value of the z parameter of the muscarinic-sensitive K current in the previous timestep)


%   Calculation of the membrane voltage for the conditionally-bursting pyramidal neuron
%   The formulation and the parameters follow the model in Golomb, D., Yue, C. and Yaari, Y., 2006. Contribution of persistent Na+ current and M-type K+ current to somatic bursting in CA1 pyramidal cells: combined experimental and modeling study. Journal of neurophysiology, 96(4), pp.1912-1926. 
    
    Cm = 1; % Membrane Capacitance uF/cm^2
    ENa = 55; % mV Na reversal potential
    EK = -90; % mV K reversal potential
    El = -70; % mV Leakage reversal potential

    gbarNa = 35; % mS/cm^2 Na conductance
    gbarK = 6; % mS/cm^2 K conductance
    gbarA = 1.4; % mS/cm^2 A K current conductance
    gbarM = 1; % mS/cm^2 M K current conductance
    gbarl = 0.05; % mS/cm^2 Leakage conductance
    gbarNaP = 0.3; % 0 - 0.41 mS/cm^2 NaP conductance

    theta_p = zeros(size(v_p));
    theta_p(modes == 0) = -41;
    if numel(theta_p(modes == 1)) > 0
        theta_p(modes == 1) = -41 - astro_shift(Ca_e);
    end

    
    m_Na = m_inf_Na(v_p);
    h_Na_c = h_Na_p + dt*((h_inf_Na(v_p) - h_Na_p)./tau_h_Na(v_p));
    p_NaP = p_inf_NaP(v_p, theta_p);
    n_K_c = n_K_p + dt*((n_inf_K(v_p) - n_K_p)./tau_n_K(v_p));
    a_A = a_inf_A(v_p);
    tau_b_A = 15;
    b_A_c = b_A_p + dt*((b_inf_A(v_p) - b_A_p)/tau_b_A);
    tau_z_M = 75;
    z_M_c = z_M_p + dt*((z_inf_M(v_p) - z_M_p)/tau_z_M);

    
    gNa = gbarNa*m_Na.^3.*h_Na_c;
    gNaP = gbarNaP*p_NaP;
    gK = gbarK*n_K_c.^4;
    gA = gbarA*a_A.^3.*b_A_c;
    gM = gbarM*z_M_c;
    gl = gbarl;
    
    
    INa = gNa.*(v_p - ENa);
    INaP = gNaP.*(v_p - ENa);
    IK = gK.*(v_p - EK);
    IA = gA.*(v_p - EK);
    IM = gM.*(v_p - EK);
    Il = gl*(v_p - El);
    
    
    v_c = v_p + dt*((1/Cm)*(I_p-(INa+IK+Il+INaP+IA+IM)));
end

