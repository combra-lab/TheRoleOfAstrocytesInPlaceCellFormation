function shift = astro_shift(Ca_e)
%   astro_shift : Astrocyte-induced shift of the activation curve of the 
%                 persistent Na channel by extracellular calcium concentration 
% param : Ca_e (extracellular calcium concentration ) 

% Shifting of the activation threshold of the activation curve of the 
% persistent Na channel by extracellular calcium concentration. 

    shift = zeros(length(Ca_e), 1);
    shift(Ca_e <= 0.6) = 6;
    ind1 = Ca_e > 0.6; 
    ind2 = Ca_e < 1.2;
    ind = and(ind1, ind2);
    shift(ind) = -20*(Ca_e(ind) - 1);

end
