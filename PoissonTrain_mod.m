function [train] = PoissonTrain_mod(dt, dur, fr, mod)
%   PosissonTrain_mod : Generate a space-modulated Poisson spike train
% param(scalar) : dt (timestep)
% param(scalar) : dur (train duration)
% param(scalar) : fr (average firing rate)
% param(scalar) : mod (modulation factor)

% Generate a Poisson spike train of duration dur with a firing rate of fr, 
% whose firing rate is modulated as a function of the location.
    prob = rand(1, dur/dt);
    train = zeros(1, dur/dt);
    train(prob < fr*dt*0.001*mod) = 1;
end

