function [train] = PoissonTrain(dt, dur, fr)
%   PosissonTrain : Generate a Poisson spike train

% param(scalar) : dt (timestep)
% param(scalar) : dur (train duration)
% param(scalar) : fr (average firing rate)

% Generate a Poisson spike train of duration dur with a firing rate of fr.
    prob = rand(1, dur/dt);
    train = zeros(1, dur/dt);
    train(prob < fr*dt*0.001) = 1;
end

