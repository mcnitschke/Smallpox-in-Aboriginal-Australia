function [S,E,I,R,TE] = Spread(S,E,I,R,TE,Amat,dt,Bmat)
% Adaptive tau-leaping approximation for disease spread.
% P is the propensity
% Tau leap algorithm....
% Use the adaptive tau to deep tau optimal
% TE - total unique exposed individuals

M = size(S,1);  % Number of patches
SumI = sum(I);

t = 0;

while t < dt

    P = Bmat .* S .* SumI ./ Amat;

    tau = adaptive_tau(S, P, dt - t);
    New_E = poissrnd(P .* tau);

    % Update compartments
    S = S - New_E;
    E = E + New_E;

    % Update time
    t = t + tau;
    TE = TE + sum(New_E(:));
    S(S < 0) = 0;
end

end
