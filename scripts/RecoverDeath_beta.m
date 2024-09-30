function[S,E,I,R,Death]=...
    RecoverDeath_beta(S,E,I,R,Death,alpha,delta,gamma)
%I->R & R->D
%Created 3 May, 2024
Inew_E = binornd(E,alpha);   % E--> I1
I = I + Inew_E;
E = E - Inew_E;

% Deathnew_I = binornd(I,delta); % I3--> D
% I = I - Deathnew_I;
% Death = Death + Deathnew_I;

End_I = binornd(I,gamma); % I4-->R
% Deathnew = binornd(End_I,delta);
% Rnew_I = End_I - Deathnew;
Death = Death + End_I;
R = R + End_I;
I = I - End_I;
end