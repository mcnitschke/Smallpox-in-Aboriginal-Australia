function[S,E,I,R]=RecoverDeath(S,E,I,R,alpha,gamma)
%I->R & R->D
%Created 3 May, 2024
Inew_E = binornd(E,alpha);   % E--> I
I = I + Inew_E;
E = E - Inew_E;

End_I = binornd(I,gamma); % I-->R
R = R + End_I;
I = I - End_I;
end
