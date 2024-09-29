function[S,E,I,R] = ...
   Movement1Stage(S,E,I,R,MuRho1,MuRho2,MuRho3,MuRho4,MuRho5,Q)
%Created 3 May, 2024
% All individuals same...but four Infectious stages
%Preallocate 
Destination = Q./sum(Q,2);
% Movers in and out of region
S_movers = binornd(S,MuRho1);
E_movers = binornd(E,MuRho1);
I_movers = binornd(I,MuRho2);
R_movers = binornd(R,MuRho1);

%Note: Returning individuals are easy as we just determine if they "return"
%to their resident patches.  However, we need an extra step when they leave
%since we first decide to leave, then those must be divided into their
%destinations according to a set distribution (multinomial distribiton).

%%%% reshuffle populations %%%%%%
% Determine where they go...
    vs = diag(S_movers); %number of visitors in each patch
    ve = diag(E_movers);
    vi = diag(I_movers);
    vr = diag(R_movers);
    
    S_visit = mnrnd(vs,Destination);  % where visitors go...
    E_visit = mnrnd(ve,Destination);
    I_visit = mnrnd(vi,Destination);
    R_visit = mnrnd(vr,Destination);
    
    S_return = diag(sum((S_movers-diag(diag(S_movers))),2));
    E_return = diag(sum((E_movers-diag(diag(E_movers))),2));
    I_return = diag(sum((I_movers-diag(diag(I_movers))),2));
    R_return = diag(sum((R_movers-diag(diag(R_movers))),2));

% Put them in place...
S = S - S_movers + S_visit + S_return;
E = E - E_movers + E_visit + E_return;
I = I - I_movers + I_visit + I_return;
R = R - R_movers + R_visit + R_return;

end