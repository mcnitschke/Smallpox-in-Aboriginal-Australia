function[IPP,TI,TNumPatches,TimeEnd,TotalUniqueE]=Middle(NR,E0,A,M,L,mu,rho,alpha,gamma,beta,Q,PopDist,a,b);
% Created May, 3 2024
% Simplified global model ....No age structure....4 stages of Infectious
% mu - probability of moving out
% sigma - probability of visiting a patch (distributed over available)
% rho- probability of return in
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = 1; % time step (days)
%Initialize Populations %%%%%%%%%%%%%%%

Bmat = repmat(beta,M,M);
Amat = repmat(A,1,M);

TTNewS = zeros(NR,L);
TTNewE = zeros(NR,L);
TTNewI = zeros(NR,L);
TTNewR = zeros(NR,L);
TTNewD = zeros(NR,L);

TPP_S = zeros(L,M,NR);
TPP_E = zeros(L,M,NR);
TPP_I = zeros(L,M,NR);
TPP_R = zeros(L,M,NR);
TPP_D = zeros(L,M,NR);
TNumPatches = zeros(NR,L);

for j = 1:NR;

S = zeros(M,M);
E = S;
I = S;
R = S;
Death = S;

TNewS = zeros(1,L);
TNewE = zeros(1,L);
TNewI = zeros(1,L);
TNewR = zeros(1,L);
TNewD = zeros(1,L);

PPNewS = zeros(L,M);
PPNewE = zeros(L,M);
PPNewI = zeros(L,M);
PPNewR = zeros(L,M);
PPNewD = zeros(L,M);
%%%%%%%Initial Seed Population %%%%%%%%%%%%%%%%%%%%%%%%%%%%
S = diag(PopDist);
E = diag(E0); 
%%%%%%%%%%%%%%Sigma & Rho --visit & return rates%%%%%%%%%%%%%%%%%%%%%%
MuRho1 = rho(1)*ones(M,M); % Healthy 
MuRho2 = rho(2)*ones(M,M); % I1
MuRho3 = rho(3)*ones(M,M); % I2
MuRho4 = rho(4)*ones(M,M); % I3
MuRho5 = rho(5)*ones(M,M); % I4

out_1 = mu(1,:); % move out rate --healthy 
out_2 = mu(2,:); % move out rate --I1
out_3 = mu(3,:); % move out rate --I2
out_4 = mu(4,:); % move out rate --I3
out_5 = mu(5,:); % move out rate --I4

    MuRho1(:,:)=MuRho1(:,:)-diag(diag(MuRho1))+diag([out_1]);
    MuRho2(:,:)=MuRho2(:,:)-diag(diag(MuRho2))+diag([out_2]);
    MuRho3(:,:)=MuRho3(:,:)-diag(diag(MuRho3))+diag([out_3]);
    MuRho4(:,:)=MuRho4(:,:)-diag(diag(MuRho4))+diag([out_4]);
    MuRho5(:,:)=MuRho5(:,:)-diag(diag(MuRho5))+diag([out_5]);  
    TE = 0;
for t =1:L
t
 [S,E,I,R,TE]=Spread(S,E,I,R,TE,Amat,dt,Bmat);

 [S,E,I,R] = RecoverDeath(S,E,I,R,alpha,gamma);

 [S,E,I,R] = Movement(S,E,I,R,MuRho1,MuRho2,Q);

NumPatches(t,:) = nnz(sum(I));

PPNewS(t,:) = sum(S);
PPNewE(t,:) = sum(E);
PPNewI(t,:) = sum(I);
PPNewR(t,:) = sum(R);

TNewS(t) = sum(sum(S));
TNewE(t) = sum(sum(E));
TNewI(t) = sum(sum(I));
TNewR(t) = sum(sum(R));

if (sum(sum(R))>=1) && (sum(sum(I))<1)
    break
end

end
TotalUniqueE(j) = TE;
TimeEnd(j) = t;
TNumPatches(j,1:length(NumPatches)) =NumPatches';

TI(j,1:length(TNewI)) = TNewI; 
SPPI = sum(PPNewI);
SPPI(SPPI~=0)=1;
IPP(j,:) = SPPI;

end
end
