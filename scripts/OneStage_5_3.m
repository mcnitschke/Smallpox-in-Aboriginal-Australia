M = 592;  %Total number of patches
E0 = zeros(M,1); 

load('Q.mat')  % Q - probability of movement from node to node.
load('Tindale.mat')
load('Populations.mat')
A = Tindale.area(:,1);  % area in square miles
L = 5000;

mu = zeros(5,M);
NR = 10; %Number of Realisations for each parameter combination

E0(5) = 1;  %Initial infected population and location
PopDist = round(Populations)';
rho = [1/5,0,0,0,0]; 
alpha = 1/12;    
gamma = 1/18;  
beta = .334;
% avec = [10 15 20 25 30];      %regular intra-patch interactions (trade)
% bvec = [25 50 75 100 125];     %less regular intra-patch interactions (gatherings)
avec = [50];
bvec = [125];
tau = 1/100;

for i = 1:length(bvec);
    for j = 1:length(avec);
    a = avec(j);
    b = bvec(i);

for k = 1:M  %movement depends on population
    mu(1,k) = (1/365)*(a/PopDist(k))+(1/730)*(b/PopDist(k));
end
[PPI,TI,TNumPatches,TimeEnd,TotalUniqueE]=...
ABM1StageC(NR,E0,A,M,L,mu,rho,alpha,gamma,beta,tau,Q,PopDist,a,b);

MeanTIMatrix(j,:) = mean(TI); %1xt vector
FrequencyPP(j,:) = sum(PPI);
MeanTimeEndVec(j) = median(TimeEnd);
MeanTUniqueIVec(j) = mean(TotalUniqueE);
MaxPatchesVec(j) = median(max(TNumPatches,[],2));

end
MedianTimeEndMat(i,:) = MeanTimeEndVec;
MeanTUniqueIMat(i,:) = MeanTUniqueIVec;
MedianMaxSpread(i,:) = MaxPatchesVec;
MeanTIArray(:,:,i) = MeanTIMatrix;
FrequencyPPArray(:,:,i) = FrequencyPP;
end