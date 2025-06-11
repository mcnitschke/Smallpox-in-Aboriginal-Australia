M = 592;                         % Total number of patches
E0 = zeros(M,1); 

load('Q_Estimate.mat')  
load('Tindale.mat')
load('Population_Estimate.mat')
load('PatchProbs.mat')

Populations = SPM*1.18;          % Population estimates
Q = Q_Estimate;                  % Raw probability between-patch movement 
A = Tindale.area(:,1);           % Area of each patch
L = 100;                        % Maximum time horizon          
tic
mu = zeros(5,M);
NR = 1;                          % Number of Realisations
E0(518) = 1;                     % Initial infected population and location
PopDist = round(Populations)';
rho     = [1/5,0,0,0,0]; 
alpha   = 1/12;       
gamma   = 1/18;  
Pref    = PatchProbs; % Node probabilities
Coordinates = Tindale.coords(:,1:2);
s       = 6;
ft      = 0.25;
fc      = 0.8;
nu      = 0.75;
epsilon = 0.01;
lambda  = 0.1;
beta    = .334;
avec    = [8];
bvec    = [1];

Ht = ((ft*(Populations/s)).^nu)./Populations;
Hc = ((fc*(Populations/s)).^nu)./Populations;

for i = 1:length(bvec)
    for j = 1:length(avec)
    a = avec(j);
    b = bvec(i);
    mu = (a/365)*Ht + (b/730)*Hc; 
    Q_new = DistanceDecay(Q, Coordinates, lambda, epsilon,Pref); 

[PPI,TI,TNumPatches,TimeEnd,TotalUniqueE]=...
Middle(NR,E0,A,M,L,mu,rho,alpha,gamma,beta,Q_new,PopDist);

MeanTIMatrix(j,:) = mean(TI); 
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
toc