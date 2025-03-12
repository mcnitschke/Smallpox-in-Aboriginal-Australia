M = 592;                         % Total number of patches
E0 = zeros(M,1); 

load('Q_Estimate.mat')  
load('Tindale.mat')
load('Population_Estimate.mat')

Populations = SPM*1.18;          % Population estimates
Q = Q_Estimate;                  % Raw probability between-patch movement 
A = Tindale.area(:,1);           % Area of each patch
L = 10000;                       % Maximum time horizon          

mu = zeros(5,M);
NR = 10;                         % Number of Realisations

E0(518) = 1;                     % Initial infected population and location
PopDist = round(Populations)';
rho = [1/5,0,0,0,0]; 
alpha = 1/12;       
gamma = 1/18;  
beta = .334;
avec = [50];
bvec = [125];

for i = 1:length(bvec);
    for j = 1:length(avec);
    a = avec(j);
    b = bvec(i);
    mu(1,:) = (1/365)*(a./PopDist)+(1/730)*(b./PopDist); 

[PPI,TI,TNumPatches,TimeEnd,TotalUniqueE]=...
Middle(NR,E0,A,M,L,mu,rho,alpha,gamma,beta,Q,PopDist,a,b);

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
