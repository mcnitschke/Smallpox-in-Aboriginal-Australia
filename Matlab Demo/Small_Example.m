M = 4;                         % Total number of patches
E0 = zeros(M,1); 

load('A.mat')  
load('Preference.mat')
load('Q_small.mat')
load('Populations.mat')
load('Coordinates.mat')

Populations = round(Populations);

L = 5000;                        % Maximum time horizon          
NR = 1;                          % Number of Realisations
E0(2) = 1;                     % Initial infected population and location
rho     = [1/5,0,0,0,0]; 
alpha   = 1/12;       
gamma   = 1/18;  
s       = 6;
ft      = 0.25;
fc      = 0.8;
nu      = 0.75;
epsilon = 0.0;
lambda  = 0.1;
beta    = .334;
a    = 8;
b    = 1;

Ht = ((ft*(Populations/s)).^nu)./Populations;
Hc = ((fc*(Populations/s)).^nu)./Populations;


    mu = (a/365)*Ht + (b/730)*Hc; 
    Q_new = DistanceDecay(Q, Coordinates, lambda, epsilon,Preference); 

[PPI,TE,TI,TR,TNumPatches,TimeEnd,TotalUniqueE]=...
Middle_Example(NR,E0,A,M,L,mu,rho,alpha,gamma,beta,Q_new,Populations);

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [0 0 12 4]);
subplot(1,2,1)
plot(TE(1:TimeEnd),'LineWidth',1)
hold on
plot(TI(1:TimeEnd),'LineWidth',1)
plot(TR(1:TimeEnd),'LineWidth',1)
xlabel('time (days)')
ylabel('total individuals')
legend('Exposed','Infectious','Recovered','location','northwest')
title('Disease Dynamics')
subplot(1,2,2)
plot(TNumPatches(1:TimeEnd),'LineWidth',1)
title('Patches with Infectious Individuals')
xlabel('time (days)')
ylabel('total patches (max 4)')
