function[S,E,I,R]=Spread_beta(S,E,I,R,A,dt,Bmat,tau)
% created May 3, 2024

M = size(S,1);  %Number of patches

SumI = repmat(sum(I),M,1);
Amat = repmat(A,1,M);
P = Bmat.*S.*SumI./Amat;

t = 0;

while t<dt
New_E = poissrnd(P.*tau); 

S = S - New_E;
E = E + New_E;
t = t + tau;
end
S(S<0)=0;  % To avoid overcounting...from tau-leap approx
end
