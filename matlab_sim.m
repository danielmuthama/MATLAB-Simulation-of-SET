q=1.602e-19; % electronic charge (C)
kb=1.381e-23; % Boltzman constant (J/K) % Definition of Device parameters
c1=1.0e-20; % tunnel capacitor C1 (F)
c2=2.1e-19; % tunnel capacitor C2 (F)
cg=1.0e-18; % gate capacitor Cg (F)
ctotal=c1+c2+cg; % total capacitance (F)
mega=1000000; % definition of mega=106
r1=15*mega; % tunnel resistance R1 (Ohm)
r2=250*mega; % tunnel resistance R2 (Ohm)
%PART ONE COMPLETE

Vg=0; % gate voltage (V)
q0=0; % background charge q0 is assumed to be zero
temp=10; % temperature T (K)
vmin=-0.5; % drain voltage minimum Vmin (V)
vmax=0.5; % drain voltage maximum Vmax (V)
NV=1000; % number of grid from Vmin to Vmax
dV=(vmax-vmin)/NV; % drain voltage increment of each grid point

for iv=1:NV % loop start for drain voltage
%PART TWO COMPLETE
V(iv)=vmin+iv*dV; % drain voltage in each grid point%
Nmin=-20; % minimum number of N (charge number in dot)
Nmax=20; % maximum number of N (charge number in dot)
for ne=1:Nmax-Nmin % loop start for N 
n=Nmin+ne; % N charge number in dot
dF1p=q/ctotal*(0.5*q+(n*q-q0)-(c2+cg)*V(iv)+cg*Vg); 
dF1n=q/ctotal*(0.5*q-(n*q-q0)+(c2+cg)*V(iv)-cg*Vg); 
dF2p=q/ctotal*(0.5*q-(n*q-q0)-c1*V(iv)-cg*Vg); 
dF2n=q/ctotal*(0.5*q+(n*q-q0)+c1*V(iv)+cg*Vg); 
% Noted that loop end for N is located after calculation of 
%PART THREE COMPLETE

if dF1p<0 
    T1p(ne)=1/(r1*q*q)*(-dF1p)/(1-exp(dF1p/(kb*temp))); %? positive in equation
else
    T1p(ne)=1e-1; %  positive is assumed to be very small
end
if dF1n<0 
    T1n(ne)=1/(r1*q*q)*(-dF1n)/(1-exp(dF1n/(kb*temp))); % ? negative in equation
else
    T1n(ne)=1e-1; %  negative is assumed to be very small
end
if dF2p<0 
    T2p(ne)=1/(r2*q*q)*(-dF2p)/(1-exp(dF2p/(kb*temp))); % ? positive in equation (26b) 
else 
    T2p(ne)=1e-1; %  positive is assumed to be very small
end
if dF2n<0 
    T2n(ne)=1/(r2*q*q)*(-dF2n)/(1-exp(dF2n/(kb*temp))); % negative equation
else 
    T2n(ne)=1e-1; %  negative is assumed to be very small
end
end
%PART FOUR COMPLETE

p(1)=0.001; % ?(Nmin) is assumed to be 0.01
p(Nmax-Nmin)=0.001; % ?(Nmax) is assumed to be 0.01Sixth, normalization of ? is done. Here, ?????????? is calculated. 

sum=0; % sum=0 is initial value to calculate ?
for ne=2:Nmax-Nmin 
    p(ne)=p(ne-1)*(T2n(ne-1)+T1p(ne-1))/(T2p(ne)+T1n(ne)); % calculation of ?(N) in equation (28)% The conditions below are used to avoid divergence of Matlab calculation
if p(ne)>1e250 
p(ne)=1e250; 
end
if p(ne)<1e-250 
    p(ne)=1e-250; 
end
% --------------------- 
sum=sum+p(ne); 
end
for ne=2:Nmax-Nmin 
    p(ne)=p(ne)/sum; % Normalization in equation
end
%PART FIVE COMPLETE

sumI=0; % sumI=0 is initial condition for current calculation
for ne=2:Nmax-Nmin 
    sumI=sumI+p(ne)*(T2p(ne)-T2n(ne)); 
end
I(iv)=q*sumI; % I in equation of drain voltage loop
end
plot(V,I); % plot of I vs V
for iv=1:NV-1 
dIdV(iv)=(I(iv+1)-I(iv))/dV; % calculation of dIdV
end
figure; 
plot(V(1,1:NV-1),dIdV); % plot of dIdV vs V
%PART SIX COMPLETE