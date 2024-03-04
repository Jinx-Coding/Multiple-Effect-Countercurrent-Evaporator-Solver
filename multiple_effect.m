
% Multiple Effect Evaporator Solver %
% Giovanni Correra 03/2024 %

clear variables
close all
clc

OPTIONS = optimset('Display','off','MaxIter',1e20,'MaxFunEvals',1e20...
    ,'TolFun',1e-10,'Algorithm','levenberg-marquardt');

% ----------------------------- Data ----------------------------------- % 

TF = 95 + 273.15; % (K) %
omF = 0.05; % (kg/kg) %
TS = 110 + 273.15; % (K) %
omf = 0.3; % (kg/kg) %
Tf = 80 + 273.15; % (K) %
F = 98*1000/3600; % (kg/s) %

% ---------------------------- Solution -------------------------------- %

N = 4;
n0 = first_guess(TF,Tf,F,N,omF,omf);
n = fsolve(@(n) parameters(n,F,TF,omF,TS,omf,Tf,N),n0,OPTIONS);
[T,V,L,om,A,U,S,dT,np] = data(n,F,omF,TS,omf,Tf,TF,N);

SE = sum(V)/S(N);
Atot = N*A;

Pop = zeros(1,N);
for i = 1 : N

    P0 = saturation(T(i));
    P = fsolve(@(P) boiling_point(P,T(i)-273.15,om(i)),P0,OPTIONS);
    Pop(i) = P;

end

% ------------------------- Post - processing -------------------------- %

postprocessing(V,L,om,A,U,S,np,F,omF,TS,TF,N,SE,Atot,T,Pop,dT);

% ---------------------------- Functions ------------------------------- %

function n0 = first_guess(TF,Tf,F,N,omF,omf)

n01 = zeros(1,N-1);
n02 = zeros(1,N-1);

Lf = F*omF/omf;
L0 = (F-Lf)/(N);
n02(1) = F-L0;

for i = 1 : N-1

    n01(i) = TF + i*(Tf-TF)/(N-1);
    n02(i) = F - i*L0;
    
end

n0 = [n01,n02];

end

function np = parameters(n,F,TF,omF,TS,omf,Tf,N)

% Mass Balances %

Lf = F*omF/omf;
L = [n(1,N:2*(N-1)),Lf];
F = [F,L(1,1:N-1)];
V = zeros(1,N);
om = zeros(1,N);

for i = 1 : N
    V(i) = F(i) - L(i); 
end

om(1) = omF*F(1)/L(1);

for i = 1 : N-1
    om(i+1) = om(i)*F(i)/L(i);
end

om(N) = omf;

% Energy Balances %

T = [n(1,1:N-1),Tf];
U = zeros(1,N);

for i = 1 : N
    U(i) = (4.086*(T(i)-273.15)+72.6)/om(i);
end

hL = zeros(1,N);
hV = zeros(1,N);

for i = 1 : N
    hL(i) = (4.184-2.9337*om(i))*(T(i)-273.15);
    hV(i) = exp((64.87678+11.76476*(log(647.096/T(i)))^0.35 ...
           -11.94431*(647.096/T(i))^2+6.29015*(647.096/T(i))^3 ...
           -0.99893*(647.096/T(i))^4)^0.5);
end

hF0 = (4.184-2.9337*omF)*(TF-273.15);
hF = [hF0,hL(1,1:N-1)];

Ts = [T(1,2:N),TS];

lambda = zeros(1,N);

for i = 1 : N
    lambda(i) = 2500.304 -2.2521025*(Ts(i)-273.15) - ...
                 0.021465847*(Ts(i)-273.15)^1.5 + ...
                 3.1750136e-4*(Ts(i)-273.15)^2.5 - ...
                 2.8607959e-5*(Ts(i)-273.15)^3;
end

S = zeros(1,N);
S(N) = (L(N)*hL(N)+V(N)*hV(N)-L(N-1)*hL(N-1))/lambda(N);

for i = 1 : N-1
    S(i) = V(i+1);
end

A = S(N)*lambda(N)*1000/(U(N)*(Ts(N)-T(N)));

% Solving equations %

np1 = zeros(1,N-1);
np2 = zeros(1,N-1);

for i = 1 : (N-1)

    np1(i) = F(i)*hF(i) + S(i)*lambda(i) - L(i)*hL(i) - V(i)*hV(i); 
    np2(i) = U(i)*A*(Ts(i)-T(i))*1e-3 - S(i)*lambda(i);

end

np = [np1,np2];

end

function [T,V,L,om,A,U,S,dT,np] = data(n,F,omF,TS,omf,Tf,TF,N)

% Mass Balances %

Lf = F*omF/omf;
L = [n(1,N:2*(N-1)),Lf];
F = [F,L(1,1:N-1)];
V = zeros(1,N);
om = zeros(1,N);

for i = 1 : N
    V(i) = F(i) - L(i); 
end

om(1) = omF*F(1)/L(1);

for i = 1 : N-1
    om(i+1) = om(i)*F(i)/L(i);
end

om(N) = omf;

% Energy Balances %

T = [n(1,1:N-1),Tf];
U = zeros(1,N);

for i = 1 : N
    U(i) = (4.086*(T(i)-273.15)+72.6)/om(i);
end

hL = zeros(1,N);
hV = zeros(1,N);

for i = 1 : N
    hL(i) = (4.184-2.9337*om(i))*(T(i)-273.15);
    hV(i) = exp((64.87678+11.76476*(log(647.096/T(i)))^0.35 ...
           -11.94431*(647.096/T(i))^2+6.29015*(647.096/T(i))^3 ...
           -0.99893*(647.096/T(i))^4)^0.5);
end

hF0 = (4.184-2.9337*omF)*(TF-273.15);
hF = [hF0,hL(1,1:N-1)];

Ts = [T(1,2:N),TS];

lambda = zeros(1,N);

for i = 1 : N
    lambda(i) = 2500.304 -2.2521025*(Ts(i)-273.15) - ...
                 0.021465847*(Ts(i)-273.15)^1.5 + ...
                 3.1750136e-4*(Ts(i)-273.15)^2.5 - ...
                 2.8607959e-5*(Ts(i)-273.15)^3;
end

S = zeros(1,N);
S(N) = (L(N)*hL(N)+V(N)*hV(N)-L(N-1)*hL(N-1))/lambda(N);

for i = 1 : N-1
    S(i) = V(i+1);
end

A = S(N)*lambda(N)*1000/(U(N)*(Ts(N)-T(N)));

% Solving equations %

np1 = zeros(1,N-1);
np2 = zeros(1,N-1);

for i = 1 : (N-1)

    np1(i) = F(i)*hF(i) + S(i)*lambda(i) - L(i)*hL(i) - V(i)*hV(i); 
    np2(i) = U(i)*A*(Ts(i)-T(i))*1e-3 - S(i)*lambda(i);

end

dT = zeros(1,N);

for i = 1 : N

    dT(i) = Ts(i)-T(i);

end

np = [np1,np2];

end

function P = saturation(T)

T = T-273.15;
P = 0.000135*T^3 - 0.00552*T^2 + 0.175918*T + 0.207651; % (kPa) %

end

function f = boiling_point(P,T,om)

Ts = 18.536*log(P) + 5.5052; % (C) %
f = T - Ts - (0.175*om^1.11)*exp(3.86*om)*(1000*P/101325)^0.2898;

end

function postprocessing(V,L,om,A,U,S,np,F,omF,TS,TF,N,SE,Atot,T,P,dT)

fprintf(['zero functions : ' repmat(' %e ',1,numel(np)) ' [-]\n'],np);
fprintf('n. effects = %.0f [-]',N)
fprintf('  A = %.2f [m2]',A)
fprintf('  A,tot = %.2f [m2]\n', Atot)
fprintf('F = %.2f [kg/s]   ',F)
fprintf('  T,F = %.2f [°C]',TF-273.15)
fprintf('  om,F = %.2f %% [-]\n',omF*100)
fprintf('S = %.2f [kg/s]',S(N))
fprintf('  T,S = %.2f [°C]',TS-273.15)
fprintf('  Steam economy = %.2f\n',SE)

for i = 1 : N

    fprintf('V,'); fprintf('%d',i); fprintf(' = %.2f [kg/s]  ',V(i))
    fprintf('L,'); fprintf('%d',i); fprintf(' = %.2f [kg/s]  ',L(i))
    fprintf('om,'); fprintf('%d',i); fprintf(' = %.2f %% [-]  ',om(i)*100)
    fprintf('T,'); fprintf('%d',i); fprintf(' = %.2f [°C]  ',T(i)-273.15)
    fprintf('dT,'); fprintf('%d',i); fprintf(' = %.2f [°C]   ',dT(i))
    fprintf('U,'); fprintf('%d',i); fprintf(' = %.2f [W/m2K]  ',U(i))
    fprintf('P,'); fprintf('%d',i); fprintf(' = %.2f [kPa]  ',P(i))
    fprintf('\n')

end

end


