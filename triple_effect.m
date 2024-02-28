
% Triple Effect Evaporator Solver %
% Giovanni Correra 02/2024 %

clear variables
close all
clc

OPTIONS = optimset('Display','on','MaxIter',1e20,'MaxFunEvals',1e20...
    ,'TolFun',1e-10,'Algorithm','levenberg-marquardt');

% Data % 

TF = 45 + 273.15; % (K) %
omF = 0.1; % (kg/kg) %
TS = 100 + 273.15; % (K) %
om3 = 0.3; % (kg/kg) %
T3 = 80 + 273.15; % (K) %
F = 1000/36; % (kg/s) %

% Solution %

n0 = [TF+(T3-TF)/3,TF+2*(T3-TF)/3,0.75*F,0.5*F];
n = fsolve(@(n) parameters(n,F,TF,omF,TS,om3,T3),n0,OPTIONS);
[V1,V2,V3,L3,om1,om2,A,U1,U2,U3,S,np] = data(n,F,omF,TS,om3,T3,TF);
P0 = [saturation(n(1)),saturation(n(2)),saturation(T3)];
P1 = fsolve(@(P1) boiling_point(P1,n(1)-273.15,om1),P0(1),OPTIONS);
P2 = fsolve(@(P2) boiling_point(P2,n(2)-273.15,om2),P0(2),OPTIONS);
P3 = fsolve(@(P3) boiling_point(P3,T3-273.15,om3),P0(3),OPTIONS);

% Post Processing %

fprintf('A = %.2f (m2)      ',A); fprintf(' F = %.2f (kg/s)',F)
fprintf('   om,F = %.4f (-)  ',omF); fprintf('T,F = %.2f (°C)',TF-273.15)
fprintf('   S = %.2f (kg/s)\n',S)
fprintf('L1 = %.2f (kg/s)   ',n(3)); fprintf('V1 = %.2f (kg/s)   ',V1)
fprintf('om,1 = %.4f (-)  ',om1); fprintf('T,1 = %.2f (°C)',n(1)-273.15)
fprintf('   U1 = %.2f (W/m2K)',U1); fprintf('   P1 = %.2f (kPa)\n',P1)
fprintf('L2 = %.2f (kg/s)   ',n(4)); fprintf('V2 = %.2f (kg/s)   ',V2)
fprintf('om,2 = %.4f (-)  ',om2); fprintf('T,2 = %.2f (°C)',n(2)-273.15)
fprintf('   U2 = %.2f (W/m2K)',U2); fprintf('   P2 = %.2f (kPa)\n',P2)
fprintf('L3 = %.2f (kg/s)   ',L3); fprintf('V3 = %.2f (kg/s)   ',V3)
fprintf('om,3 = %.4f (-)  ',om3); fprintf('T,3 = %.2f (°C)',T3-273.15)
fprintf('   U3 = %.2f (W/m2K)',U3); fprintf('   P3 = %.2f (kPa)\n',P3)


% Functions %

function np = parameters(n,F,TF,omF,TS,om3,T3)

T1 = n(1);
T2 = n(2);
L1 = n(3);
L2 = n(4);

% Mass Balances %

L3 = F*omF/om3;
V3 = L2 - L3;
V1 = F - L1;
om1 = omF*F/L1;
V2 = L1 - L2;
om2 = om1*L1/L2;

% Energy Balances %

lambdaS = -(5.609e-3)*TS^2 + 1.355*TS + 2537;
U3 = (4.086*(T3-273.15)+72.6)/om3;
hL3 = (4.184-2.9337*om3)*(T3-273.15);
hL2 = (4.184-2.9337*om2)*(T2-273.15);
hV3 = exp((64.87678+11.76476*(log(647.096/T3))^0.35 ...
           -11.94431*(647.096/T3)^2+6.29015*(647.096/T3)^3 ...
           -0.99893*(647.096/T3)^4)^0.5);
S = (L3*hL3+V3*hV3-L2*hL2)/lambdaS;
A = S*lambdaS*1000/(U3*(TS-T3));

% Eq 1 %

lambda2 = -(5.609e-3)*T2^2 + 1.355*T2 + 2537;
hF = (4.184-2.9337*omF)*(TF-273.15);
hL1 = (4.184-2.9337*om1)*(T1-273.15);
hV1 = exp((64.87678+11.76476*(log(647.096/T1))^0.35 ...
           -11.94431*(647.096/T1)^2+6.29015*(647.096/T1)^3 ...
           -0.99893*(647.096/T1)^4)^0.5);

np(1) = V2*lambda2 + F*hF - L1*hL1 - V1*hV1;

% Eq 2 %

U1 = (4.086*(T1-273.15)+72.6)/om1;
hV2 = exp((64.87678+11.76476*(log(647.096/T2))^0.35 ...
           -11.94431*(647.096/T2)^2+6.29015*(647.096/T2)^3 ...
           -0.99893*(647.096/T2)^4)^0.5);

np(2) = (U1*A*(T2-T1)/1000 - V2*lambda2);

% Eq 3 %

lambda3 = -(5.609e-3)*T3^2 + 1.355*T3 + 2537;

np(3) = V3*lambda3 + L1*hL1 - L2*hL2 - V2*hV2;

% Eq 4 %

U2 = (4.086*(T2-273.15)+72.6)/om2;

np(4) = (U2*A*(T3-T2)/1000 - V3*lambda3);

end

function [V1,V2,V3,L3,om1,om2,A,U1,U2,U3,S,np] = data(n,F,omF,TS,om3,T3,TF)

T1 = n(1);
T2 = n(2);
L1 = n(3);
L2 = n(4);

% Mass Balances %

L3 = F*omF/om3;
V3 = L2 - L3;
V1 = F - L1;
om1 = omF*F/L1;
V2 = L1 - L2;
om2 = om1*L1/L2;

% Energy Balances %

lambdaS = -(5.609e-3)*TS^2 + 1.355*TS + 2537;
U3 = (4.086*(T3-273.15)+72.6)/om3;
hL3 = (4.184-2.9337*om3)*(T3-273.15);
hL2 = (4.184-2.9337*om2)*(T2-273.15);
hV3 = exp((64.87678+11.76476*(log(647.096/T3))^0.35 ...
           -11.94431*(647.096/T3)^2+6.29015*(647.096/T3)^3 ...
           -0.99893*(647.096/T3)^4)^0.5);
S = (L3*hL3+V3*hV3-L2*hL2)/lambdaS;
A = S*lambdaS*1000/(U3*(TS-T3));

% Eq 1 %

lambda2 = -(5.609e-3)*T2^2 + 1.355*T2 + 2537;
hF = (4.184-2.9337*omF)*(TF-273.15);
hL1 = (4.184-2.9337*om1)*(T1-273.15);
hV1 = exp((64.87678+11.76476*(log(647.096/T1))^0.35 ...
           -11.94431*(647.096/T1)^2+6.29015*(647.096/T1)^3 ...
           -0.99893*(647.096/T1)^4)^0.5);

np(1) = V2*lambda2 + F*hF - L1*hL1 - V1*hV1;

% Eq 2 %

U1 = (4.086*(T1-273.15)+72.6)/om1;
hV2 = exp((64.87678+11.76476*(log(647.096/T2))^0.35 ...
           -11.94431*(647.096/T2)^2+6.29015*(647.096/T2)^3 ...
           -0.99893*(647.096/T2)^4)^0.5);

np(2) = (U1*A*(T2-T1)/1000 - V2*lambda2);

% Eq 3 %

lambda3 = -(5.609e-3)*T3^2 + 1.355*T3 + 2537;

np(3) = V3*lambda3 + L1*hL1 - L2*hL2 - V2*hV2;

% Eq 4 %

U2 = (4.086*(T2-273.15)+72.6)/om2;

np(4) = (U2*A*(T3-T2)/1000 - V3*lambda3);


end

function P = saturation(T)

T = T-273.15;
P = 0.000135*T^3 - 0.00552*T^2 + 0.175918*T + 0.207651; % (kPa) %

end

function f = boiling_point(P,T,om)

Ts = 18.536*log(P) + 5.5052; % (C) %
f = T - Ts - (0.175*om^1.11)*exp(3.86*om)*(1000*P/101325)^0.2898;

end