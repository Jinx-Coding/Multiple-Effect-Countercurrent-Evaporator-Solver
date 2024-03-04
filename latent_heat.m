
close all
clear variables
clc

TS = linspace(280,400.15,201);
lambdaS = zeros(1,length(TS));
hV = zeros(1,length(TS));

for i = 1 : length(TS)

    lambdaS(i) = 2500.304 -2.2521025*(TS(i)-273.15) - ...
                 0.021465847*(TS(i)-273.15)^1.5 + ...
                 3.1750136e-4*(TS(i)-273.15)^2.5 - ...
                 2.8607959e-5*(TS(i)-273.15)^3;

    hV(i) = exp((64.87678+11.76476*(log(647.096/TS(i)))^0.35 ...
           -11.94431*(647.096/TS(i))^2+6.29015*(647.096/TS(i))^3 ...
           -0.99893*(647.096/TS(i))^4)^0.5);
    
end

figure(1)
plot(TS-273.15,hV,'b')
hold on
plot(TS-273.15,lambdaS,'r')
legend('Steam Entalpy','Latent heat')
xlabel('T [C]')
ylabel('h [kJ/kg]')
grid on
box on


