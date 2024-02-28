# Triple Effect Evaporator Solver #
# Giovanni Correra 02/2024 #

import numpy
from scipy.optimize import fsolve

# Data #

TF = 45 + 273.15  # (K) #
omF = 0.1  # (kg/kg) #
TS = 100 + 273.15  # (K) #
om3 = 0.3  # (kg/kg) #
T3 = 80 + 273.15  # (K) #
F = 1000/36  # (kg/s) #

# Functions #

# Main function #


def parameters(n, feed, feedT, feedOm, steamT, finalOm, maxT):
    T1 = n[0]
    T2 = n[1]
    L1 = n[2]
    L2 = n[3]

    # Mass Balances #

    L3 = feed * feedOm / finalOm
    V3 = L2 - L3
    V1 = feed - L1
    om1 = feedOm * feed / L1
    V2 = L1 - L2
    om2 = om1 * L1 / L2

    # Energy Balances #

    lambdaS = -5.609e-3 * steamT ** 2 + 1.355 * steamT + 2537
    U3 = (4.086 * (maxT - 273.15) + 72.6) / finalOm
    hL3 = (4.184 - 2.9337 * finalOm) * (maxT - 273.15)
    hL2 = (4.184 - 2.9337 * om2) * (T2 - 273.15)
    hV3 = numpy.exp(numpy.sqrt(64.87678 + 11.76476 * (numpy.log(647.096 / maxT)) ** 0.35
                               - 11.94431 * (647.096 / maxT) ** 2 + 6.29015 * (647.096 / maxT) ** 3
                               - 0.99893 * (647.096 / maxT) ** 4))
    S = (L3 * hL3 + V3 * hV3 - L2 * hL2) / lambdaS
    A = S * lambdaS * 1000 / (U3 * (steamT - maxT))

    # Eq 1 #

    lambda2 = -5.609e-3 * T2 ** 2 + 1.355 * T2 + 2537
    hF = (4.184 - 2.9337 * feedOm) * (feedT - 273.15)
    hL1 = (4.184 - 2.9337 * om1) * (T1 - 273.15)
    hV1 = numpy.exp(numpy.sqrt(64.87678 + 11.76476 * (numpy.log(647.096 / T1)) ** 0.35
                               - 11.94431 * (647.096 / T1) ** 2 + 6.29015 * (647.096 / T1) ** 3
                               - 0.99893 * (647.096 / T1) ** 4))

    Z1 = V2 * lambda2 + feed * hF - L1 * hL1 - V1 * hV1

    # Eq 2 #

    U1 = (4.086 * (T1 - 273.15) + 72.6) / om1
    hV2 = numpy.exp(numpy.sqrt(64.87678 + 11.76476 * (numpy.log(647.096 / T2)) ** 0.35
                               - 11.94431 * (647.096 / T2) ** 2 + 6.29015 * (647.096 / T2) ** 3
                               - 0.99893 * (647.096 / T2) ** 4))

    Z2 = (U1 * A * (T2 - T1) / 1000 - V2 * lambda2)

    # Eq 3 #

    lambda3 = -5.609e-3 * maxT ** 2 + 1.355 * maxT + 2537

    Z3 = V3 * lambda3 + L1 * hL1 - L2 * hL2 - V2 * hV2

    # Eq 4 #

    U2 = (4.086 * (T2 - 273.15) + 72.6) / om2

    Z4 = (U2 * A * (maxT - T2) / 1000 - V3 * lambda3)

    Z = numpy.array((Z1, Z2, Z3, Z4))

    return Z


# Return values function #


def data(n, feed, feedOm, steamT, finalOm, finalT, feedT):

    T1 = n[0]
    T2 = n[1]
    L1 = n[2]
    L2 = n[3]

    # Mass balances #

    L3 = feed * feedOm / finalOm
    V3 = L2 - L3
    V1 = feed - L1
    om1 = feedOm * feed / L1
    V2 = L1 - L2
    om2 = om1*L1/L2

    # Energy Balances #

    lambdaS = -5.609e-3 * steamT ** 2 + 1.355 * steamT + 2537
    U3 = (4.086 * (finalT - 273.15) + 72.6) / finalOm
    hL3 = (4.184 - 2.9337 * finalOm) * (finalT - 273.15)
    hL2 = (4.184-2.9337*om2)*(T2-273.15)
    hV3 = numpy.exp((64.87678 + 11.76476 * (numpy.log(647.096 / finalT)) ** 0.35
                     - 11.94431 * (647.096 / finalT) ** 2 + 6.29015 * (647.096 / finalT) ** 3
                     - 0.99893 * (647.096 / finalT) ** 4) ** 0.5)
    S = (L3*hL3+V3*hV3-L2*hL2)/lambdaS
    A = S*lambdaS*1000/(U3 * (steamT - finalT))

    # Eq 1 #

    lambda2 = -5.609e-3*T2**2 + 1.355*T2 + 2537
    hF = (4.184 - 2.9337 * feedOm) * (feedT - 273.15)
    hL1 = (4.184-2.9337*om1)*(T1-273.15)
    hV1 = numpy.exp((64.87678+11.76476*(numpy.log(647.096/T1))**0.35
                     - 11.94431*(647.096/T1)**2+6.29015*(647.096/T1)**3
                     - 0.99893*(647.096/T1)**4)**0.5)

    Z1 = V2 * lambda2 + feed * hF - L1 * hL1 - V1 * hV1

    # Eq 2 #

    U1 = (4.086*(T1-273.15)+72.6)/om1
    hV2 = numpy.exp((64.87678+11.76476*(numpy.log(647.096/T2))**0.35
                     - 11.94431*(647.096/T2)**2+6.29015*(647.096/T2)**3
                     - 0.99893*(647.096/T2)**4)**0.5)

    Z2 = (U1*A*(T2-T1)/1000 - V2*lambda2)

    # Eq 3 #

    lambda3 = -5.609e-3 * finalT ** 2 + 1.355 * finalT + 2537

    Z3 = V3*lambda3 + L1*hL1 - L2*hL2 - V2*hV2

    # Eq 4 #

    U2 = (4.086*(T2-273.15)+72.6)/om2

    Z4 = (U2 * A * (finalT - T2) / 1000 - V3 * lambda3)

    Z = numpy.array((Z1, Z2, Z3, Z4))

    return Z, V1, V2, V3, L3, om1, om2, A, U1, U2, U3, S


# Saturation Pressure first guess #

def saturation(tmp):

    tmp = tmp - 273.15
    P = 0.000135*tmp**3 - 0.00552*tmp**2 + 0.175918*tmp + 0.207651  # (kPa) #

    return P


# Tomato Paste saturation pressure #

def boiling_point(P, T, om):

    Ts = 18.536*numpy.log(P) + 5.5052  # (C) #
    f = T - Ts - (0.175*om**1.11)*numpy.exp(3.86*om)*(1000*P/101325)**0.2898

    return f


# Solution #

sol0 = numpy.array((TF+(T3-TF)/3, TF+2*(T3-TF)/3, 0.75*F, 0.5*F))
sol = fsolve(lambda n: parameters(n, F, TF, omF, TS, om3, T3), sol0)
temp1 = float(sol[0])
temp2 = float(sol[1])
liquid1 = float(sol[2])
liquid2 = float(sol[3])
system_parameters = data(sol, F, omF, TS, om3, T3, TF)
zeros = system_parameters[0]
P0 = numpy.array((saturation(temp1), saturation(temp2), saturation(T3)))
pressure1 = (fsolve(lambda P1: boiling_point(P1, temp1-273.15, system_parameters[5]), P0[0]))
pressure2 = (fsolve(lambda P2: boiling_point(P2, temp2-273.15, system_parameters[6]), P0[1]))
pressure3 = (fsolve(lambda P3: boiling_point(P3, T3-273.15, om3), P0[2]))

# Post Processing #

print('\n')
print('zero functions = ', end="")
print(zeros)
print('F  = ', end="")
print('%.2f' % F, end=" [kg/s]")
print('   S  = ', end="")
print('%.2f' % system_parameters[11], end=" [kg/s]")
print('   omF = ', end="")
print('%.3f' % omF, end=" [kg/kg]")
print('   TF = ', end="")
print('%.2f' % (TF-273.15), end=" [C]")
print('   A  = ', end="")
print('%.2f' % system_parameters[7], end=" [m2]\n")
print('L1 = ', end="")
print('%.2f' % liquid1, end=" [kg/s]")
print('   V1 = ', end="")
print('%.2f' % system_parameters[1], end=" [kg/s]")
print('   om1 = ', end="")
print('%.3f' % system_parameters[5], end=" [kg/kg]")
print('   T1 = ', end="")
print('%.2f' % (temp1-273.15), end=" [C]")
print('   P1 = ', end="")
print('%.2f' % float(pressure1[0]), end="  [kPa]")
print('   U1 = ', end="")
print('%.2f' % system_parameters[8], end=" [W/m2K]\n")
print('L2 = ', end="")
print('%.2f' % liquid2, end=" [kg/s]")
print('   V2 = ', end="")
print('%.2f' % system_parameters[2], end=" [kg/s]")
print('   om2 = ', end="")
print('%.3f' % system_parameters[6], end=" [kg/kg]")
print('   T2 = ', end="")
print('%.2f' % (temp2-273.15), end=" [C]")
print('   P2 = ', end="")
print('%.2f' % float(pressure2[0]), end="  [kPa]")
print('   U2 = ', end="")
print('%.2f' % system_parameters[9], end=" [W/m2K]\n")
print('L2 =  ', end="")
print('%.2f' % system_parameters[4], end=" [kg/s]")
print('   V2 = ', end="")
print('%.2f' % system_parameters[3], end=" [kg/s]")
print('   om3 = ', end="")
print('%.3f' % om3, end=" [kg/kg]")
print('   T3 = ', end="")
print('%.2f' % (T3-273.15), end=" [C]")
print('   P3 = ', end="")
print('%.2f' % float(pressure3[0]), end="  [kPa]")
print('   U2 = ', end="")
print('%.2f' % system_parameters[10], end=" [W/m2K]\n")