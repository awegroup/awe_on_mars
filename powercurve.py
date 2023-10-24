#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 22:37:27 2023

Power curve computed following the three-regimes strategy.

The accuracy of the wind speeds delimiting the three regimes is set by the
wind_speed_delta parameter.

@author: Roland Schmehl
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import optimize as op
from pylab import np
mpl.rcParams['font.family'] = "Open Sans"
mpl.rcParams.update({'font.size': 14})
mpl.rcParams['figure.figsize'] = 10, 5.625
mpl.rc('xtick', labelsize=14)
mpl.rc('ytick', labelsize=14)
mpl.rcParams['pdf.fonttype'] = 42 # Output Type 3 (Type3) or Type 42 (TrueType)

# Environmental properties
atmosphere_density        =  0.01    # kg/**3
wind_speed_min            =  0.      # m/s
wind_speed_max            =  50.     # m/s
wind_speed_delta          =  0.5     # m/s     ! accuracy of wind speed regimes

# Kite properties
kite_planform_area       =  200.    # m**2
kite_lift_coefficient_out =  0.71    # -
kite_drag_coefficient_out =  0.14    # -
kite_lift_coefficient_in  =  0.71    # -
kite_drag_coefficient_in  =  0.07    # -

# Tether properties
nominal_tether_force      =  5100.   # N
tether_drag_coefficient   =  1.1     # -
tether_diameter           =  0.00484 # m

# Generator properties
nominal_generator_power   =  77000.  # W

# Operational parameters
elevation_angle_out       =  25.     # deg
elevation_angle_in        =  60.     # deg
reeling_speed_min_limit   = -21.     # m/s
reeling_speed_max_limit   =   8.     # m/s

# Derived properties
lift_to_drag_ratio_out = kite_lift_coefficient_out / kite_drag_coefficient_out
E2 = lift_to_drag_ratio_out**2
cosine_beta_out  = np.cos(np.radians(elevation_angle_out))
cosine_beta_in   = np.cos(np.radians(elevation_angle_in))
force_factor_out = kite_lift_coefficient_out * np.sqrt(1+1/E2) * (1+E2)
power_factor_ideal = force_factor_out * cosine_beta_out**3 * 4/27

wind_speed_range = wind_speed_max - wind_speed_min
num_wind_speeds  = int(wind_speed_range/wind_speed_delta + 1)
wind_speed  = np.linspace(wind_speed_min, wind_speed_max, num_wind_speeds)
cycle_power = []
power_ideal = []

# Objective function for the three wind speed domains
def objective_function_1(x):
    f_out    = x[0]
    f_in     = x[1]
    a        = 1 - 2*f_in*cosine_beta_in + f_in**2
    gamma_in = kite_lift_coefficient_in * np.sqrt(a/(1 - cosine_beta_in**2))
    return -((cosine_beta_out - f_out)**2 - (gamma_in / force_factor_out) * \
                a) * (f_in*f_out) / (f_in - f_out)

def objective_function_2(x, mu_F, f_nF):
    f_in     = x[0]
    a        = 1 - 2*f_in*cosine_beta_in + f_in**2
    b        = (mu_F - 1) * cosine_beta_out + f_nF
    gamma_in = kite_lift_coefficient_in * np.sqrt(a/(1 - cosine_beta_in**2))
    return -(((cosine_beta_out - f_nF) / mu_F)**2  \
             - (gamma_in / force_factor_out) * a)  \
             * f_in*b/(mu_F*f_in-b)

def objective_function_3(x, mu_P, f_nP):
    f_in     = x[0]
    a        = 1 - 2*f_in*cosine_beta_in + f_in**2
    gamma_in = kite_lift_coefficient_in * np.sqrt(a/(1 - cosine_beta_in**2))
    return -(((cosine_beta_out - f_nP) / mu_P)**2  \
             - (gamma_in / force_factor_out) * a)  \
             * f_in*f_nP/(mu_P*f_in-f_nP)

print("num_wind_speeds = ", num_wind_speeds)

###############################################################################

# Constant properties
rho = atmosphere_density
S   = kite_planform_area
CLo = kite_lift_coefficient_out

# Initialize wind speed regimes
wind_speed_regime      = 1
wind_speed_force_limit = 0
wind_speed_power_limit = 0

# Loop over wind speed range
for v_w in wind_speed:

    # Dynamic pressure
    q  = 0.5 * atmosphere_density * v_w**2

    # Wind power density
    Pw = q*v_w

    # Reeling factor limits
    f_max = reeling_speed_max_limit / v_w
    f_min = reeling_speed_min_limit / v_w

    # Unconstrained operation
    if wind_speed_regime == 1:

        starting_point = (0.001, -0.001)
        bounds         = ((0.001,  f_max), (f_min, -0.001),)

        optimisation_result = op.minimize(objective_function_1, \
                                          starting_point,       \
                                          bounds=bounds,        \
                                          method='SLSQP')

        # Reeling factors
        f_out = optimisation_result['x'][0]
        f_in  = optimisation_result['x'][1]

        # Normalized cycle power
        p_c = -objective_function_1 ([f_out, f_in])

        # Tether force and mechanical power during reel out
        Ft_out = q * kite_planform_area * force_factor_out * \
                     (cosine_beta_out - f_out)**2
        P_out  = Ft_out * v_w * f_out

        if Ft_out > nominal_tether_force:
            wind_speed_regime      = 2
            wind_speed_force_limit = v_w
            f_nF = f_out
            print()

    # Constrained tether force
    if wind_speed_regime == 2:

        mu_F = v_w / wind_speed_force_limit

        starting_point = (-0.001)
        bounds         = ((f_min, -0.001),)

        optimisation_result = op.minimize(objective_function_2, \
                                          starting_point,       \
                                          args=(mu_F, f_nF),    \
                                          bounds=bounds,        \
                                          method='SLSQP')

        # Reeling factors
        f_out = (cosine_beta_out * (mu_F - 1) + f_nF)/mu_F
        f_in  = optimisation_result['x'][0]

        # Normalized cycle power
        p_c = -objective_function_2 ([f_in], mu_F, f_nF)

        # Tether force and mechanical power during reel out
        Ft_out = q * kite_planform_area * force_factor_out * \
                     (cosine_beta_out - f_out)**2
        P_out  = Ft_out * v_w * f_out

        if P_out > nominal_generator_power:
            wind_speed_regime      = 3
            wind_speed_power_limit = v_w
            f_nP = f_out
            print()

    # Constrained tether force and generator power
    if wind_speed_regime == 3:

        mu_P  = v_w / wind_speed_power_limit
        f_out = f_nP / mu_P

        # Reduce force factor to meet tether force limit
        force_factor_out = nominal_tether_force / (q * kite_planform_area \
                           * (cosine_beta_out - f_out)**2)

        starting_point = (-0.001)
        bounds         = ((f_min, -0.001),)

        optimisation_result = op.minimize(objective_function_3, \
                                          starting_point,       \
                                          args=(mu_P, f_nP),    \
                                          bounds=bounds,        \
                                          method='SLSQP')

        # Reeling factors
        f_in  = optimisation_result['x'][0]

        # Normalized cycle power
        p_c = -objective_function_3 ([f_in], mu_P, f_nP)

        # Tether force and mechanical power during reel out
        Ft_out = q * kite_planform_area * force_factor_out * \
                     (cosine_beta_out - f_out)**2
        P_out  = Ft_out * v_w * f_out

    print("{:4.1f}".format(v_w),    \
          "{:5.3f}".format(f_out),  \
          "{:5.3f}".format(f_in),   \
          "{:5.0f}".format(Ft_out), \
          "{:6.0f}".format(P_out),  \
          "{:4.1f}".format(v_w * f_out), \
          "{:5.2f}".format(force_factor_out))

#    print("K_a    =", "{:5.3f}".format(viking_1["speedofsound"]/ref),"   ", \
#                      "{:5.3f}".format(arsia_north["speedofsound"]/ref))

    cycle_power.append(p_c * force_factor_out * kite_planform_area * Pw)
    power_ideal.append(power_factor_ideal * kite_planform_area * Pw)

power_min = np.min(power_ideal)
power_max = np.max(power_ideal)

plt.figure()
#plt.tight_layout()
plt.xlabel(r"Wind speed [m/s]")
plt.ylabel(r"Mechanical power [kW]")
plt.title('Power curve')
plt.xlim([0, 40])
plt.ylim([0, 40])
plt.vlines(wind_speed_force_limit, 0, 40, colors='k', linestyles='solid')
plt.vlines(wind_speed_power_limit, 0, 40, colors='r', linestyles='solid')
plt.plot(wind_speed, np.asarray(power_ideal)/1000, 'r', linestyle=':', label=r"$f_{\mathrm{opt}}$")
plt.plot(wind_speed, np.asarray(cycle_power)/1000, 'b', linestyle='-', label=r"$f_{\mathrm{opt}}$")
plt.savefig("powercurve.svg")