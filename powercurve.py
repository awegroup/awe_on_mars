#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 22:37:27 2023

Power curve computed following the three-regimes strategy

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
wind_speed_delta          =  1.      # m/s

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

wind_speed_range = wind_speed_max - wind_speed_min
num_wind_speeds  = int(wind_speed_range/wind_speed_delta + 1)
wind_speed       = np.linspace(wind_speed_min, wind_speed_max, num_wind_speeds)
mechanical_power = []

def objective_function_1(x):
    f_out    = x[0]
    f_in     = x[1]
    a        = 1 - 2*f_in*cosine_beta_in + f_in**2
    gamma_in = kite_lift_coefficient_in * np.sqrt(a/(1 - cosine_beta_in**2))
    return -(((cosine_beta_out - f_out)**2 - (gamma_in / force_factor_out) * \
                a) * (f_in*f_out) / (f_in - f_out))

print("num_wind_speeds = ", num_wind_speeds)

###############################################################################

# Constant properties
rho = atmosphere_density
S   = kite_planform_area
CLo = kite_lift_coefficient_out


# Loop over wind speed range
for v_w in wind_speed:

    # Dynamic pressure
    q  = 0.5 * atmosphere_density * v_w**2

    # Wind power density
    Pw = q*v_w

    # Reeling factor limits
    f_max = reeling_speed_max_limit / v_w
    f_min = reeling_speed_min_limit / v_w

    starting_point = (0.001, -0.001)
    bounds = ((0.001,  f_max),
              (f_min, -0.001),)

    optimisation_result = op.minimize(objective_function_1, starting_point,
                                      bounds=bounds, method='SLSQP')

    # Reeling factors
    f_out = optimisation_result['x'][0]
    f_in = optimisation_result['x'][1]

    Ft_out = q * kite_planform_area * force_factor_out * \
             (cosine_beta_out - f_out)**2

    print(v_w, f_out, f_in, Ft_out)

    mechanical_power.append(Pw)

mechanical_power_min = np.min(mechanical_power)
mechanical_power_max = np.max(mechanical_power)

plt.figure()
#plt.tight_layout()
plt.xlabel(r"Wind speed [m/s]")
plt.ylabel(r"Mechanical power [W]")
plt.title('Power curve')
plt.xlim([wind_speed_min, wind_speed_max])
plt.ylim([mechanical_power_min, mechanical_power_max])
plt.plot(wind_speed, mechanical_power, 'r', linestyle=':', label=r"$f_{\mathrm{opt}}$")
plt.savefig("powercurve.svg")