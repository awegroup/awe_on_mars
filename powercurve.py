#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 22:37:27 2023

Power curve computed following the three-regimes strategy.

In this implementation, a constant lift-to-drag ratio E is prescribed during
reel-in, which means that the elevation angle during reel-in is computed as a
dependent variable.

@author: Roland Schmehl
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import optimize as op
from pylab import np
mpl.rcParams['font.family'] = "Open Sans"
mpl.rcParams.update({'font.size': 18})
mpl.rcParams['figure.figsize'] = 10, 5.625
mpl.rc('xtick', labelsize=16)
mpl.rc('ytick', labelsize=16)
mpl.rcParams['pdf.fonttype'] = 42 # Output Type 3 (Type3) or Type 42 (TrueType)

# Environmental properties
atmosphere_density        =  0.01    # kg/**3
wind_speed_min            =  1.      # m/s
wind_speed_max            =  45.     # m/s
wind_speed_delta          =  0.1     # m/s

# Kite properties
kite_planform_area       =  200.     # m**2
kite_lift_coefficient_out =  0.71    # -
kite_drag_coefficient_out =  0.14    # -
kite_lift_coefficient_in  =  0.39    # -
kite_drag_coefficient_in  =  0.39    # -

# Tether properties
nominal_tether_force      =  5100.   # N
tether_drag_coefficient   =  1.1     # -
tether_diameter           =  0.00484 # m
tether_length_max         =  385.    # m
tether_length_min         =  240.    # m

# Generator properties
nominal_generator_power   =  77000.  # W

# Operational parameters
elevation_angle_out       =  25.     # deg
reeling_speed_min_limit   = -21.     # m/s
reeling_speed_max_limit   =   8.     # m/s

# Derived properties
rm_out = 0.5 * (tether_length_min + tether_length_max)
CD_out = kite_drag_coefficient_out + 0.25 * tether_drag_coefficient \
         * tether_diameter * rm_out / kite_planform_area
E2out  = (kite_lift_coefficient_out / CD_out)**2
E2in   = (kite_lift_coefficient_in  / kite_drag_coefficient_in)**2
cosine_beta_out  = np.cos(np.radians(elevation_angle_out))
force_factor_out = kite_lift_coefficient_out * np.sqrt(1+1/E2out) * (1+E2out)
force_factor_in  = kite_lift_coefficient_in  * np.sqrt(1+1/E2in)
power_factor_ideal = force_factor_out * cosine_beta_out**3 * 4/27

wind_speed_range = wind_speed_max - wind_speed_min
num_wind_speeds  = int(wind_speed_range/wind_speed_delta + 1)
wind_speed  = np.linspace(wind_speed_min, wind_speed_max, num_wind_speeds)

# Lists
reeling_factor_out = []
reeling_factor_in  = []
tether_force_out   = []
tether_force_in    = []
power_out          = []
power_in           = []
cycle_power        = []
power_ideal        = []
elevation_angle_in = []


# Objective function for the three wind speed domains
def objective_function_1(x):
    f_out    = x[0]
    f_in     = x[1]
    return -((cosine_beta_out - f_out)**2 \
             - (force_factor_in / force_factor_out) \
               * (np.sqrt(1 + E2in*(1 - f_in**2)) - f_in)**2/(1 + E2in)) \
            * (f_in*f_out) / (f_in - f_out)

def objective_function_2(x, mu_F, f_nF):
    f_in     = x[0]
    b        = (mu_F - 1) * cosine_beta_out + f_nF
    return -(((cosine_beta_out - f_nF) / mu_F)**2  \
             - (force_factor_in / force_factor_out) \
               * (np.sqrt(1 + E2in*(1 - f_in**2)) - f_in)**2/(1 + E2in))  \
            * f_in*b/(mu_F*f_in-b)

def objective_function_3(x, mu_P, f_nP):
    f_in     = x[0]
#    return -(((cosine_beta_out - f_nP) / mu_P)**2  \
#             - (force_factor_in / force_factor_out) \
#               * (np.sqrt(1 + E2in*(1 - f_in**2)) - f_in)**2/(1 + E2in))  \
#            * f_in*f_nP/(mu_P*f_in-f_nP)
    return -((cosine_beta_out - f_nP/mu_P)**2  \
             - (force_factor_in / force_factor_out) \
               * (np.sqrt(1 + E2in*(1 - f_in**2)) - f_in)**2/(1 + E2in))  \
            * f_in*f_nP/(mu_P*f_in-f_nP)

print("num_wind_speeds = ", num_wind_speeds)

###############################################################################

# Initialize wind speed regimes
wind_speed_regime      = 1
wind_speed_force_limit = 0
wind_speed_power_limit = 0
print("Wind speed regime 1")

# Loop over wind speed range
for v_w in wind_speed:

    # Dynamic pressure
    q  = 0.5 * atmosphere_density * v_w**2

    # Wind power density
    P_w = q*v_w

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
                                          method='COBYLA')

        # Reeling factors
        f_out = optimisation_result['x'][0]
        f_in  = optimisation_result['x'][1]

        # Normalized cycle power
        p_c = -objective_function_1 ([f_out, f_in])

        # Tether force during reel-out
        Ft_out = q * kite_planform_area * force_factor_out  \
                   * (cosine_beta_out - f_out)**2

        if Ft_out > nominal_tether_force:
            wind_speed_regime = 2

            # Determine precise value of v_w,F by interval bisection
            v_b  = v_w
            v_a  = v_w - wind_speed_delta
            c    = 0.5 * atmosphere_density * kite_planform_area \
                       * force_factor_out * (cosine_beta_out - f_out)**2
            nmax = 100
            eps  = 0.1
            for i in range(nmax):
                v  = (v_a + v_b)/2
                Ft = c * v**2
                if Ft > nominal_tether_force:
                    v_b = v
                else:
                    v_a = v
                if abs(Ft-nominal_tether_force) < eps:
                    break
            else:
                print("!!! search v_w,F stopped after nmax=", nmax, "iterations")
                print("--> increase nmax and rerun")

            wind_speed_force_limit = v
            f_nF  = f_out # works because f_out is constant in regime 1

            print()
            print("Wind speed regime 2 with v_n,F at", "{:5.2f}".format(wind_speed_force_limit))
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
                                          method='COBYLA')

        # Reeling factors
        f_out = (cosine_beta_out * (mu_F - 1) + f_nF)/mu_F
        f_in  = optimisation_result['x'][0]

        # Normalized cycle power
        p_c = -objective_function_2 ([f_in], mu_F, f_nF)

        # Tether force and mechanical power during reel out
        Ft_out = q * kite_planform_area * force_factor_out * \
                     (cosine_beta_out - f_out)**2

        # Mechanical power during reel out
        P_out  = Ft_out * v_w * f_out

        if P_out > nominal_generator_power:
            wind_speed_regime = 3

            # Determine precise value of v_w,P by interval bisection
            v_b  = v_w
            v_a  = v_w - wind_speed_delta
            c    = 0.5 * atmosphere_density * kite_planform_area \
                       * force_factor_out
            nmax = 100
            eps  = 1
            for i in range(nmax):
                v  = (v_a + v_b)/2
                mu = v / wind_speed_force_limit
                f  = (cosine_beta_out * (mu - 1) + f_nF)/mu
                P  = c * (cosine_beta_out - f)**2 * v**3 * f
                if P > nominal_generator_power:
                    v_b = v
                else:
                    v_a = v
                if abs(P-nominal_generator_power) < eps:
                    break
            else:
                print("!!! search v_w,P stopped after nmax=", nmax, "iterations")
                print("--> increase nmax and rerun")

            wind_speed_power_limit = v
            f_nP = f

            print()
            print("Wind speed regime 3 with v_n,P at", "{:5.2f}".format(wind_speed_power_limit))
            print()

    # Constrained tether force and generator power
    if wind_speed_regime == 3:

        mu_P  = v_w / wind_speed_power_limit
        f_out = f_nP / mu_P

        # Reduce force factor to comply with tether force limit
        force_factor_out = nominal_tether_force / (q * kite_planform_area \
                           * (cosine_beta_out - f_out)**2)

        # Alternative strategy to depower: increasing the elevation angle
#        cosine_beta_out = np.sqrt(nominal_tether_force / (q \
#                          * kite_planform_area * force_factor_out)) + f_out


        starting_point = (-0.001)
        bounds         = ((f_min, -0.001),)

        optimisation_result = op.minimize(objective_function_3, \
                                          starting_point,       \
                                          args=(mu_P, f_nP),    \
                                          bounds=bounds,        \
                                          method='COBYLA')

        # Reeling factors
        f_in  = optimisation_result['x'][0]

        # Normalized cycle power
        p_c = -objective_function_3 ([f_in], mu_P, f_nP)

    # Tether force
    Ft_out = q * kite_planform_area * force_factor_out  \
               * (cosine_beta_out - f_out)**2
    Ft_in  = q * kite_planform_area * force_factor_in  \
               * (np.sqrt(1 + E2in*(1 - f_in**2)) - f_in)**2/(1 + E2in)

    # Mechanical power
    P_out  = Ft_out * v_w * f_out
    P_in   = Ft_in  * v_w * f_in

    # Elevation angle reel-in-phase
    beta_in = np.arccos((np.sqrt(1 + E2in*(1 - f_in**2)) \
                                   + f_in*E2in)/(1 + E2in))

    print("{:4.1f}".format(v_w),    \
          "{:5.3f}".format(f_out),  \
          "{:5.3f}".format(f_in),   \
          "{:5.0f}".format(Ft_out), \
          "{:5.0f}".format(Ft_in),  \
          "{:6.0f}".format(P_out),  \
          "{:6.0f}".format(P_in),   \
          "{:4.1f}".format(v_w * f_out), \
          "{:4.1f}".format(v_w * f_in), \
          "{:5.2f}".format(force_factor_out), \
          "{:5.2f}".format(force_factor_in), \
          "{:4.1f}".format(np.degrees(beta_in)))

    reeling_factor_out.append(f_out)
    reeling_factor_in.append(f_in)
    tether_force_out.append(Ft_out)
    tether_force_in.append(Ft_in)
    power_out.append(P_out)
    power_in.append(P_in)
    cycle_power.append(p_c * force_factor_out * kite_planform_area * P_w)
    power_ideal.append(power_factor_ideal * kite_planform_area * P_w)
    elevation_angle_in.append(np.degrees(beta_in))

print()
cp = np.array(cycle_power)
print("> Rated power: ", "{:6.0f}".format(max(cp)), "W at ", \
      "{:4.1f}".format(wind_speed[cp.argmax()]), "m/s")

power_min = np.min(power_ideal)
power_max = np.max(power_ideal)

fig, ax1 = plt.subplots()
ax1.set(xlabel=r"Wind speed, m/s", ylabel=r"Mechanical power, kW")
ax1.set_xlim([0, 50])
ax1.set_ylim([0, 80])
#ax1.grid()
ax1.vlines(wind_speed_force_limit, 0, 100, colors='k', linestyles=':')
ax1.vlines(wind_speed_power_limit, 0, 100, colors='k', linestyles=':')
ax1.annotate("1",(15,40), ha="center", va="center", bbox={"boxstyle" : "circle", "color":"white", "ec" : "k"})
ax1.annotate("2",(30,40), ha="center", va="center", bbox={"boxstyle" : "circle", "color":"white", "ec" : "k"})
ax1.annotate("3",(37,40), ha="center", va="center", bbox={"boxstyle" : "circle", "color":"white", "ec" : "k"})
ax1.annotate(r"$v_{\mathrm{n,F}}$",(24.5,-2.5), annotation_clip=False, ha="center", va="center")
ax1.annotate(r"$v_{\mathrm{n,P}}$",(35.,-2.5), annotation_clip=False, ha="center", va="center")
ax1.plot(wind_speed,  np.asarray(power_ideal)/1000, 'k', linestyle=':', label=r"$P_{\mathrm{opt}}$")
ax1.plot(wind_speed,  np.asarray(cycle_power)/1000, 'b', linestyle='-', label=r"$P_{\mathrm{c}}$")
ax1.plot(wind_speed,  np.asarray(power_out)/1000, 'g', linestyle='--', label=r"$P_{\mathrm{o}}$")
ax1.plot(wind_speed, -np.asarray(power_in)/1000, 'r', linestyle='--', label=r"$-P_{\mathrm{i}}$")
ax1.legend(facecolor="white", edgecolor="white")
fig.savefig("powercurve.svg")

fig, ax1 = plt.subplots()
ax1.set(xlabel=r"Wind speed, m/s", ylabel=r"Reeling factor")
ax1.set_xlim([0, 50])
ax1.set_ylim([0, 1.5])
#ax1.grid()
ax1.vlines(wind_speed_force_limit, 0, 100, colors='k', linestyles=':')
ax1.vlines(wind_speed_power_limit, 0, 100, colors='k', linestyles=':')
ax1.annotate(r"$v_{\mathrm{n,F}}$",(24.5,-0.045), annotation_clip=False, ha="center", va="center")
ax1.annotate(r"$v_{\mathrm{n,P}}$",(35,-0.045), annotation_clip=False, ha="center", va="center")
ax1.plot(wind_speed,  np.asarray(reeling_factor_out), 'g', linestyle='--', label=r"$f_{\mathrm{o}}$")
ax1.plot(wind_speed, -np.asarray(reeling_factor_in), 'r', linestyle='--', label=r"$-f_{\mathrm{i}}$")
ax2 = ax1.twinx()
ax2.set(ylabel=r"Elevation angle, deg")
ax2.set_ylim([0, 140])
ax2.plot(wind_speed,  np.asarray(elevation_angle_in), 'b', linestyle='-', label=r"$\beta_{\mathrm{i}}$")
fig.legend(facecolor="white", edgecolor="white", loc="upper right", bbox_to_anchor=(1,1), bbox_transform=ax1.transAxes)
fig.savefig("operations.svg")
