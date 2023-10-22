#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 20:50:10 2023

Two locations on Mars are investigated and compared to sea-level standard (SLS)
conditions on Earth:

1. Viking 1    at  22.27° N, 312.05° E and z = -3,600 m (altitude above aeroid)
2. Arsia North at -3.062° N, 236.07° E and z =  4,550 m (altitude above aeroid)

Data is extracted from the Mars Climate Database (MCD) v6.1
https://www-mars.lmd.jussieu.fr/mcd_python/

Since we are interested in local average data we set time coordinates
Ls=all and Lt=all, and set averaging=diurnal. This will generate 2d-data files
with identical columns. By extracting e.g. the first column we have the
variation of the diurnal (daily) average of the solar longitudes.
Each column of the 2d-data files contains 25 data points, for solar longitudes
0 to 360. We use data at 10 m above ground level.
The results of the MCD online requests for the two locations are collected in
a separate data file, where also the averaging over the solar longitudes is
performed. The resulting averages are used as input in the Python script below.

Wikipedia: Sea-level standard (SLS) conditions.
https://en.wikipedia.org/wiki/Standard_sea-level_conditions

https://www.eea.europa.eu/publications/europes-changing-climate-hazards-1/wind/wind-mean-wind-speed

@author: Roland Schmehl
"""
import numpy as np
import sys

# Sea-level standard (SLS) conditions from https://en.wikipedia.org/wiki/Standard_sea-level_conditions
earth_sls =	{
  "name":          "Earth SLS",
  "density":        1.225,     # kg/m^3 @SLS
  "viscosity":      1.8e-5,    # Pa*s   @SLS
  "gravity":        9.8,       # m/s^2  @SLS
  "speedofsound":   343,       # m/s    @SLS and T = 293.15 K (= 20 °C)
  "windspeed_low":  3,         # m/s    @10 m annual mean for European land area
  "windspeed_high": 10,        # m/s    @100 m annual mean for North Sea region
}

viking_1 =	{
  "name":          "Viking_1",
  "density":        0.020,     # kg/m^3 Williams (2020)
  "viscosity":      1.09e-5,   # Pa*s
  "gravity":        3.7,       # m/s^2
  "speedofsound":   233,       # m/s
  "windspeed_low":  5,         # m/s    @10 m
  "windspeed_high": 7,         # m/s    @100 m
  "k_vw":           np.array([2, 3]),
}

arsia_north = {
  "name":          "Arsia_North",
  "density":        0.010,     # kg/m^3 MCD
  "viscosity":      1.04e-5,   # Pa*s
  "gravity":        3.7,       # m/s^2
  "speedofsound":   228,       # m/s
  "windspeed_low":  15,        # m/s    @10 m
  "windspeed_high": 18,        # m/s    @100 m
  "k_vw":           np.array([2, 3]),
}

print("Table 6.2: Scaling relations...")
print()
print("         Viking_1  Arsia_North")
ref = earth_sls["density"]
print("K_rho  =", "{:6.4f}".format(viking_1["density"]/ref),"  ", \
                  "{:6.4f}".format(arsia_north["density"]/ref))
ref = earth_sls["viscosity"]
print("K_mu   =", "{:5.3f}".format(viking_1["viscosity"]/ref),"   ", \
                  "{:5.3f}".format(arsia_north["viscosity"]/ref))
ref = earth_sls["gravity"]
print("K_g    =", "{:5.3f}".format(viking_1["gravity"]/ref),"   ", \
                  "{:5.3f}".format(arsia_north["gravity"]/ref))
ref = earth_sls["speedofsound"]
print("K_a    =", "{:5.3f}".format(viking_1["speedofsound"]/ref),"   ", \
                  "{:5.3f}".format(arsia_north["speedofsound"]/ref))
ref = earth_sls["windspeed_low"]
print("K_vw   =", "{:5.3f}".format(viking_1["windspeed_low"]/ref),"   ", \
                  "{:5.3f}".format(arsia_north["windspeed_low"]/ref), "    (@10 m & low vw@Earth)")
ref = earth_sls["windspeed_high"]
print("K_vw   =", "{:5.3f}".format(viking_1["windspeed_high"]/ref),"   ", \
                  "{:5.3f}".format(arsia_north["windspeed_high"]/ref), "    (@100 m & high vw@Earth)")

for site in [viking_1, arsia_north]:

    name  = site["name"]
    k_rho = site["density"]/earth_sls["density"]
    k_mu  = site["viscosity"]/earth_sls["viscosity"]
    k_g   = site["gravity"]/earth_sls["gravity"]
    k_a   = site["speedofsound"]/earth_sls["speedofsound"]
    k_vw  = site["k_vw"]

    print()
    print("Table 6.3:", name)
    print("K_rho  =", "{:6.4f}".format(k_rho))

    # Scaling factor wind speed (preset)
    sys.stdout.write("K_vw   = ")
    for k in k_vw:
        sys.stdout.write("{:<10.0f}".format(k))
    print("(preset, i.e. not from above)")

    # Scaling factor planform area
    k_s = 1/(k_rho*k_vw**3)
    sys.stdout.write("K_S    = ")
    for k in k_s:
        sys.stdout.write("{:<10.1f}".format(k))
    print()

    # Scaling factor planform span
    k_b = 1/np.sqrt(k_rho*k_vw**3)
    sys.stdout.write("K_b    = ")
    for k in k_b:
        sys.stdout.write("{:<10.2f}".format(k))
    print()

    # Scaling factor tether force
    k_f = 1/k_vw
    sys.stdout.write("K_F    = ")
    for k in k_f:
        sys.stdout.write("{:<10.3f}".format(k))
    print()

    # Scaling factor tether diameter
    k_d = 1/np.sqrt(k_vw)
    sys.stdout.write("K_d    = ")
    for k in k_d:
        sys.stdout.write("{:<10.3f}".format(k))
    print()

    # Scaling factor membrane thickness
    k_t = np.sqrt(k_rho*k_vw)
    sys.stdout.write("K_t    = ")
    for k in k_t:
        sys.stdout.write("{:<10.3f}".format(k))
    print()

    # Scaling factor kite mass
    k_m = 1/np.sqrt(k_rho*k_vw**5)
    sys.stdout.write("K_m    = ")
    for k in k_m:
        sys.stdout.write("{:<10.3f}".format(k))
    print()

    # Scaling factor gravitational force
    k_fg = np.sqrt(k_g**2/(k_rho*k_vw**5))
    sys.stdout.write("K_Fg   = ")
    for k in k_fg:
        sys.stdout.write("{:<10.3f}".format(k))
    print()

    # Scaling factor launching easiness
    k_nu = np.sqrt(k_rho*k_vw**3/k_g**2)
    sys.stdout.write("K_nu   = ")
    for k in k_nu:
        sys.stdout.write("{:<10.3f}".format(k))
    print()

    # Scaling factor gravitation material strain
    k_sigg = np.sqrt(k_g**2/(k_rho*k_vw**3))
    sys.stdout.write("K_sigg = ")
    for k in k_sigg:
        sys.stdout.write("{:<10.3f}".format(k))
    print()

    # Scaling factor Mach number
    k_ma = k_vw/k_a
    sys.stdout.write("K_Ma   = ")
    for k in k_ma:
        sys.stdout.write("{:<10.3f}".format(k))
    print()

    # Scaling factor Reynolds number
    k_re = np.sqrt(k_rho/(k_vw*k_mu**2))
    sys.stdout.write("K_Re   = ")
    for k in k_re:
        sys.stdout.write("{:<10.3f}".format(k))
    print()

    # Scaling factor turning performance
    k_turn = k_vw**2
    sys.stdout.write("K_turn = ")
    for k in k_turn:
        sys.stdout.write("{:<10.0f}".format(k))
    print()

    # Turning gravitational importance
    eps_g = k_sigg
    sys.stdout.write("eps_g  = ")
    for e in eps_g:
        sys.stdout.write("{:<10.3f}".format(e))
    print()

