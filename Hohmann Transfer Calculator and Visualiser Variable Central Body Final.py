# This calculator is for calculating the delta-V's required for a Hohmann Transfer between
# two coplanar circular orbits, as well as plotting the trajectory. Also only accurate when
# orbit radii is within the selected central body's Sphere of Influence (SOI) with respect to the sun,
# Calculation for SOI is also shown below. Use orbit radii only within this SOI for accurate results.
# Method inspired by Hale
# Written by John Matthew Atok

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.patches import Arc

# Celestial Body Values, all taken from NSSDC (NASA) body fact sheets.
# mu = Gravitational parameter, GM (km^3/s^2)
# r = Volumetric mean radius (km)
# m = Mass (kg)
# a =  Semi-major axis (km)

celestial_bodies = {
    
    "Sun": {
        "mu": 132712e6,
        "r": 695700,
        "m": 1988400e24,
        "colour": 'yellow'

    },

    "Mercury": {
        "mu": 0.022032e6,
        "r": 2439.7,
        "m": 0.33010e24,
        "a": 57.909e6,
        "colour": 'dimgrey'
    },

    "Venus":{
        "mu": 0.32486e6,
        "r": 6051.8,
        "m": 4.8673e24,
        "a": 108.210e6,
        "colour": 'lightsalmon'
    },


    "Earth": {
        "mu": 0.39860e6,
        "r": 6371.000,
        "m": 5.9722e24,
        "a": 149.598e6,
        "colour": 'blue'	
    },

    "Moon":{
        "mu": 0.00490e6,
        "r": 1737.4,
        "m": 0.07346e24,
        "a": 0.3844e6,                  # Semi-major axis orbiting Earth (km)
        "colour": 'silver'
    },

    "Mars":{
        "mu": 0.042828e6,
        "r": 3389.5,
        "m": 0.64169e24,
        "a": 227.956e6,
        "colour": 'darkred'
    },

    "Jupiter":{
        "mu": 126.687e6,
        "r": 69911,
        "m": 1898.13e24,
        "a": 778.479e6,
        "colour": 'chocolate'

    },

    "Saturn":{
        "mu": 37.931e6,
        "r": 58232,
        "m": 568.32e24,
        "a": 1432.041e6,
        "colour": 'darkorange'
    },

    "Uranus":{
        "mu": 5.7940e6,
        "r": 25362,
        "m": 86.811e24,
        "a": 2867.043e6,
        "colour": 'paleturquoise'
    },

    "Neptune":{
        "mu": 6.8351e6,
        "r": 24622,
        "m": 102.409e24,
        "a": 4514.953e6,
        "colour": 'royalblue'
    },

    "Pluto":{
        "mu": 0.000870e6,	
        "r": 1188,
        "m": 0.01303e24,
        "a": 5869.656e6,
        "colour": 'lightslategrey'
    }
}

# Inputs

central_body = "Earth"                                                 # Selects the central body from the "celestial_bodies" above
alt1 = 3000                                                             # Altitude of initial orbit (km)
alt2 = 10000                                                 # Altitude of final orbit (km)

# Load constants for selected central body

constants = celestial_bodies[central_body]                         # Pulls defined constants from "celestial_bodies"
mu = constants['mu']                                                   # Sets mu
r = constants['r']                                                     # Sets r
m = constants['m']                                                     # Sets m
a = constants['a']                                                     # Sets a

# SOI with respect to central body

if central_body == "Moon":

    SOI = a * (m/celestial_bodies["Earth"]["m"])**0.4

    print(f"SOI of {central_body} with respect to the Earth is {SOI} km")                         # SOI of Moon wrt the Earth (km)

else:

    SOI = a * (m/celestial_bodies["Sun"]["m"])**0.4                                               # SOI wrt the Sun (km)

    print(f"SOI of {central_body} with respect to the Sun is {SOI} km")

# Initial Orbit

r1 = r + alt1                                                         # Radius of initial orbit (km)
v1 = np.sqrt(mu/r1)                                                   # Velocity at initial orbit (km/s)
E1 = -mu / (2 * r1)                                                   # Energy (km^2/s^2)

print(f"Velocity at INITIAL orbit is {v1} km/s")
print(f"Energy at INITIAL orbit is {E1} km^2/s^2")

# Final Orbit

r2 = r + alt2                                                         # Radius of initial orbit (km)
v2 = np.sqrt(mu/r2)                                                   # Velocity at initial orbit (km/s)
E2 = -mu / (2 * r2)                                                   # Energy (km^2/s^2)

print(f"Velocity at FINAL orbit is {v2} km/s")
print(f"Energy at FINAL orbit is {E2} km^2/s^2")

# Transfer Trajectory

a_t = (r1 + r2)/2                                                       # Semi-major axis of transfer ellipse (km)
E_t = -mu / (2 * a_t)                                                 # Energy (km^2/s^2)
v_tp = np.sqrt(2 * (E_t + (mu/r1)))                                   # Velocity at periapsis of transfer ellipse (km/s)
H_t = r1 * v_tp                                                         # Specific angular momentum (km^2/s)
e_t = np.sqrt(1 + ((2 * E_t * H_t**2)/mu**2))                         # Eccentricity of transfer ellipse
v_ta = H_t/r2                                                           # Velocity at apoapsis of transfer ellipse (km/s)

print(f"Semi-major axis of TRANSFER ellipse is {a_t} km")
print(f"Energy of TRANSFER ellipse is {E_t} km^2/s^2")
print(f"Velocity at perigee of transfer orbit is {v_tp} km/s")
print(f"Velocity at apogee of transfer orbit is {v_ta} km/s")
print(f"Specific angular momentum at perigee of transfer orbit is {H_t} km^2/s")
print(f"Eccentricity of transfer orbit is {e_t}")

# Delta-V Calculations

dV1 = v_tp - v1                                                         # Delta-V required for transfer ellipse insertion
dV2 = v2 - v_ta                                                         # Delta-V required for final orbit insertion
dV_T = abs(dV1) + abs(dV2)                                              # Total Delta V required for Hohmann (magnitude of values, ignore signs)

print(f"Transfer ellipse insertion (dV1) delta-V is {dV1} km/s")
print(f"Final orbit insertion (dV2) delta-V is {dV2} km/s")
print(f"Total delta-V is {dV_T} km/s")

# Time of Flight (TOF) Calculations

P = (2 * np.pi) * np.sqrt((a_t**3)/mu)                                # Full period of transfer ellipse
TOF_s = P/2                                                             # Half period, transfer time (seconds)
TOF_m = TOF_s / 60                                                      # Transfer time (minutes)
TOF_h = TOF_m / 60                                                      # Transfer time (hours)
TOF_d = TOF_h / 24                                                      # Transfer time (days)
TOF_y = TOF_d / 360

print(f"Transfer time is {TOF_s} seconds")
print(f"Transfer time is {TOF_m} minutes")
print(f"Transfer time is {TOF_h} hours")
print(f"Transfer time is {TOF_d} days")
print(f"Transfer time is {TOF_y} years")

# Axes 

fig, ax = plt.subplots(figsize=(10,10))
ax.set_aspect('equal')
ax.set_facecolor('black')

# Plotting Central Body

colour = constants["colour"]

cb_plot = plt.Circle((0,0), r, color = colour, label = central_body)
ax.add_patch(cb_plot)

# Plot Initial Orbit

orb1 = plt.Circle((0,0), r1, color = 'green', 
                  fill = False, 
                  label = 'Initial Orbit')

ax.add_patch(orb1)

# Plot Final Orbit

orb2 = plt.Circle((0,0), r2, color = 'yellow',
                   fill = False, 
                   label = 'Final Orbit')

ax.add_patch(orb2)

# Plot Transfer Ellipse
                                     
b_t = np.sqrt(r1 * r2)                                                      # Semi-minor axis of transfer ellipse (km)
ell_centre_x = (r2 - r1) / 2
tran_ell = Arc(xy=(-ell_centre_x, 0), 
                   width = 2 * a_t, 
                   height = 2 * b_t, 
                   edgecolor = 'purple', 
                   facecolor = 'none',
                   linestyle = '--',
                   theta1 = 0,
                   theta2 = 180, 
                   linewidth = 2,
                   label = 'Hohmann Transfer')

ax.add_patch(tran_ell)

# Plot graph

if alt2 > alt1:

    plot_lim = r2  * 1.1

else:

    plot_lim = r1  * 1.1

ax.set_xlim(-plot_lim, plot_lim)
ax.set_ylim(-plot_lim, plot_lim)

ax.set_title('Hohmann Transfer Visual', color = 'white')
ax.set_xlabel('X (km)')
ax.set_ylabel('Y (km)')
ax.legend()
ax.grid(True, linestyle = ':', color = 'gray', alpha = 0.5)

if r2 > SOI:
    print(f"WARNING: Final orbital radius (r2 = {r2} km) is greater than the calculated SOI ({SOI} km)")

plt.show()
          
# 