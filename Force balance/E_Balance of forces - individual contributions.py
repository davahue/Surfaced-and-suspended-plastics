# -*- coding: utf-8 -*-
"""
Created on Tue Aug 16, 2022

@author: Daniel Valero & Antonio Moreno-Rodenas
"""

import numpy as np
import matplotlib.pyplot as plt

# cylindrical vase geom parameters

Weight_dry_measured = 2.09E-3 # Cup_PP100 kg
Volume_total_measured = 2303.46E-9 # m3


## CUP_PP_100 definition ------------------------------------------------------
# Mass, volume
Mass = 2.09E-3 # Cup_PP_100 [kg]
Volume_total = 2303.46E-9 # [m3]

# Geometrical parameters
R = 55.2E-3 /2 # Cylinder radius [m]
H = 83E-3 # Cylinder height [m]
d = Volume_total/(np.pi * R**2 + 2*np.pi*R*H) # Thickness material [m]


rho_w = 998 # kg / m3
rho_p = Weight_dry_measured/Volume_total_measured # kg / m3

g = 9.81 # m/s^2

# bubbles
N_b = 2 # number of bubbles
d_b = 1E-3 # diameter bubbles m

# surface tension
s_t = 0.0728 # surface tension N/m
angle = 115-90 # contact angle, only the excess over 90 degrees. For instance, for contact angle of 95°, angle = 5.

print("Thickness:", d)


s = np.linspace(0, 2*R, 10000) # submergence depth [origing from lowest vase point]

#### Area circular segment
A_cs = R**2 * np.arccos(1 - s / R) - (R - s) * np.sqrt(R**2 - (R - s)**2)

#### chord circular segment
c_chord = 2*np.sqrt(2*R*s - s**2)

#### Perimeter circular segment exterior
l_cs = np.array([(R * 2 * np.arcsin(np.sqrt((R - si*0.5)*8*si)/(2*R))) if (R - si) > 0 else (2* np.pi * R - R * 2 * np.arcsin(np.sqrt((R - (2*R - si)*0.5)*8*(2*R - si))/(2*R))) for si in s ] )

### circular segment https://en.wikipedia.org/wiki/Circular_segment

print('Sanity check')
## Checking area circular segment
print('total area estimated             {}'.format(A_cs[-1]))
print('total area of circle             {}'.format(np.pi * R**2))
print('half area of circle estimated    {}'.format(R**2 * np.arccos(1 - R / R) - (R - R) * np.sqrt(R**2 - (R - R)**2)))
print('half area of circle              {}'.format(0.5 * np.pi * R**2))

## Checking ext length circular segment
print('total perimeter circle estimated {}'.format(l_cs[-1]))
print('total perimeter circle           {}'.format(2 * np.pi * R))
print('half perimeter estimated         {}'.format([(R * 2 * np.arcsin(np.sqrt((R - si*0.5)*8*si)/(2*R))) if (R - si) > 0 else (2* np.pi * R - R * 2 * np.arcsin(np.sqrt((R - (2*R - si)*0.5)*8*(2*R - si))/(2*R))) for si in [R]][0] ))
print('half perimeter circle            {}'.format(np.pi * R))

fig = plt.figure(figsize=(10,3))
ax1 = fig.add_subplot(131)
ax1.plot(s/(2*R), A_cs)
ax1.set_xlabel('Submergence degree')
ax1.set_ylabel('Area base (m2)')

ax2 = fig.add_subplot(132)
ax2.plot(s/(2*R), l_cs, label = 'circ segment wet perim')
ax2.plot(s/(2*R), c_chord, label = 'chord wet base')

ax2.set_xlabel('Submergence degree')
ax2.set_ylabel('Length m')

plt.legend()
plt.tight_layout()


V_submerged = A_cs * d + l_cs * H * d # total volume submerged vase in m3
V_total = np.pi * R**2 * d + 2*np.pi*R * H * d # total volume vase in m3

print('Check difference of volumes when fully submerged = {}'.format(V_submerged[-1] - V_total))

Volume_total_measured - V_total

fig = plt.figure(figsize=(10,3))
ax1 = fig.add_subplot(121)
ax1.plot(s/(2*R), V_submerged)

ax1.set_xlabel('submergence degree')
ax1.set_ylabel('Volume vase submerged mm3')

plt.tight_layout()

#### FORCES BALANCE

## Ft = F_bouyancy + F_bubbles + F_surfacetension - F_weight

F_bouyancy = rho_w * g * V_submerged
F_weight = rho_p * g * V_total
F_bubbles = rho_w * g * N_b * 4/3 * np.pi * (d_b/2)**3

### Proper surface tension (geometry dependent, wetting angle given by material properties, PP-water contact angle)

F_surface_bottom = s_t * np.sin(np.deg2rad(angle)) * 2*c_chord

F_surface_sidesExt = [s_t * 2 * H * np.sin(np.arcsin((si - R)/R) + np.deg2rad(angle)) if si > R else s_t * 2 * H * - np.cos(np.arccos((R - si)/R) + np.deg2rad(angle))  for si in s] ### upper and lower contribution s > R
F_surface_sidesInt = [ s_t * 2 * H * -np.cos(np.deg2rad(angle) + np.arccos((si - R)/R)) if si > R else  s_t * 2 * H * np.sin(np.arcsin((R - si)/R) + np.deg2rad(angle)) for si in s]


F_surfacetension = F_surface_bottom + F_surface_sidesExt + F_surface_sidesInt


plt.figure()

plt.plot(F_surface_sidesExt)
plt.plot(F_surface_sidesInt)
plt.plot(np.array(F_surface_sidesInt) + np.array(F_surface_sidesExt))

s[np.argmin(np.abs(F_bouyancy + F_bubbles + F_surfacetension - F_weight))]/(2*R) # submergence level equilibrium


plt.figure(figsize = (10, 5))

plt.plot(s/(2*R), F_bouyancy, '-', color = 'k', label = 'Bouyancy')
plt.plot(s/(2*R), F_bouyancy + F_bubbles, '--', color = 'C0', label = 'Bouyancy + bubbles')
plt.plot(s/(2*R), F_bouyancy + F_bubbles + F_surfacetension, color = 'C0', label = 'Bouyancy + bubbles + surf tension')

plt.axvline(s[np.argmin(np.abs(F_bouyancy + F_bubbles + F_surfacetension - F_weight))]/(2*R), color = 'r', label = 'Submergence equilibrium')
plt.axhline(F_weight, color = 'C1', label = 'Weight')

#plt.xlim(0,1)
plt.ylabel('F (N)')
plt.legend()

plt.figure(figsize = (10, 5))

plt.plot(s/(2*R), F_bouyancy, '-', color = 'k', label = 'Bouyancy')
plt.plot(s/(2*R), F_bouyancy + F_bubbles, '--', color = 'C0', label = 'Bouyancy + bubbles')
plt.plot(s/(2*R), F_bouyancy + F_bubbles + F_surfacetension, color = 'C0', label = 'Bouyancy + bubbles + surf tension')

plt.axvline(s[np.argmin(np.abs(F_bouyancy + F_bubbles + F_surfacetension - F_weight))]/(2*R), color = 'r', label = 'Submergence equilibrium')
plt.axhline(F_weight, color = 'C1', label = 'Weight')

#plt.xlim(0,1)
plt.ylabel('F (N)')
plt.legend()


# Figure ----------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(4,3.5))


Bou_C = F_bouyancy/(F_bouyancy + F_bubbles + F_surfacetension + F_weight)
St_C = F_surfacetension/(F_bouyancy + F_bubbles + F_surfacetension + F_weight)
Bub_C = F_bubbles/(F_bouyancy + F_bubbles + F_surfacetension + F_weight)
W_C = F_weight/(F_bouyancy + F_bubbles + F_surfacetension + F_weight)

col = ["firebrick", "forestgreen", "k", "lightskyblue"] # "goldenrod"
lbls = ["$F_{\\rho}$ (weight)", "$F_{\\rho}$ (buoyancy)", "$F_{b}$", "$F_{\\sigma}$"]

plt.fill_between(s/(2*R), y1 = np.zeros(len(s)), y2 = W_C, color = col[0], 
                 label = lbls[0])
plt.fill_between(s/(2*R), y1 = W_C, y2 = Bou_C + W_C, color = col[1],
                 label = lbls[1])
plt.fill_between(s/(2*R), y1 = Bou_C + W_C, y2 = Bou_C + W_C + Bub_C, 
                 color = col[2], label = lbls[2])
plt.fill_between(s/(2*R), y1 = Bou_C + W_C + Bub_C, y2 = Bou_C + W_C + Bub_C + St_C,
                 color = col[3], label = lbls[3])

plt.vlines(s[np.argmin(np.abs(F_bouyancy + F_bubbles + F_surfacetension - F_weight))]/(2*R), 0, 1, 'k', ls="--")
plt.hlines(0.5, 0, 1, 'k', linestyle = '--')

plt.xlim(0,1)
plt.ylim(0,1)

# plt.ylabel(r'$\frac{|F_i|}{\sum |F_i|}$ (-)')
plt.ylabel('$\mid F_i \mid / \sum \mid F_i \mid$ (-)')
plt.xlabel('Submergence, s/(2R) (-)')

plt.legend(loc = 3)

plt.tight_layout()

plt.savefig('E_contact115_nb2_db_1.pdf', dpi=600, format='pdf',  bbox_inches='tight')
plt.savefig('E_contact115_nb2_db_1.png', dpi=600, format='png',  bbox_inches='tight')
plt.savefig('E_contact115_nb2_db_1.svg', dpi=600, format='svg',  bbox_inches='tight')



# # FIGURE ----------------------------------------------------------------------

# fig, ax = plt.subplots(figsize=(4,3.5))


# Bou_C = F_bouyancy/(F_bouyancy + F_bubbles  + F_weight)
# St_C = F_surfacetension/(F_bouyancy + F_bubbles  + F_weight)
# Bub_C = F_bubbles/(F_bouyancy + F_bubbles  + F_weight)
# W_C = F_weight/(F_bouyancy + F_bubbles  + F_weight)

# col = ["firebrick", "forestgreen", "k", "lightskyblue"] # "goldenrod"
# lbls = ["$F_{\\rho}$ (weight)", "$F_{\\rho}$ (buoyancy)", "$F_{b}$", "$F_{\\sigma}$"]

# plt.fill_between(s/(2*R), y1 = np.zeros(len(s)), y2 = W_C, color = col[0], alpha = 0.8, label = 'Weight')
# plt.fill_between(s/(2*R), y1 = W_C, y2 = Bou_C + W_C, color = col[1], alpha = 0.8, label = 'Bouyancy')
# plt.fill_between(s/(2*R), y1 = Bou_C + W_C, y2 = Bou_C + W_C + Bub_C, color = 'k', alpha = 0.8, label = 'Bubbles')

# plt.vlines(s[np.argmin(np.abs(F_bouyancy + F_bubbles - F_weight))]/(2*R), 0, 1, 'k', ls="--")
# plt.hlines(0.5, 0, 1, 'k', linestyle = '--')

# plt.xlim(0,1)
# plt.ylim(0,1)

# plt.ylabel(r'$\frac{|F_i|}{\sum |F_i|}$ (-)')
# plt.xlabel('Submergence, s/(2R) (-)')

# plt.legend(loc = 3)

# plt.tight_layout()