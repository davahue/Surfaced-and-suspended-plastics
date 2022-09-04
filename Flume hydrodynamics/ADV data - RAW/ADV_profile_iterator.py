# -*- coding: utf-8 -*-
"""
Created on Sun May 16, 2021

@author: Daniel Valero (d.valero@un-ihe.org)
---------------------------
This script iterates over the ADV folders and compiles the data.

"""
# Import libraries ----------------
import numpy as np
from scipy import stats
#import pandas as pd
from matplotlib import pyplot as plt
import glob

# Import in-house libraries -------
import faux_input as fi

# Filtering:
SNRthreshold = 10
CORthreshold = 80
# folders def
folders = ['Upstream', 'Middle', 'Downstream',  'Downstream-fine']
x_xsec = [0, 3.4, 6.8] # m, 0 is where the rough domain starts. -------------- # RECHECK WITH BIRUK
valid_fmt = '*.vna'
info_file = 'info.txt'
PI = 0.20; PI_lbl = '$\Pi$ = ' + str(PI); # this is compatible with measurements of
# Nezu and Rodi 1986, and Castro-Orgaz 2009 in free surface boundary layers (developed and developing).

# Walk over the directories
for fol in folders:
    subfolders = glob.glob(fol+'/*/')
    for subfol in subfolders:        
        
        if not glob.glob(subfol+info_file):
            print('No info file in folder: ', subfol) # if there is no info.txt file
                                                      # then we don't analyse it
        else:
            print('Info file found in folder:', subfol)          
            files = glob.glob(subfol+valid_fmt) # retrieves all files with valid format
            [q_step,Hmax,ks] =  fi.read_info(subfol+info_file)
            # Analyse files...
            [Z,Vxm,Vym,Vzm,
             vxvxn,vyvyn,vzvzn,
             vxvy,vxvz,vyvz,
             Tx,Ty,Tz,
             acc_rate,rot] = fi.ADV_profile_rot(subfol,SNRthreshold,CORthreshold,filt='GN2002W2003DV')
            ufs = np.nanmax(Vxm) # max non-nan value
            d_pos = np.nanargmax(Vxm) # First guess of delta (Z[d_pos]), to avoid fitting log-law in all the profile.
            [ushear_logwakelaw, delta_logwakelaw, res] = fi.find_ushear_delta_logwakelaw(Z[0:d_pos+1],Vxm[0:d_pos+1],ufs,PI,ks)
            # These two (ushear_logwakelaw, delta_logwakelaw) are the ones
            # that should be more robust.
            # Estimate q:
            q = fi.q_estimate(ushear_logwakelaw, delta_logwakelaw, ufs, Hmax, ks, PI)
            # Store the name of the folder analysed:
            subfol_analy = subfol
            
                                   
"""
From here downwards, just ad-hoc plotting for checkouts
"""
# Multi-Plotting ------------------------------------------------------------
# Mean velocities -----------------------------------------------------------

fig, axs = plt.subplots(3,figsize=(4,8), sharex=True, sharey=True)

title = 'Downstream cross-section, q-step = ' + str(q_step)
fig.suptitle(title)

axs[0].vlines(x=0,ymin=0, ymax=1, colors='k', linestyles='solid')
axs[0].vlines(x=0,ymin=Z.max()/delta_logwakelaw, ymax=1, colors='k', linestyles='dashed')
axs[0].scatter(Vxm/ufs, Z/delta_logwakelaw)
axs[0].set(xlabel='$\overline{V_x}/u_{fs}$ (-)', ylabel='$z/\delta$ (-)')
axs[1].vlines(x=0,ymin=0, ymax=1, colors='k', linestyles='solid')
axs[1].vlines(x=0,ymin=Z.max()/delta_logwakelaw, ymax=1, colors='k', linestyles='dashed')
axs[1].scatter(Vym/ufs, Z/delta_logwakelaw)
axs[1].set(xlabel='$\overline{V_y}/u_{fs}$ (-)', ylabel='$z/\delta$ (-)')
axs[2].vlines(x=0,ymin=0, ymax=1, colors='k', linestyles='solid')
axs[2].vlines(x=0,ymin=Z.max()/delta_logwakelaw, ymax=1, colors='k', linestyles='dashed')
axs[2].scatter(Vzm/ufs, Z/delta_logwakelaw)
axs[2].set(xlabel='$\overline{V_z}/u_{fs}$ (-)', ylabel='$z/\delta$ (-)')

fig.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()

# Reynolds shear stresses ---------------------------------------------------

fig, axs = plt.subplots(3,figsize=(4,8), sharex=True, sharey=True)

#title = 'Downstream cross-section, q = ' + str(q) + ' m$^2$/s'
fig.suptitle(title)

axs[0].vlines(x=0,ymin=0, ymax=1, colors='k', linestyles='solid')
axs[0].vlines(x=0,ymin=Z.max()/delta_logwakelaw, ymax=1, colors='k', linestyles='dashed')
axs[0].scatter(vxvy/ushear_logwakelaw**2, Z/delta_logwakelaw)
axs[0].set(xlabel='$\overline{v_x v_y}/u_*^2$ (-)', ylabel='$z/\delta$ (-)')
axs[1].vlines(x=0,ymin=0, ymax=1, colors='k', linestyles='solid')
axs[1].vlines(x=0,ymin=Z.max()/delta_logwakelaw, ymax=1, colors='k', linestyles='dashed')
axs[1].scatter(vxvz/ushear_logwakelaw**2, Z/delta_logwakelaw)
axs[1].set(xlabel='$\overline{v_x v_z}/u_*^2$ (-)', ylabel='$z/\delta$ (-)')
axs[2].vlines(x=0,ymin=0, ymax=1, colors='k', linestyles='solid')
axs[2].vlines(x=0,ymin=Z.max()/delta_logwakelaw, ymax=1, colors='k', linestyles='dashed')
axs[2].scatter(vyvz/ushear_logwakelaw**2, Z/delta_logwakelaw)
axs[2].set(xlabel='$\overline{v_y v_z}/u_*^2$ (-)', ylabel='$z/\delta$ (-)')

fig.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()

# Reynolds normal stresses --------------------------------------------------

fig, axs = plt.subplots(3,figsize=(4,8), sharex=True, sharey=True)

fig.suptitle(title)

axs[0].vlines(x=0,ymin=0, ymax=1, colors='k', linestyles='solid')
axs[0].vlines(x=0,ymin=Z.max()/delta_logwakelaw, ymax=1, colors='k', linestyles='dashed')
axs[0].scatter(vxvxn/ushear_logwakelaw**2, Z/delta_logwakelaw)
axs[0].set(xlabel='$\overline{v_x v_x}/u_*^2$ (-)', ylabel='$z/\delta$ (-)')
axs[1].vlines(x=0,ymin=0, ymax=1, colors='k', linestyles='solid')
axs[1].vlines(x=0,ymin=Z.max()/delta_logwakelaw, ymax=1, colors='k', linestyles='dashed')
axs[1].scatter(vyvyn/ushear_logwakelaw**2, Z/delta_logwakelaw)
axs[1].set(xlabel='$\overline{v_y v_y}/u_*^2$ (-)', ylabel='$z/\delta$ (-)')
axs[2].vlines(x=0,ymin=0, ymax=1, colors='k', linestyles='solid')
axs[2].vlines(x=0,ymin=Z.max()/delta_logwakelaw, ymax=1, colors='k', linestyles='dashed')
axs[2].scatter(vzvzn/ushear_logwakelaw**2, Z/delta_logwakelaw)
axs[2].set(xlabel='$\overline{v_z v_z}/u_*^2$ (-)', ylabel='$z/\delta$ (-)')

fig.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()

# Multiplots are finished by now ---------------------------------------------

# -------------------------------------
# Mean velocity plot - fitted parameters -------------------------------------
Ulog = fi.logwakelaw(Z[Z<delta_logwakelaw], delta_logwakelaw, ks, PI, ushear_logwakelaw)

plt.figure()
plt.title(title+' - Mean Velocity')
plt.plot(Ulog/ufs, Z[Z<delta_logwakelaw]/delta_logwakelaw, 'k-',label='Log-wake law, fitted parameters')
plt.scatter(Vxm/ufs, Z/delta_logwakelaw,label='ADV data')
plt.xlabel('$\overline{V_x}/u_{fs}$ (-)')
plt.ylabel('$Z/ \delta$  (-)')
plt.grid()
plt.legend(loc='best')
plt.show()

# Velocity defect law -------------------------------------------------------
plt.figure()

defect1 = fi.veldefect_logwakelaw(delta_logwakelaw, PI, Z)

plt.loglog(Z/delta_logwakelaw, defect1,  'r-',label='Log-wake law, fitted profile')

plt.scatter( Z/delta_logwakelaw,(ufs-Vxm)/ushear_logwakelaw,label='ADV data')

plt.xlabel('$Z/ \delta$')
plt.ylabel('$(u_{fs}-\overline{V_x})/u_{*}$')
plt.grid()
plt.legend(loc='best')
plt.show()

# Mean velocity profile with dimensions --------------------------------------

plt.figure()
plt.scatter(Vxm, Z, label='ADV data')

ZZ = np.linspace(Z.min(),Z.max(), 1000)

plt.semilogx(fi.loglaw(ZZ[ZZ<delta_logwakelaw], delta_logwakelaw, ks, ushear_logwakelaw) , 
           ZZ[ZZ<delta_logwakelaw], c='k', label='Log-wake law')

plt.xlabel('$\overline{V_x}$ (m/s)')
plt.ylabel('$Z$ (m)')
plt.legend(loc='best')
plt.tight_layout()
plt.show()

# Let's check the shear velocities --------------------------------

aa = -vxvz[Z<delta_logwakelaw]
bb = Z[Z<delta_logwakelaw]
# Siegel slope for a robust estimation of the data slope:
[medslope, medintercept] = stats.siegelslopes(aa[~np.isnan(aa)], bb[~np.isnan(aa)]) # (-vxvz[Z<delta_logwakelaw], Z[Z<delta_logwakelaw])

fig = plt.figure()
ax = fig.add_subplot(111)
vxvy_fit = medintercept + medslope * Z
delta_uxuz = -medintercept/medslope # BL thickness based on the uxuz law.

ax.vlines(x=0,ymin=0, ymax=1, colors='k', linestyles='solid')
ax.vlines(x=0,ymin=Z.max()/delta_logwakelaw, ymax=1, colors='k', linestyles='dashed')
ax.plot(vxvy_fit/medintercept, Z/delta_logwakelaw, c='steelblue',
        label='Robust fit to $\overline{v_x v_z}, \, z< \delta$' )
# theorical shear
Yshear = (1-Z/delta_logwakelaw)
ax.plot(Yshear[Z<delta_logwakelaw], Z[Z<delta_logwakelaw]/delta_logwakelaw, c='k',
        label='Linear distribution$' )

ax.scatter(-vxvz/medintercept, Z/delta_logwakelaw,
           label='$\overline{v_x v_z}$, dimless with $u_*$ \n (own)' )
ax.scatter(-vxvz/ushear_logwakelaw**2, Z/delta_logwakelaw,
           label='$\overline{v_x v_z}$, dimless with $u_*$ \n (from mean velocity gradient)' )

plt.xlabel('$\overline{v_x v_z}/u_*^2$ (-)')
plt.ylabel('$z/\delta$ (-)')
plt.legend(loc='best')
plt.show()

# -------------------------------------------------------------
# Missalignment of the ADV
plt.figure()
plt.scatter(180*rot/np.pi,Z/delta_logwakelaw)
plt.xlabel('Misalignment of ADV with streamwise component (deg)')
plt.ylabel('$z/\delta$ (-)')
plt.grid()

plt.show()

# Acceptance rates
plt.figure()
plt.scatter(acc_rate, Z/delta_logwakelaw)
plt.xlabel('Acceptance rate (%)')
plt.ylabel('$z/\delta$ (-)')
plt.grid()

plt.show()
# -------------------------------------------------------------

ushear_uxuz = np.sqrt(medintercept)
print("q: ", q)
print("ufs: ", ufs)

print("u* from log-wake law profile: ", ushear_logwakelaw)
print("u* from uxuz: ", ushear_uxuz)
print("u* ratio: ", ushear_logwakelaw/ushear_uxuz)

print("delta from log-wake law profile: ", delta_logwakelaw)
print("delta from uxuz: ", delta_uxuz)

print("Median angle of rotation before correction: ", np.nanmedian(180*rot/np.pi))
print("Median angle of rotation after correction: ", np.nanmedian(180*np.arctan2(Vym,Vxm)/np.pi))

print("Median acceptance rate (%):", np.median(acc_rate))

# Finally, save the data! ---------------------------------------------------
# Implement a data-saving routine
filename = subfol_analy.replace('\\', '_')

# Measurement depth
np.savetxt(filename+'Z.txt',Z)
# Velocities
np.savetxt(filename+'Vxm.txt',Vxm)
np.savetxt(filename+'Vym.txt',Vym)
np.savetxt(filename+'Vzm.txt',Vzm)
np.savetxt(filename+'rot.txt',rot)
# Reynolds normal stresses
np.savetxt(filename+'vxvxn.txt',vxvxn)
np.savetxt(filename+'vyvyn.txt',vyvyn)
np.savetxt(filename+'vzvzn.txt',vzvzn)
# Reynolds shear stresses
np.savetxt(filename+'vxvy.txt',vxvy)
np.savetxt(filename+'vxvz.txt',vxvz)
np.savetxt(filename+'vyvz.txt',vyvz)
# Integral timescales
np.savetxt(filename+'Tx.txt',Tx)
np.savetxt(filename+'Ty.txt',Ty)
np.savetxt(filename+'Tz.txt',Tz)

# Other turbulent quantities
np.savetxt(filename+'ushear.txt',[ushear_logwakelaw, ushear_uxuz])
np.savetxt(filename+'delta.txt',[delta_logwakelaw, delta_uxuz])

# Other hydrodynamic quantities
np.savetxt(filename+'q_Hmax_ks_ufs_PI.txt',[q, Hmax, ks, ufs, PI])

# Other ADV and analysis params 
np.savetxt(filename+'ADV_checks.txt',[SNRthreshold, CORthreshold])
np.savetxt(filename+'acc_rate.txt',acc_rate)

print("Finished analysis:", filename)

