# -*- coding: utf-8 -*-
"""
Created on Wed Nov  3 15:06:56 2021

This code analyzes particles concentrations extracted from 3D trajectory 
reconstructions.

@author: D. Valero
"""

#import os, glob
import numpy as np
import scipy.stats as stats
import pandas as pd
from matplotlib import pyplot as plt

# In-house libraries
import Concentration_profiles_faux as cp

# Constants definition -------------------------------------------------------
# Geom
BedLev = 36.3# cm
WaterLev = BedLev-27.8 # cm
Hmax = BedLev-WaterLev
# Hydro
kappa = 0.41

# Choose a file --------------------------------------------------------------
# folder = 'Masks_trajectories'
# filename = 'Mask_V5_d5.csv'

folder = 'Cups_trajectories'
filename = 'Cup_PP_98_def_V3_d5.csv'

# folder = 'Films_trajectories'
# filename = 'Film_HDPE_15_V5_d5.csv'

path = folder + '\\' + filename

# Load data ------------------------------------------------------------------
data = pd.read_csv(path)
print('Accessed:',path)
# Access shear velocity and rising velocity
ushear, w = cp.get_flow_sample_chars(filename)
beta = w/(kappa*ushear)

# Start ANALYSIS -------------------------------------------------------------
# Retrieve stats from trajectory data
plane = 55 # cm, plane = 55 cm for midplane of the analysis window.
[tp, xp, yp, zp, vpx, vpy, vpz, surfp] = cp.extract_particle_stats(data, plane)
zp = zp + WaterLev # From here on, zero at the free surface (i.e.: zp = z-h)

for i in range(0,len(zp)):
    if "HDPE_15_" in filename:
        if zp[i] < -6.29: # Max dimension (lambda)
            surfp[i] = 0
    elif "HDPE_100_" in filename:
        if zp[i] < -16.6: # Max dimension (lambda)
            surfp[i] = 0
    elif "Cup_PP_100" in filename:
        if zp[i] < -10.8: # Max dimension (lambda)
            surfp[i] = 0
    elif "Cup_PP_98" in filename:
        if zp[i] < -11.4: # Max dimension (lambda)
            surfp[i] = 0
    elif "Cup_PP_50" in filename:
        if zp[i] < -10.8: # Max dimension (lambda)
            surfp[i] = 0
    elif "Cup_PP_05" in filename:
        if zp[i] < -4.82: # Max dimension (lambda)
            surfp[i] = 0

a_surfaced_video = np.min(zp[surfp==1])

# Initialize KS-variables
zpsorted = np.sort(zp)
Nsamples = np.zeros(len(zpsorted))
st = np.zeros(len(zpsorted))
pv = np.zeros(len(zpsorted))
# Conduct KS tests
for i in range(5, len(zpsorted)-1):
    ai = zpsorted[i]-0.0001 # So we always include the last sample
    Nsamples[i] = len(zpsorted[0:i])
    ZRouse_true = -cp.inv_sampling_Rouse(-ai,1,Hmax,beta,10000)
    [st[i],pv[i]] = stats.ks_2samp(zpsorted[0:i], ZRouse_true, alternative='two-sided')

plt.figure()
plt.hist(ZRouse_true,bins=200,cumulative=True, histtype='step',color='k', density=True)
plt.hist(zpsorted,bins=50,cumulative=True, histtype='step', density=True);
plt.xlabel('$z-H$ (cm)')
plt.ylabel('Cummulative Density Function')

plt.show()

# Find KS-based a_surfaced
significance_indexes = np.where(np.diff(np.sign(pv-1E-3)))[0]
i_ = significance_indexes[1]+1
if i_ < len(zpsorted)-1:
    a_surfaced_KS = zpsorted[i_]
else:
    a_surfaced_KS = np.nan
# We ignored the first values due to insufficient number of samples. Take care of it:
st[Nsamples==0] = np.nan; pv[Nsamples==0] = np.nan; Nsamples[Nsamples==0] = np.nan; 
# Check if this is statistically significant --------------------------------
# PLOTS
cp.plot_significance(Nsamples, st, pv, zpsorted, a_surfaced_video, a_surfaced_KS, i_,
                      BedLev, WaterLev, filename)



# ---------------------------------------------------------------------------
# Visualization: vertical count profiles
# Experimental data -----------------------
plt.figure(figsize=(3.0,3.0))
# Settings ----------------------------------
plt.axhline(y=0.0, color='b', linestyle='-', label='Free surface')
plt.axhline(y=0.0+a_surfaced_video, color='r', linestyle='-.', label='$- a$')
if np.isnan(a_surfaced_KS):
    print("KS test did not detect any a_KS.")
else:
    plt.axhline(y=0.0+a_surfaced_KS, color='r',  lw=3, linestyle=':', label='$- a_{KS}$')
plt.axhline(y=-BedLev+WaterLev, color='k', linestyle='--', label='Roughness crests')

plt.ylabel('$z - H$ (cm)') # , fontsize=8
plt.xlabel('count (-)') # , fontsize=8
plt.tight_layout()


# To be used in the Rouse profile:
a = a_surfaced_video

# Define a bin size, and calculate the corresponding Ca for Rouse
binsize = 1.0 # 1.00 # 2.5 # cm
bins=np.arange(-BedLev+WaterLev, 0+binsize, binsize)
[hist_count, hist_bins, axs] = plt.hist(zp,bins,color='k', lw=1, ls="-",alpha=1, histtype="step",
         orientation="horizontal", label="Total particles\' C.o.G. \ncount")
plt.hist(zp[surfp==1],bins,color='firebrick',alpha=0.5, histtype="stepfilled",
         orientation="horizontal", label='Surfaced particles\' C.o.G. \ncount')

# Total surfaced and suspended transport -----------
# from 0 to a: surfaced. Beyond: suspended.
z = (hist_bins[0:-1]+hist_bins[1:])/2
binlist = hist_bins[1:]<a # lower side of the bin is below of "a".
Csus = np.sum(hist_count[binlist]) # Where binlist==TRUE, suspended transport.

try:
    Ca = cp.get_Ca(z, -a, Hmax, beta, Csus)
    # Get the Rouse profile with physically-based beta and fit Ca.
    C = cp.Rouse(-z, -a, Ca, Hmax, beta)

    plt.scatter(C, z, lw=2, zorder=5, marker="+",c="sandybrown", 
                label="Rouse profile, $\\beta=w/\\kappa u_*$ ") 
except Exception:
    print("Could not retrieve a Ca value for Rouse profile plotting purposes.")
    C = np.zeros_like(z)
    Ca = 0

print("--------------------------------------")
print("Suspended particles counted:", int(Csus))
print("Suspended particles with Rouse:", round(np.nansum(C)))
print("--------------------------------------")

# Save plot ------------------------------------------------------------------
plotname = filename[:-4]+'_CountProfile'
plt.savefig(plotname+'.pdf', dpi=600, format='pdf',  bbox_inches='tight')
plt.savefig(plotname+'.png', dpi=600, format='png',  bbox_inches='tight')
plt.savefig(plotname+'.svg', dpi=600, format='svg',  bbox_inches='tight')

plt.legend(loc='lower right') #, fontsize=8
plt.savefig(plotname+'_legend.pdf', dpi=600, format='pdf',  bbox_inches='tight')
plt.savefig(plotname+'_legend.png', dpi=600, format='png',  bbox_inches='tight')
plt.savefig(plotname+'_legend.svg', dpi=600, format='svg',  bbox_inches='tight')

# Save data ------------------------------------------------------------------
# VECTORS:
    # midplane
np.savetxt(filename[:-4]+'_tp.txt',tp)
np.savetxt(filename[:-4]+'_xp.txt',xp)
np.savetxt(filename[:-4]+'_yp.txt',yp)
np.savetxt(filename[:-4]+'_zp.txt',zp)
np.savetxt(filename[:-4]+'_vpx.txt',vpx)
np.savetxt(filename[:-4]+'_vpy.txt',vpy)
np.savetxt(filename[:-4]+'_vpz.txt',vpz)
np.savetxt(filename[:-4]+'_surfp.txt',surfp)
    # KS test
np.savetxt(filename[:-4]+'_Nsamples.txt',Nsamples)
np.savetxt(filename[:-4]+'_st.txt',st)
np.savetxt(filename[:-4]+'_pv.txt',pv)

# SCALARS:
scalars = [plane, beta, binsize, Ca, a_surfaced_video, a_surfaced_KS]
np.savetxt(filename[:-4]+'_scalars.txt',scalars)  

print("a_KS: ", -a_surfaced_KS)
print("a_cam: ", -a_surfaced_video)
