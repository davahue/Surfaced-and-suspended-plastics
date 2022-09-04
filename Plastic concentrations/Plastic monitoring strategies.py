# -*- coding: utf-8 -*-
"""
Created on Wed Jan 12 15:21:14 2022

This function tests different riverine monitoring strategies and reports the
uncertainty and bias based on our laboratory experiences.

@author: Valero
"""

import os
import numpy as np
from scipy import stats
from scipy.interpolate import interp1d
from scipy import integrate
# import pandas as pd
from matplotlib import pyplot as plt

import Concentration_profiles_faux as cp


def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())

def relmedabse(predictions, targets):
    return 100*np.median(abs(predictions-targets)/targets)

def relmeane(predictions, targets):
    return 100*np.mean((predictions-targets)/targets)

def relerror(predictions, targets):
    return (100*(predictions-targets)/targets)

def absrelerror(predictions, targets):
    return np.abs(100*(predictions-targets)/targets)

def r_val(predictions, targets):
    slope, intercept, r_value, p_value, std_err = stats.linregress(targets, predictions)
    return r_value

def error_report(predicted, targets):
    
    print("r2:", r_val(predicted, targets)**2)
    # print("RMSE:", rmse(predicted, targets))
    print("Relative Mean Err: ", relmeane(predicted, targets))
    print("Relative Median Abs Err: ", relmedabse(predicted, targets))
    print("----------------------------")
    
    return

filelist = os.listdir()

files_zp = []
files_surf = []
files_scalars = []

ushear = [] # shear velocity
w = [] # rising velocity
lmax = []
lortho = []
M = []
rhop = []

kappa = 0.41 # -
Hmax = 0.278 # m
sigma = 0.072 # N/m
rhow = 998 # kg/m3
g = 9.81 # m/s2


for file in filelist:
    if "_zp" in file:
        files_zp.append(file)
        us, ws = cp.get_flow_sample_chars(file)
        ushear.append(us)
        w.append(ws)
        lmaxi, lorthoi, Mi, rhopi = cp.get_sample_geomchars(file)
        lmax.append(lmaxi)
        lortho.append(lorthoi)
        M.append(Mi)
        rhop.append(rhopi)
        
    elif "_surfp" in file:
        files_surf.append(file)
    elif "_scalars" in file:
        files_scalars.append(file)

# list to np.array ----
ushear = np.asarray(ushear); w = np.asarray(w);
lmax = np.asarray(lmax)/1000 # retrieved in mm
lortho = np.asarray(lortho)/1000 # retrieved in mm;
M = np.asarray(M); rhop = np.asarray(rhop)

a = []; a_surfaced_video = []; a_surfaced_KS = [] # a-level estimations
Ntotal = [] # Total count of particles
Nsus = [] # Count particles in the suspended layer
Nsurf = [] # Count in the surface layer, suspended or touching the free surface
Nsurfsurf = [] # in the surface layer, suspended or touching the free surface
Nsurfsus = [] # in the surface layer, yet travelling through suspension
beta = []
Ca = []

# Str. 2
Np_str2 = [] # Estimation of plastic budget with 20 % depth visibility

# Str. 5
dz = Hmax/10 # 

CX = [] # for the estimation of the total budget including suspension
CX_factor_X = 1-0.50 # 1 minus the depth percentage at measurement

# Str. 6
CX2 = []
CX_factor_X2 = 1-0.50 # 1 minus the depth percentage at measurement

# For Rouse-based corrections (Str. 7)
CX_factor_XRouse7 = 1-0.50
CXRouse7 = []
Ksus7 = []


# For Rouse-based corrections (Str. 8)
CX_factor_XRouse8 = 1-0.60
CXRouse8 = []
Ksus8 = []

# For plotting:
marker = []
label = []

for i in range(0,len(files_zp)):
    
    [mi, li] = cp.get_marker_label(files_surf[i])
    marker.append(mi)
    label.append(li)
    zp = np.loadtxt(files_zp[i])/100
    surfp = np.loadtxt(files_surf[i])

    [plane, betai, binsize, Cai, a_surfaced_videoi, a_surfaced_KSi] = np.loadtxt(files_scalars[i])
    beta.append(betai) # betai = wi/(kappa*usheari)
    Ca.append(Cai)   

    # Compute budgets --------------------------------------------------------
    
    ai = -a_surfaced_videoi/100
    a_surfaced_video.append(ai)
    a_surfaced_KS.append(-a_surfaced_KSi/100)
    a.append(ai)
    
    Ntotal.append(len(zp)) # Total plastics count
    Nsus.append(np.sum(zp<-ai)) # Plastics below level a (a positive now)
    Nsurf.append(np.sum(zp>=-ai))
    # Nsurf.append(Ntotal[i]-Nsus[i])
    
    # Within the surface transport layer:
    Nsurfsurf.append(np.sum(surfp==1)) # Plastics identified as touching the free surface
    Nsurfsus.append(Nsurf[i] - Nsurfsurf[i])
    
    COUNT3 = np.sum((zp>=-0.20*Hmax) * (surfp==0))
    Np_str2.append(COUNT3 + np.sum(surfp==1))
    
    # Now, considering suspension
    
    zp[zp>=0] = 0.0 # Anything over the free surface (floating) would be trapped by surface skimming
    
    binno = round((0.00+Hmax)/dz)
    aa, bb = np.histogram(zp, bins=binno, range=[-Hmax, 0.00])
    bbb = (bb[:-1] + bb[1:])/2 # This means, a bin per 10 % of the flow depth.
                                # Thus, suspension sampling is representative 
                                # of the count in a 10 % water column depth slot.
    f = interp1d(bbb, aa) # Want a continuous function to interpolate at X % of the water depth
    # For surface observation
    CX.append(f(-CX_factor_X*Hmax))
    # For skimmed flow in the surface
    CX2.append(f(-CX_factor_X2*Hmax))
    
    
    # For Rouse-based correction:
        # Str. 8
    CXRouse8.append(f(-CX_factor_XRouse8*Hmax))
        
    zz = np.linspace(-0.1999999*Hmax, -Hmax, 100)    
    CC = cp.Rouse(-zz, 0.1999999*Hmax, 1, Hmax, betai)
    CCinterp = interp1d(zz, CC)
    Ksusi = integrate.trapz(CC,x=-zz)/(integrate.trapz(-zz,x=-zz)*CCinterp(-CX_factor_XRouse8*Hmax))
    Ksus8.append(Ksusi)
    
    # Str. 7
    CXRouse7.append(f(-CX_factor_XRouse7*Hmax))
    zz7 = np.linspace(-0.1*Hmax, -Hmax, 100)
    CC7 = cp.Rouse(-zz7, 0.1*Hmax, 1, Hmax, betai)
    CCinterp7 = interp1d(zz7, CC7)
    Ksusi = integrate.trapz(CC7,x=-zz7)/(integrate.trapz(-zz7,x=-zz7)*CCinterp(-CX_factor_XRouse7*Hmax))
    Ksus7.append(Ksusi)
    
    
Ntotal = np.asarray(Ntotal); Nsus = np.asarray(Nsus); Nsurf = np.asarray(Nsurf)   
Nsurfsurf = np.asarray(Nsurfsurf); Nsurfsus = np.asarray(Nsurfsus)
a = np.asarray(a); a_surfaced_video = np.asarray(a_surfaced_video); a_surfaced_KS = np.asarray(a_surfaced_KS)
beta = np.asarray(beta); Ca = np.asarray(Ca)

# Dimensionless relationships ------------------------------------------------
AUX1 = 8*sigma/(rhow*ushear**2)
AUX2 = (lmax+lortho)/(lmax*lortho)
Gamma = AUX1*AUX2

# For Str. 3 -- Eq. 5-based correction
a1 = 6.661; b1 = 0.6432; c1 = 2.37
fxy = a1*(np.exp(-beta/b1) + np.exp(-Gamma/c1))+1
Np_est_Eq5 = fxy*Nsurfsurf


# Strategy 1
Np_str1 = Nsurfsurf
# Strategy 2
Np_str2 = np.asarray(Np_str2)

# Strategy 4: Visual observation of the upper 20 % and estimation of suspension with Eq. (6)
fraction_sus = np.exp(-beta/1.61)*np.exp(-Gamma/38.08)
Np_str4 = Np_str2 / (1-fraction_sus)

# Strategy 5: Count at the free surface and sampling point in suspension (uniform).
CX = np.asarray(CX)
Np_str5 = Np_str1 + CX*binno
# Strategy 6: Visual observation of the upper 20 % and sampling point in suspension (uniform). 
CX2 = np.asarray(CX2)
Np_str6 = Np_str2 + 0.80*CX2*binno # We assume that 20 % of the water column 
                                   # is surface transport and 80 % of the 
                                   # water column is suspension (with 
                                   # uniform concentration).

# Strategy 8: skimming 20 % flow and Rouse profile assumption at suspension point
CXRouse8 = np.asarray(CXRouse8)
Np_str8 = Np_str2 + CXRouse8*Ksus8

CXRouse7 = np.asarray(CXRouse7)
Np_str7 = Np_str1 + CXRouse7*Ksus7

# Plot the different monitoring strategies --------------------------------

plt.figure(figsize=(7,4))

# Strategy 1: Visual observation at the free surface is representative of all plastic budget

for i in range(0,len(beta)):
    plt.scatter(beta[i], relerror(Np_str1[i], Ntotal[i]),
                s = 30, lw=2, facecolors='brown', edgecolor='r',
                marker=marker[i])
plt.scatter(beta[i], relerror(Np_str1[i], Ntotal[i]),
            s = 30, lw=2,facecolors='brown', edgecolor='r',
            marker=marker[i], label='Str. 1') # : Visual \nobservation with \nlimitted visibility

print("ERROR REPORT -- Strategy 1")
error_report(Np_str1, Ntotal)

# Strategy 2: Visual observation of the free surface which includes 20 % depth visibility (equivalent to skimming the surface layer)

for i in range(0,len(beta)):
    plt.scatter(beta[i], relerror(Np_str2[i], Ntotal[i]),
                s=30, lw=2,marker=marker[i], facecolors='orange',
            edgecolor='darkorange')

plt.scatter(beta[i], relerror(Np_str2[i], Ntotal[i]),
            s=30, lw=2,marker=marker[i], facecolors='orange',
            edgecolor='darkorange', label='Str. 2') # : Visual \nobservation with \nimproved visibility

print("ERROR REPORT -- Strategy 2")
error_report(Np_str2, Ntotal)

# Strategy 3: Visual observation at the free surface corrected by the expected fraction of surfaced plastics (Eq. (5))

for i in range(0, len(beta)):
    plt.scatter(beta[i], relerror(Np_est_Eq5[i], Ntotal[i]), marker=marker[i],
                facecolors='none', edgecolor='k')
plt.scatter(beta[i], relerror(Np_est_Eq5[i], Ntotal[i]), marker=marker[i],
                facecolors='none', edgecolor='k',label="Str. 3")
print("ERROR REPORT -- Strategy 3")
error_report(Np_est_Eq5, Ntotal)

# Strategy 4: Visual observation of the upper 20 % and estimation of suspension with Eq. (6)

for i in range(0, len(beta)):
    plt.scatter(beta[i], relerror(Np_str4[i], Ntotal[i]), marker=marker[i],
                facecolors='beige', edgecolor='darkkhaki')
plt.scatter(beta[i], relerror(Np_str4[i], Ntotal[i]), marker=marker[i],
                facecolors='beige', edgecolor='darkkhaki',label="Str. 4") # : Str B \nand estimation \nof suspension (Eq. 6)
print("ERROR REPORT -- Strategy 4")
error_report(Np_str4, Ntotal)

# Strategy 5: Count at the free surface and point in suspension (uniform).

for i in range(0, len(beta)):
    plt.scatter(beta[i], relerror(Np_str5[i], Ntotal[i]), marker=marker[i],
                facecolors='0.85', edgecolor='k')
plt.scatter(beta[i], relerror(Np_str5[i], Ntotal[i]), marker=marker[i],
                facecolors='0.85', edgecolor='k',label="Str. 5")
print("ERROR REPORT -- Strategy 5")
error_report(Np_str5, Ntotal)

# Strategy 6: Visual observation of the upper 20 % and sampling point in suspension (uniform). 

for i in range(0, len(beta)):
    plt.scatter(beta[i], relerror(Np_str6[i], Ntotal[i]), marker=marker[i],
                facecolors='0.65', lw=1, edgecolor='k')
plt.scatter(beta[i], relerror(Np_str6[i], Ntotal[i]), marker=marker[i],
                facecolors='0.65', lw=1, edgecolor='k',label="Str. 6")
print("ERROR REPORT -- Strategy 6")
error_report(Np_str6, Ntotal)

# Strategy 7: Visual observation of the upper 20 % and sampling point in suspension (Rouse, starting at 20 % depth). 

for i in range(0, len(beta)):
    plt.scatter(beta[i], relerror(Np_str7[i], Ntotal[i]), marker=marker[i],
                facecolors='0.45', lw=1, edgecolor='k')
plt.scatter(beta[i], relerror(Np_str7[i], Ntotal[i]), marker=marker[i],
                facecolors='0.45', lw=1, edgecolor='k',label="Str. 7")
print("ERROR REPORT -- Strategy 7")
error_report(Np_str7, Ntotal)

# Strategy 8: Visual observation of the upper 20 % and sampling point in suspension (Rouse). 

for i in range(0, len(beta)):
    plt.scatter(beta[i], relerror(Np_str8[i], Ntotal[i]), marker=marker[i],
                facecolors='turquoise', lw=2, edgecolor='forestgreen')
plt.scatter(beta[i], relerror(Np_str8[i], Ntotal[i]), marker=marker[i],
                facecolors='turquoise', lw=2, edgecolor='forestgreen',label="Str. 8") # : Str B \nand sampling \nof suspension (Rouse \nprofile)
print("ERROR REPORT -- Strategy 8")
error_report(Np_str8, Ntotal)

# Null error line
plt.hlines(y = 0,xmin=beta.min(),xmax=beta.max(), color='b')

# Wide-range overview:
plt.semilogx()
# plt.yscale('symlog')

plt.ylabel('Plastic budget relative error (%)')
plt.xlabel('$\\beta $ (-)')
plt.legend(loc='best',ncol=2)

plt.savefig('Plastic_monitoring_rel_error.pdf', dpi=600, format='pdf',  bbox_inches='tight')
plt.savefig('Plastic_monitoring_rel_error.png', dpi=600, format='png',  bbox_inches='tight')
plt.savefig('Plastic_monitoring_rel_error.svg', dpi=600, format='svg',  bbox_inches='tight')

plt.show()

# ----------------------------------------------------------------------------

number_of_bins = 20

# Data sets

labels = ["Str. 1", "Str. 2", "Str. 3", "Str. 4", "Str. 5", "Str. 6", "Str. 7", "Str. 8"]
data_sets = [relerror(Np_str1, Ntotal),
             relerror(Np_str2, Ntotal),
             relerror(Np_est_Eq5, Ntotal), # Str. 3
             relerror(Np_str4, Ntotal),
             relerror(Np_str5, Ntotal),
             relerror(Np_str6, Ntotal),
             relerror(Np_str7, Ntotal),
             relerror(Np_str8, Ntotal)]


hist_range = (np.min(data_sets), np.max(data_sets))
binned_data_sets = [
    np.histogram(d, range=hist_range, bins=number_of_bins)[0]
    for d in data_sets
]

binned_maximums = np.max(binned_data_sets, axis=1)
x_locations = np.arange(0, sum(binned_maximums), np.max(binned_maximums))
x_locations = np.arange(0, len(binned_maximums)*np.max(binned_maximums), np.max(binned_maximums))

# The bin_edges are the same for all of the histograms
bin_edges = np.linspace(hist_range[0], hist_range[1], number_of_bins + 1)
centers = 0.5 * (bin_edges + np.roll(bin_edges, 1))[:-1]
heights = np.diff(bin_edges)

# Start plotting
fig, ax = plt.subplots()
# Color examples
# FC = ['brown', 'orange', 'beige', 'turquoise']
# EC = ['r', 'darkorange', 'darkkhaki', 'forestgreen']
i = 0
for x_loc, binned_data in zip(x_locations, binned_data_sets):
    lefts = x_loc - 0.5 * binned_data
    ax.barh(centers, binned_data, height=heights, left=lefts) # , facecolor=FC[i], edgecolor=EC[i])
    i = i+1
ax.set_xticks(x_locations)
ax.set_xticklabels(labels)

#ax.set_ylim(-100, 50)
ax.set_ylabel("Rel. error (%)")
# ax.set_xlabel("")
plt.hlines(y = 0,xmin=x_locations.min(),xmax=x_locations.max(), color='b')

plt.show()

# ----------------------------------------------------------------------------
# Boxplot for all strategies

plt.figure(figsize=(3.,2.0))

plt.hlines(0, xmin=0.5, xmax=8.5, color='b')
plt.boxplot(data_sets, notch=True, sym ="+")

plt.xlabel("Str.")
plt.ylabel("Rel. error (%)")

plt.savefig('Plastic_monitoring_boxplot.pdf', dpi=600, format='pdf',  bbox_inches='tight')
plt.savefig('Plastic_monitoring_boxplot.png', dpi=600, format='png',  bbox_inches='tight')
plt.savefig('Plastic_monitoring_boxplot.svg', dpi=600, format='svg',  bbox_inches='tight')

