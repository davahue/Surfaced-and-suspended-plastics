# -*- coding: utf-8 -*-
"""
Created on Wed Jan 12 15:21:14 2022

@author: Valero
"""

import os
import numpy as np
from scipy import stats
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.ticker as mticker

import Concentration_profiles as cp

def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())


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
    beta.append(betai)
    Ca.append(Cai)   

    # Compute budgets --------------------------------------------------------
    
    ai = -a_surfaced_videoi/100
    a_surfaced_video.append(ai)
    a_surfaced_KS.append(-a_surfaced_KSi/100)
    a.append(ai)
    
    Ntotal.append(len(zp)) # Total particle count
    Nsus.append(np.sum(zp<-ai)) # Particles below level a (a positive now)
    Nsurf.append(np.sum(zp>=-ai))
    
    # Within the surface transport layer:
    Nsurfsurf.append(np.sum(surfp==1)) # Particles identified as touching the free surface
    Nsurfsus.append(Nsurf[i] - Nsurfsurf[i])
    
    bins= np.linspace(-0.25,0,40)
    plt.hist(zp, density = True, bins=bins, color='k', histtype='step', alpha=0.15, orientation='horizontal')
    plt.hist(zp[surfp==0]/ai, density = True, color='b', alpha=0.15, bins=bins,orientation='horizontal');

Ntotal = np.asarray(Ntotal); Nsus = np.asarray(Nsus); Nsurf = np.asarray(Nsurf)   
Nsurfsurf = np.asarray(Nsurfsurf); Nsurfsus = np.asarray(Nsurfsus)
a = np.asarray(a); a_surfaced_video = np.asarray(a_surfaced_video); a_surfaced_KS = np.asarray(a_surfaced_KS)
a_surfaced_KS[a_surfaced_KS==1E-5] = np.nan # 0.01 mm correspond to not detected.
beta = np.asarray(beta); Ca = np.asarray(Ca)

# Hydrodynamic dimensionless relationships -----------------------------------
# plastic-based inverse Weber number
AUX1 = 8*sigma/(rhow*ushear**2)
AUX2 = (lmax+lortho)/(lmax*lortho)
Gamma = AUX1*AUX2

# plastic-based Bond number
Deltarho = (rhow-rhop)/rhop
Lambda = Deltarho*g*M/(sigma*(lmax+lortho))
# ----------------------------------------------------------------------------
# Concentrations -------------------------------------------------------------

Csurf = Nsurf/ Ntotal
Csus = Nsus/ Ntotal
Ctotal = Ntotal/ Ntotal
Csurfsus = Nsurfsus/ Ntotal
Csurfsurf = Nsurfsurf/ Ntotal


# # ----------------------------------------------------------------------------

plt.figure(figsize=(4.5,3.5))
plt.scatter(Gamma, Lambda, c=beta, marker='s', edgecolor='k')
plt.colorbar(label='$\\beta$ (-)')
plt.xlabel("$\\Gamma$ (-)")
plt.ylabel("$\\Lambda$ (-)")

plt.semilogx()
plt.semilogy()

plt.tight_layout()

plt.savefig('Gamma_Lambda_beta_map.pdf', dpi=600, format='pdf',  bbox_inches='tight')
plt.savefig('Gamma_Lambda_beta_map.png', dpi=600, format='png',  bbox_inches='tight')
plt.savefig('Gamma_Lambda_beta_map.svg', dpi=600, format='svg',  bbox_inches='tight')

# -----

plt.figure(figsize=(4.5,3.5))
plt.scatter(Gamma, Lambda, c=Nsus/Ntotal, marker='s', edgecolor='k')
plt.colorbar(label='$C_{\\beta}$ (-)') # $N_{\\beta}/N_{p}$
plt.clim(0,1)

plt.xlabel("$\\Gamma$ (-)")
plt.ylabel("$\\Lambda$ (-)")

plt.semilogx()
plt.semilogy()

plt.tight_layout()

plt.savefig('Gamma_Lambda_Cbeta_map.pdf', dpi=600, format='pdf',  bbox_inches='tight')
plt.savefig('Gamma_Lambda_Cbeta_map.png', dpi=600, format='png',  bbox_inches='tight')
plt.savefig('Gamma_Lambda_Cbeta_map.svg', dpi=600, format='svg',  bbox_inches='tight')

# -----
plt.figure(figsize=(4.5,3.5))
plt.scatter(Gamma, Lambda, c=Nsurfsurf/Ntotal, marker='s', edgecolor='k')
plt.colorbar(label='$C_{\\Gamma, surf}$ (-)')
plt.clim(0,1)

plt.xlabel("$\\Gamma$ (-)")
plt.ylabel("$\\Lambda$ (-)")

plt.semilogx()
plt.semilogy()

plt.tight_layout()

plt.savefig('Gamma_Lambda_Csurfsurf_map.pdf', dpi=600, format='pdf',  bbox_inches='tight')
plt.savefig('Gamma_Lambda_Csurfsurf_map.png', dpi=600, format='png',  bbox_inches='tight')
plt.savefig('Gamma_Lambda_Csurfsurf_map.svg', dpi=600, format='svg',  bbox_inches='tight')


data = {'Filename': files_zp,
        '$\\Gamma$ (-)': Gamma,
        '$\\Lambda$ (-)':Lambda,
        '$\\beta$ (-)': beta,
        '$a/H$ (-)': a/Hmax,
        '$a$ (m)': a,
        '$l_{max}$ (mm)': lmax,
        '$l_{ortho}$ (mm)': lortho,
        '$M$ (g)': M,
        '$\\rho_p$ (kg/m3)': rhop,
        '$w$ (m/s)': w,
        '$u_*$ (m/s)': ushear,
        '$N_p$ (-)': Ntotal,
        '$N_{\\beta} (-)$': Nsus,
        '$N_{\\Gamma} (-)$': Nsurf,
        '$N_{\\Gamma,sus} (-)$': Nsurfsus,
        '$N_{\\Gamma,surf} (-)$': Nsurfsurf,
        '$N_{\\Gamma}/N_p$ (-)': Nsurf/Ntotal,
        '$C_{p}$ (-)': Ctotal,
        '$C_{\\Gamma} (-)$': Csurf,
        '$C_{\\beta}$ (-)': Csus,
        '$C_{\\Gamma,sus}$ (-)': Csurfsus,
        '$C_{\\Gamma,surf}$ (-)': Csurfsurf
        }


df = pd.DataFrame(data=data)
df.to_csv("00_Data_summary.csv")

# Start of surface transport -------------------------------------------------

L = np.sqrt(lmax*lmax + lortho*lortho)

plt.figure(figsize=(3.4,3.15))
for i in range(0, len(Gamma)):
    plt.scatter(beta[i], a[i]/L[i], s=40, edgecolor='k', facecolor='k',
            marker=marker[i],zorder=4)

    
plt.scatter(beta[i], a[i]/L[i], s=40, edgecolor='k', facecolor='k', # s=300*a/Hmax
        marker=marker[i],label='$a/\\lambda$ (-)')

FITx = np.linspace(0, 12, 100)
FITy = 0.36+0.14*np.tanh(-FITx+3.2)
plt.plot(FITx,FITy,'b--', label='Eq. 4')


plt.xlabel('$\\beta $ (-)')
plt.ylabel('$a/\\lambda $ (-)')

plt.loglog()
#plt.semilogy()
plt.legend()
ax=plt.gca()
ax.yaxis.set_minor_formatter(mticker.ScalarFormatter())

plt.savefig('a_detection.pdf', dpi=600, format='pdf',  bbox_inches='tight')
plt.savefig('a_detection.png', dpi=600, format='png',  bbox_inches='tight')
plt.savefig('a_detection.svg', dpi=600, format='svg',  bbox_inches='tight')

"""
Print predictions accuracy:
    
"""
predicted = L*(0.36+0.14*np.tanh(-beta+3.2))
slope, intercept, r_value, p_value, std_err = stats.linregress(a, predicted)
print("Prediction of a, r2: ", r_value**2)
print("RMSE: ", rmse(predicted, a))
print("Median Abs Err: ", np.median(abs(predicted-a)))
print("----------------------------")


# ----------------------------------------------------------------------------
# *beta* is the right parameter to combine with either Gamma or Lambda. Gamma has lesser uncertainty in its determination.

plt.figure(figsize=(4.,3.5))
for i in range(0, len(Gamma)):
    plt.scatter(Gamma[i], beta[i], c=Csurfsurf[i],s=40, edgecolor='w', # s=300*a/Hmax
            marker=marker[i],cmap='viridis', vmin=0, vmax=1, zorder=9)
plt.loglog()
# -----

DG = 0.05*(np.max(Gamma)-np.min(Gamma))
bL = 0.05*(np.max(beta)-np.min(beta))

Gamma_pr = np.linspace(np.min(Gamma)-DG, np.max(Gamma)+DG, 1000)
beta_pr = np.linspace(np.min(beta)-bL, np.max(beta)+bL, 1000)
GG, bb = np.meshgrid(Gamma_pr,beta_pr)  

a1 = 6.661; b1 = 0.6432; c1 = 2.37
fxy = a1*(np.exp(-bb/b1) + np.exp(-GG/c1))+1

C_surfsurf_pred = 1/fxy


# ------

CS = plt.pcolormesh(GG,bb,C_surfsurf_pred, cmap='viridis',shading='auto',
                    vmin=0, vmax=1)
CS.set_rasterized(True)

lvls = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
CS_lines = plt.contour(GG,bb,C_surfsurf_pred, colors = 'w', linestyles = '--',
                    levels=lvls)
plt.clabel(CS_lines, lvls[1::2], inline=True, fontsize=8, fmt='%1.1f',)

plt.colorbar(CS, shrink=0.40, location='bottom',anchor=(1,0), label='$C_{\\Gamma,surf}$ (-)')

plt.xlabel("$\\Gamma$ (-)")
plt.ylabel("$\\beta$ (-)")

plt.xlim([np.min(Gamma)-DG, np.max(Gamma)+DG])
plt.ylim([np.min(beta)-bL, np.max(beta)+bL])

plt.tight_layout()

plt.savefig('Gamma_beta_Csurfsurf_raster.pdf', dpi=600, format='pdf',  bbox_inches='tight')
plt.savefig('Gamma_beta_Csurfsurf_raster.png', dpi=600, format='png',  bbox_inches='tight')
plt.savefig('Gamma_beta_Csurfsurf_raster.svg', dpi=600, format='svg',  bbox_inches='tight')


plt.show()
# ----------------------------------------------------------------------------
# Subfigure -- internal

plt.figure(figsize=(2.5,2))

fxy = a1*(np.exp(-beta/b1) + np.exp(-Gamma/c1))+1

predicted = 1/fxy


for i in range(0, len(Gamma)):
    plt.scatter(Csurfsurf[i], predicted[i], marker= marker[i], edgecolor='b',
                facecolor='None')
plt.plot([0,1],[0,1], 'k-')
plt.plot([0.1,1],[0,0.9], 'k--')
plt.plot([0,0.9],[0.1,1], 'k--')
plt.plot([0.25,1],[0,0.75], 'k:')
plt.plot([0,0.75],[0.25,1], 'k:')

plt.xlabel("$C_{\\Gamma,surf}$ (-),\n experimental")
plt.ylabel("$C_{\\Gamma,surf}$ (-),\n modelled")

plt.xlim([0, 1])
plt.ylim([0, 1])
plt.tight_layout()

plt.savefig('Gamma_beta_Csurfsurf_acc.pdf', dpi=600, format='pdf',  bbox_inches='tight')
plt.savefig('Gamma_beta_Csurfsurf_acc.png', dpi=600, format='png',  bbox_inches='tight')
plt.savefig('Gamma_beta_Csurfsurf_acc.svg', dpi=600, format='svg',  bbox_inches='tight')

plt.show()


slope, intercept, r_value, p_value, std_err = stats.linregress(Csurfsurf, predicted)
print("Prediction of Csurf,surf, r2: ")
print("Prediction of Csurf,surf, r2: ", r_value**2)
print("RMSE: ", rmse(predicted, Csurfsurf))
print("Median Abs Err: ", np.median(abs(predicted-Csurfsurf)))

# ----------------------------------------------------------------------------

plt.figure(figsize=(4.,3.5))

for i in range(0, len(Gamma)):
    plt.scatter(Gamma[i], beta[i], c=Csus[i],s=40, edgecolor='w', # s=300*a/Hmax
            marker=marker[i],cmap='viridis', vmin=0, vmax=1, zorder=9)

plt.xscale('log')
plt.yscale('log')
# -----

DG = 0.05*(np.max(Gamma)-np.min(Gamma))
bL = 0.05*(np.max(beta)-np.min(beta))

Gamma_pr = np.linspace(np.min(Gamma)-DG, np.max(Gamma)+DG, 1000)
beta_pr = np.linspace(np.min(beta)-bL, np.max(beta)+bL, 1000)
GG, bb = np.meshgrid(Gamma_pr,beta_pr)  

Gc = 38.0773685 # SOLVER optimised
bc = 1.613540437 # SOLVER optimised

C_sus_pred = np.exp(-GG/Gc)*np.exp(-bb/bc)

# ------

CS = plt.pcolormesh(GG,bb,C_sus_pred, cmap='viridis',shading='auto',
                    vmin=0, vmax=1)
CS.set_rasterized(True)

lvls = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
CS_lines2 = plt.contour(GG,bb,C_sus_pred, colors = 'w', linestyles = '--',
                    levels=lvls)
plt.clabel(CS_lines2, lvls[1::2], inline=True, fontsize=8, fmt='%1.1f')

plt.colorbar(CS, shrink=0.40, location='bottom',anchor=(1,0), label='$C_{\\beta}$ (-)')


plt.xlabel("$\\Gamma$ (-)")
plt.ylabel("$\\beta$ (-)")

plt.xlim([np.min(Gamma)-DG, np.max(Gamma)+DG])
plt.ylim([np.min(beta)-bL, np.max(beta)+bL])

plt.tight_layout()

plt.savefig('Gamma_beta_Csus_raster.pdf', dpi=600, format='pdf',  bbox_inches='tight')
plt.savefig('Gamma_beta_Csus_raster.png', dpi=600, format='png',  bbox_inches='tight')
plt.savefig('Gamma_beta_Csus_raster.svg', dpi=600, format='svg',  bbox_inches='tight')

plt.show()


# ----------------------------------------------------------------------------
# Subfigure -- internal

plt.figure(figsize=(2.5,2))
predicted = np.exp(-Gamma/Gc)*np.exp(-beta/bc)
 
for i in range(0, len(Gamma)):
    plt.scatter(Csus[i], predicted[i], marker=marker[i], edgecolor='b', facecolor='None')

plt.plot([0,1],[0,1], 'k-')
plt.plot([0.1,1],[0,0.9], 'k--')
plt.plot([0,0.9],[0.1,1], 'k--')
plt.plot([0.25,1],[0,0.75], 'k:')
plt.plot([0,0.75],[0.25,1], 'k:')

plt.xlabel("$C_{\\beta}$ (-),\n experimental")
plt.ylabel("$C_{\\beta}$ (-),\n modelled")

plt.xlim([0, 1])
plt.ylim([0, 1])
plt.tight_layout()

plt.savefig('Gamma_beta_Cbeta_acc.pdf', dpi=600, format='pdf',  bbox_inches='tight')
plt.savefig('Gamma_beta_Cbeta_acc.png', dpi=600, format='png',  bbox_inches='tight')
plt.savefig('Gamma_beta_Cbeta_acc.svg', dpi=600, format='svg',  bbox_inches='tight')

plt.show()


slope, intercept, r_value, p_value, std_err = stats.linregress(Csus, predicted)
print("Prediction of Cbeta, r2: ", r_value**2)
print("RMSE: ", rmse(predicted, Csus))
print("Median Abs Err: ", np.median(abs(predicted-Csus)))
