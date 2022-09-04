# -*- coding: utf-8 -*-
"""
Created on Tue Jan  4 18:30:37 2022

Rouse - related functions
@author: Daniel Valero
"""

import numpy as np
import scipy.interpolate as interpolate
from scipy.optimize import minimize
import scipy.stats as stats
from matplotlib import pyplot as plt


# Functions ------------------------------------------------------------------

def Rouse(z,a,Ca,Hmax,beta):
    """Rouse profile, as defined in Rouse 1961 Eq. 245."""
    eta0 = Hmax - a
    eta = z - a
    C = np.zeros(len(eta))
    
    for i in range(len(eta)):
        if eta[i] < 0: # for y < a, bed load transport.
            C[i] = np.nan
        else:            
            C[i] = Ca*((1-eta[i]/eta0)/(1+eta[i]/a))**beta
    
    return C

def invRouse(C,a,Ca,Hmax,beta):
    """
    This function receives a given probability level: C/Ca
    and returns the eta coordinate where this happens.
    Not very useful, but not wrong.
    
    Rouse distribution properties:
        a,Ca,Hmax,beta
    
    Function verified 12.12.2021 against Rouse(z,a,Ca,Hmax,beta)
    """
    # Exemplary parameters
    # Hmax = 0.278 # m, from experiment
    # a = 0.04
    # Ca = 1
    # beta = 1
    eta0 = Hmax-a
    
    # Solve:
    fun = lambda x: np.abs(C/Ca-((1-x/eta0)/(1+x/a))**beta)
    #bnds = ((0, Hmax-a))
    res = minimize(fun, Hmax/2, method='L-BFGS-B')
    z_sol = res.x[0] + a # z = eta + a
    
    return z_sol

def inv_sampling_Rouse(a,Ca,Hmax,beta,nsample):
    """
    This function uses the CDF of the Rouse profile [Rouse()] to perform 
    inverse transform sampling using a uniform distribution between [0,1] for
    the CDF.
    
    nsample: number of samples to be generated.
    
    returns: Zinv, which is the Z coordinates of random plastic samples 
    travelling in suspension following the Rouse(_,a,Ca,Hmax,beta) profile.
    
    """
    Z = np.linspace(a, Hmax, 1000)
    C = Rouse(Z,a,Ca,Hmax,beta)
    CS = np.cumsum(C)
    CDF = CS/np.max(CS)
    CDF[0] = 0 # we enforce we start at zero. Error \sim 1/len(Z)
    CDF[-1] = 1 # last point too, at 1.0.
    # Create a function to obtain the Z coordinate given a CDF value.
    f = interpolate.interp1d(CDF, Z, kind='linear')
    # now, uniform distribution to enable inverse transform sampling:
    uniform_rvs = stats.uniform.rvs(loc=0, scale=1, size=nsample)
    
    Zinv = f(uniform_rvs) # random coordinates of plastic samples flyig in suspension following the Rouse(_,a,Ca,Hmax,beta) profile.
    return Zinv


# Data analysis ---------------------------------------------------------------

def extract_particle_stats(data, plane):
    """
    data: csv-type dataset, in pandas format, including:
        FrameId, x, y, z data.
    plane: x position at which stats are extracted
    """
    
    t = np.asarray(data["FrameId"])*1/60 # 60 Hz
    x = np.asarray(data['x']); z = np.asarray(data['y']); y = np.asarray(data['z']) # beware of y-z camera convention
    try:
        surf = np.asarray(data['Surfaced'])        
    except Exception:
        print("No surfaced variable found in data file. We set all to 0!")
        surf = np.zeros(len(x), dtype=int)
    
    vx = np.gradient(x, t)/100; vy = np.gradient(y, t)/100; vz = np.gradient(z, t)/100 # cm to m.
    # when we cross x = midplane:
    ZC = np.where(np.diff(np.sign(x-plane)))[0] # From Jim Brisson in Stack Overflow: https://stackoverflow.com/questions/3843017/efficiently-detect-sign-changes-in-python
    zc = []
    for i in range(0, len(ZC)):
        if x[ZC[i]+1] - x[ZC[i]-1] > 0: # we crossed in the right direction
            zc.append(ZC[i]+1) # +1, so we factually cross the "plane"
        # else, we moved backwards
        
        
    ZC = np.asarray(zc)        
    
    return t[ZC], x[ZC], y[ZC], z[ZC], vx[ZC], vy[ZC], vz[ZC], surf[ZC]


def get_Ca(z, a, Hmax, beta, Csus):
    """
    Selects Ca in the Rouse profile in order to satisfy mass conservation (i,e.,
    same samples in experiment and predicted by Rouse).
    """
    Ca_list = np.arange(1,int(Csus)+1)
    error = np.ones(len(Ca_list))
    
    if a < 0:
        a = 1 # if a is detected as out of the free surface (possible in KS), we enforce a=0 (free surface)
    
    for i in range(0,len(Ca_list)):
        C = Rouse(-z, a, Ca_list[i], Hmax, beta)
        error[i] = Csus-np.nansum(C)
        
    # When error hits zero, same samples predicted by Rouse and the experiment.
    Err_zero_i = np.where(np.diff(np.sign(error)))[0]
    if len(Err_zero_i) == 0:
        Err_zero_i = np.nanargmin(error)
    
    try:
        Ca = Ca_list[Err_zero_i][0]
    except Exception:
        Ca = Ca_list[Err_zero_i] # there was only a value
    # except Exception:
    #     Ca = Ca_list[Err_zero_i]    
    
    return Ca


# Load dataset ---------------------------------------------------------------

def get_flow_sample_chars(filename):
    """
    Selects the corresponding shear velocity for the flow discharge and rising
    velocity for the sample corresponding to each experiment.
    """
    
    # Define experiment paremeters -----------------------------------------------
    ushear_list = [0.0213042562736631,
                   0.0320206390412612,
                   0.0443926526227695,
                   0.0535055637118659,
                   0.0612913658282144] # ushear from log-wake law fit to ADV data
    
    if "V1" in filename:
        ushear = ushear_list[0]
    elif "V2" in filename:
        ushear = ushear_list[1]
    elif "V3" in filename:
        ushear = ushear_list[2]
    elif "V4" in filename:
        ushear = ushear_list[3]
    elif "V5" in filename:
        ushear = ushear_list[4]
    else:
        ushear = np.nan
        print("Flow case not defined in file name: ", filename)
    
    w_list = [26.2,
              29.2,
              18.6, 
              11.7,
              5.4,
              2.4,
              101.8] # mm/s, converted to m/s in return call
    
    if "Cup_PP_100" in filename:
        w = w_list[0]
    elif "Cup_PP_98_def" in filename:
        w = w_list[1]
    elif "Cup_PP_50" in filename:
        w = w_list[2]
    elif "Cup_PP_05" in filename:
        w = w_list[3]
    elif "Film_HDPE_100" in filename:
        w = w_list[4]    
    elif "Film_HDPE_15" in filename:
        w = w_list[5]
    elif "Mask" in filename:
        w = w_list[6]    
    else:
        w = np.nan
        print("Sample class not defined in file name: ", filename)
    
    return ushear, w/1000

def get_sample_geomchars(filename):
    """
    Selects the corresponding sample geometrical properties.
    """
    
    lmax_list = [83, # mm
                 88,
                 83,
                 36,
                 133,
                 45,
                 171] # 
    lortho_list = [69, # mm
                   72,
                   69,
                   32,
                   100,
                   44,
                   94] # 
    
    if "Cup_PP_100" in filename:
        lmax = lmax_list[0]
        lortho = lortho_list[0]
    elif "Cup_PP_98_def" in filename:
        lmax = lmax_list[1]
        lortho = lortho_list[1]
    elif "Cup_PP_50" in filename:
        lmax = lmax_list[2]
        lortho = lortho_list[2]
    elif "Cup_PP_05" in filename:
        lmax = lmax_list[3]
        lortho = lortho_list[3]
    elif "Film_HDPE_100" in filename:
        lmax = lmax_list[4]   
        lortho = lortho_list[4]
    elif "Film_HDPE_15" in filename:
        lmax = lmax_list[5]
        lortho = lortho_list[5]
    elif "Mask" in filename:
        lmax = lmax_list[6] 
        lortho = lortho_list[6]
    else:
        lmax = np.nan
        print("Sample class not defined in file name: ", filename)
    
    return lmax, lortho

def get_label(filename):
    
    if "_d5" in filename:
        label_plot = filename.replace('_d5','')
    else:
        label_plot = filename
        
    return label_plot 

# ----------------------------------------------------------------------------

def plot_significance(Nsamples, st, pv, zpsorted, a_surfaced_video, a_surfaced_KS, i_,
                      BedLev, WaterLev, filename):
    """
    Quick plot showing how KS stat and p-value change when adding samples 
    (ordered from bottom to free surface). Two plots considered:
        
        Plot 1: KS and p-value with increasing no samples with increasing "z"
            value for "a".
        Plot 2: KS and p-value vs vertical column.
        
    """
    
    
    """
    # Quick plot to verify that sample 1 and sample 2 are comparable (for stats.ks_2samp(-) )
    index = 5
    ai = -zpsorted[index]-0.0001
    ZRouse_true = -cp.inv_sampling_Rouse(ai,1,Hmax,beta,10000)
    
    plt.hist(zpsorted[0:index],bins=50,density=True);
    plt.hist(ZRouse_true, bins=200, histtype="step",density=True)
    
    """
    # Plot 1 - KS and p-value with increasing no samples with increasing "z" value for "a".---
    plt.figure(figsize=(3.5,3.5))
    plt.loglog(Nsamples,st, marker='o', c='k', label='KS stat')
    plt.loglog(Nsamples,pv, marker='x', c='b', label='p-value')
    
    Nsamples_surf = np.where(np.diff(np.sign(zpsorted-a_surfaced_video)))[0]
    
    plt.axvline(x=Nsamples_surf[-1], color='r', linestyle='-.', label='Max. depth \nof surface particle')
    plt.axvline(x=Nsamples[i_], color='r', lw=3, linestyle=':', label='End of suspension transport')
    
    N = len(zpsorted)
    
    label_lvl = "$z - H$ = " + str(round(zpsorted[int(0.10*N)],2)) + " cm"
    plt.axvline(x=10, color='grey', linestyle='--', label=label_lvl)
    label_lvl = "$z - H$ =  = " + str(round(zpsorted[int(0.35*N)],2)) + " cm"
    plt.axvline(x=50, color='grey', linestyle='--', label=label_lvl)
    label_lvl = "$z - H$ =  = " + str(round(zpsorted[int(0.75*N)],2)) + " cm"
    plt.axvline(x=100, color='grey', linestyle='--', label=label_lvl)
    # label_lvl = "$z - H$ =  = " + str(round(zpsorted[140],2)) + " cm"
    # plt.axvline(x=140, color='grey', linestyle='--', label=label_lvl)
    
    
    plt.xlabel("Number of samples")
    plt.ylabel("KS / p-value")
    plt.ylim([1E-3, 1.0])
    
    plt.tight_layout()
    #plt.legend(loc="best")
    
    # Save plot --------------------------------------------------------------   
    plotname = filename[:-4]+'_KStest_Nsamples'
    plt.savefig(plotname+'.pdf', dpi=600, format='pdf',  bbox_inches='tight')
    plt.savefig(plotname+'.png', dpi=600, format='png',  bbox_inches='tight')
    plt.savefig(plotname+'.svg', dpi=600, format='svg',  bbox_inches='tight')
    # ------------------------------------------------------------------------
    
    # Plot 2 - KS and p-value vs vertical column ---
    plt.figure(figsize=(3.5,3.5))
    plt.semilogx(st, zpsorted, marker='o', c='k', label='KS stat')
    plt.semilogx(pv, zpsorted, marker='x', c='b', label='p-value')
    
    plt.axhline(y=0.0, color='b', linestyle='-', label='Free surface')
    plt.axhline(y=a_surfaced_video, color='r', linestyle='-.', label='Max. depth \nof surface particle')
    plt.axhline(y=a_surfaced_KS, color='r', lw=3, linestyle=':', label='End of suspension transport')
    plt.axhline(y=-BedLev+WaterLev, color='k', linestyle='--', label='Roughness crests')
    
    plt.ylabel("$z - H$ (cm)")
    plt.xlabel("KS / p-value")
    plt.xlim([1E-3, 1.0])
    
    plt.tight_layout()
    #plt.legend(loc="best")
    
    # Save plot --------------------------------------------------------------
    plotname = filename[:-4]+'_KStest_z-column'
    plt.savefig(plotname+'.pdf', dpi=600, format='pdf',  bbox_inches='tight')
    plt.savefig(plotname+'.png', dpi=600, format='png',  bbox_inches='tight')
    plt.savefig(plotname+'.svg', dpi=600, format='svg',  bbox_inches='tight')
    # ------------------------------------------------------------------------
    return
