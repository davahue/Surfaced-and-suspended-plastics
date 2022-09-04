# -*- coding: utf-8 -*-
"""
Created on Sun May 16, 2021

@author: Daniel Valero (d.valero@un-ihe.org)
---------------------------
This script provides auxiliary functions for the ADV data analysis.

"""
# Import libraries ----------------
import glob
import numpy as np
from numpy import genfromtxt
import scipy as sp
from scipy.optimize import minimize
from scipy.stats.stats import pearsonr
#from scipy import signal
from scipy import stats
#from statsmodels import robust

# In-house libraries
import faux_filter as ff


# Basic filtering functions -------------------------------------------------------
# For more advanced, call: faux_filter

def getvalid(vx):
    """
    This function returns the percentage of valid points (not filtered)
    """

    return 100*np.count_nonzero(~np.isnan(vx))/len(vx)


def loadnSNRnCORfilter(filename, SNRthreshold, CORthreshold):
    """
    This function imports the basic ADV measured velocities
    filtered after a certain SNRthreshold and a CORthreshold
    """

    #data = pd.read_csv(filename) # pandas
    data = np.loadtxt(filename)
    data = data.T # get it in the correct form.

    # This column is the time
    t = data[1]
    # This column is vx
    vx = data[4]
    vy = data[5]
    vz1 = data[6]
    vz2 = data[7]

    # the 4 following rows are SNR
    SNRx = data[12]
    SNRy = data[13]
    SNRz1 = data[14]
    SNRz2 = data[15]

    # and the following 4, the CORR
    CORx = data[16]
    CORy = data[17]
    CORz1 = data[18]
    CORz2 = data[19]

    for i in range(0,len(vx)):
        # COR filtering:
        if CORx[i] < CORthreshold:
            vx[i] = np.nan
        if CORy[i] < CORthreshold:
            vy[i] = np.nan
        if CORz1[i] < CORthreshold:
            vz1[i] = np.nan
        if CORz2[i] < CORthreshold:
            vz2[i] = np.nan
        # Now, SNR:
        if SNRx[i] < SNRthreshold:
            vx[i] = np.nan
        if SNRy[i] < SNRthreshold:
            vy[i] = np.nan
        if SNRz1[i] < SNRthreshold:
            vz1[i] = np.nan
        if SNRz2[i] < SNRthreshold:
            vz2[i] = np.nan

    return t, vx,vy,vz1,vz2


# Other functions -----------------------------------------------------------

def read_info(info_route): # - function verified 17/05/2021
    """
    This function reads the info.txt file and returns the specific discharge
    (q) and the water level (Hmax).
    """
    prof_data = genfromtxt(info_route, delimiter='\t')
    # check order of info.txt
    #q = prof_data[1,0]/0.40 # unit discharge = Q/b, b = 0.40
    q_step = prof_data[1,0]
    Hmax= prof_data[1,3]
    ks= prof_data[1,4]
    
    return q_step,Hmax,ks
    

# Profile functions ---------------------------------------------------------


def get_T_scale(t,vx): # - function verified 01/01/2022
    """
    This function computes the integral time scale of a velocity series.
    Integration goes until first zero crossing.
    Input:
        - t: time array
        - vx: velocity array
    acor function performs a linear interpolation process for NaN
    replacement.
    """
    [s,acor] = ff.acor_func(t, vx)
    # i = 0; N = len(acor)
    zero_crossings = np.where(np.diff(np.sign(acor)))[0] # Credit goes to Jim Brissom: https://stackoverflow.com/questions/3843017/efficiently-detect-sign-changes-in-python
    ZC0 = zero_crossings[0] # First zero crossing
    Tx = sp.integrate.trapz(acor[0:ZC0],s[0:ZC0])

    return Tx


def get_mean_and_turb(t,vx,vy,vz1,acc_threshold):
    # acc_threshold = 60.0 # percent of accepted data required to consider a point
    acc_rate_x = getvalid(vx)
    acc_rate_y = getvalid(vy)
    acc_rate_z = getvalid(vz1)
    acc_rate = np.min([acc_rate_x, acc_rate_y, acc_rate_z]) # if any rate is below... then reject point
    
    if acc_rate < acc_threshold:
        print("Less than {} % data accepted, measurement discarded.\n".format(acc_threshold))
        # Mean velocities
        Vx_m = np.nan # np.nanmean(vx)
        Vy_m = np.nan # np.nanmean(vy)
        Vz_m = np.nan # np.nanmean(vz1)
        v1 = np.nan # vx-Vx_m[i]
        v2 = np.nan # vy-Vy_m[i]
        v3 = np.nan # vz1-Vz_m[i]
        # Shear stresses
        vxvy_s = np.nan # np.nanmean(v1*v2)
        vxvz_s = np.nan # np.nanmean(v1*v3)
        vyvz_s = np.nan # np.nanmean(v2*v3)    
        # Normal fluctuations
        vxvx_n = np.nan # np.nanmean(v1*v1)
        vyvy_n = np.nan # np.nanmean(v2*v2)
        vzvz_n = np.nan # np.nanmean(v3*v3)
        # Integral timescales
        Tx = np.nan
        Ty = np.nan
        Tz = np.nan
        
    else:
        # Mean velocities
        Vx_m = np.nanmean(vx)
        Vy_m = np.nanmean(vy)
        Vz_m = np.nanmean(vz1)
        v1 = vx-Vx_m
        v2 = vy-Vy_m
        v3 = vz1-Vz_m
        # Shear stresses
        vxvy_s = np.nanmean(v1*v2)
        vxvz_s = np.nanmean(v1*v3)
        vyvz_s = np.nanmean(v2*v3)    
        # Normal fluctuations
        vxvx_n = np.nanmean(v1*v1)
        vyvy_n = np.nanmean(v2*v2)
        vzvz_n = np.nanmean(v3*v3)
        # Integral timescales
        Tx = get_T_scale(t,v1)
        Ty = get_T_scale(t,v2)
        Tz = get_T_scale(t,v3)
            

    return Vx_m, Vy_m, Vz_m, vxvy_s, vxvz_s, vyvz_s, vxvx_n, vyvy_n, vzvz_n, Tx, Ty, Tz, acc_rate


def ADV_profile(subfol,SNRthreshold,CORthreshold,filt='GN2002W2003DV'): # - function verified 17/05/2021
    """
    This function analyses an entire ADV profile provided that it is stored 
    in a single folder. There should be an 'info.txt' file located in this
    folder, following the structure:
        q	Dz	Hmax
        59	0.02	0.29
    This function returns the profile mean and fluctuating quantities.
    
    Filtering (filt) options:
        'GN2002W2003DV': Goring and Nikora 2002, as modified by Wahl (2003) 
        and including my own suggestions to improve robustness.
        'ROCfilt': robust outlier cutoff, based on robust variance estimation
        and the universal threshold.

    """
    
    # Vectrino files:
    valid_fmt = '*.vna'
    files = glob.glob(subfol+valid_fmt)
    
    prof_data = genfromtxt(subfol+'info.txt', delimiter='\t')
    # check order of info.txt
    #q = prof_data[1,0]/0.40 # unit discharge = Q/b, b = 0.40
    Zmin= prof_data[1,1]+prof_data[1,4]/2 # consider Zmin as Z from roughness crests
                                          # and add half the roughness height.
                                          # we consider that the virtual origin of data
                                          # is at 1/2 of ks (within the roughness).
    Dz= prof_data[1,2]
    
    if Zmin < -0.0030: # in such case, we take the second measurement,
                       # according to the users manual, the pulse covers 
                       # 3-15 mm vertically, which overlaps with bed roughness
        imin = 1
        Zmin = Dz + Zmin
    else:
        imin = 0
      
    N = len(files)-imin
    
    #Zmin = 0.01 # m, from channel bed.
    Zmax = Zmin+(N-1)*Dz

        
    Z = np.linspace(Zmin,Zmax,N)
    
    # Initialise:
    Vx_m = np.zeros(N); vxvx_n = np.zeros(N)
    Vy_m = np.zeros(N); vyvy_n = np.zeros(N)
    Vz_m = np.zeros(N); vzvz_n = np.zeros(N)
    vxvy_s = np.zeros(N); vxvz_s = np.zeros(N); vyvz_s = np.zeros(N)
    Tx = np.zeros(N); Ty = np.zeros(N); Tz = np.zeros(N)
    
    acc_rate = np.zeros(N)
    
    for i in range(N):
        k = i+imin # vectors start at 0, don't forget it.
        filename = files[k] # one file at a time
        [t, vx,vy,vz1,vz2] = loadnSNRnCORfilter(filename, SNRthreshold, CORthreshold)
        sign = np.sign(np.nanmean(vx)) # ADV looking upstream?
        vx = sign*vx
        if filt=='GN2002W2003DV':
            vx = ff.GN2002W2003_DV(vx)
            vy = ff.GN2002W2003_DV(vy)
            vz1 = ff.GN2002W2003_DV(vz1)
            vz2 = ff.GN2002W2003_DV(vz2)
        elif filt=='ROC':
            [vx,_,] = ff.ROCfilt(vx)
            [vy,_,] = ff.ROCfilt(vy)
            [vz1,_,] = ff.ROCfilt(vz1)
            [vz2,_,] = ff.ROCfilt(vz2)
            
        # Acceptance rate
        acc_threshold = 60 # at least, 60 % of data should have been accepted 
                           # if we want a meaningful estimation.
        [Vx_m[i], Vy_m[i], Vz_m[i],
         vxvy_s[i], vxvz_s[i], vyvz_s[i],
         vxvx_n[i], vyvy_n[i], vzvz_n[i],
         Tx[i], Ty[i], Tz[i], 
         acc_rate[i]] = get_mean_and_turb(t,vx,vy,vz1,acc_threshold)
        
    return (Z,Vx_m,Vy_m,Vz_m,
            vxvx_n,vyvy_n,vzvz_n,
            vxvy_s,vxvz_s,vyvz_s,
            acc_rate)


def ADV_profile_rot(subfol,SNRthreshold,CORthreshold,filt='GN2002W2003DV'): # - function verified 17/05/2021
    """
    This function is identical to: ADV_profile, but it includes a rotation
    correction across the vertical profile equal to the median of 180*np.arctan2(Vym,Vxm)/np.pi.
    
    See documentation of ADV_profile() for further information.
    
    """
    print("First pass, unrotated data... \n")
    [Z,Vxm,Vym,Vzm,vxvxn,vyvyn,vzvzn,vxvy,vxvz,vyvz, #q,ufs,delta,Hmax, # delta,ushear,
             acc_rate] = ADV_profile(subfol,SNRthreshold,CORthreshold,filt='GN2002W2003DV')
    
    rot = np.arctan2(Vym,Vxm)
    rot_med = np.nanmedian(rot)
    
    print("Median rotation angle (deg): ", 180*rot_med/np.pi)
    print("Reacquiring data and projecting back rotation angle... \n")
    del Z,Vxm,Vym,Vzm,vxvxn,vyvyn,vzvzn,vxvy,vxvz,vyvz, acc_rate
    
    
    # Vectrino files:
    valid_fmt = '*.vna'
    files = glob.glob(subfol+valid_fmt)
    
    prof_data = genfromtxt(subfol+'info.txt', delimiter='\t')
    # check order of info.txt
    #q = prof_data[1,0]/0.40 # unit discharge = Q/b, b = 0.40
    Zmin= prof_data[1,1]+prof_data[1,4]/2 # consider Zmin as Z from roughness crests
                                          # and add half the roughness height.
                                          # we consider that the virtual origin of data
                                          # is at 1/2 of ks (within the roughness).
    Dz= prof_data[1,2]
    
    if Zmin < -0.0030: # in such case, we take the second measurement,
                       # according to the users manual, the pulse covers 
                       # 3-15 mm vertically, which overlaps with bed roughness
        imin = 1
        Zmin = Dz + Zmin
    else:
        imin = 0
      
    N = len(files)-imin
    
    #Zmin = 0.01 # m, from channel bed.
    Zmax = Zmin+(N-1)*Dz

        
    Z = np.linspace(Zmin,Zmax,N)
    
    # Initialise:
    Vx_m = np.zeros(N); vxvx_n = np.zeros(N)
    Vy_m = np.zeros(N); vyvy_n = np.zeros(N)
    Vz_m = np.zeros(N); vzvz_n = np.zeros(N)
    vxvy_s = np.zeros(N); vxvz_s = np.zeros(N); vyvz_s = np.zeros(N)
    Tx = np.zeros(N); Ty = np.zeros(N); Tz = np.zeros(N)
    
    acc_rate = np.zeros(N)
    
    for i in range(N):
        k = i+imin # vectors start at 0, don't forget it.
        filename = files[k] # one file at a time
        [t, vx,vy,vz1,vz2] = loadnSNRnCORfilter(filename, SNRthreshold, CORthreshold)
        sign = np.sign(np.nanmean(vx)) # ADV looking upstream?
        vx = sign*vx
        
        # Rotation to align with X-direction
        vxrot = vx*np.cos(-rot_med) - vy*np.sin(-rot_med)
        vyrot = vx*np.sin(-rot_med) + vy*np.cos(-rot_med)
        # overwrite vx, vy
        vx = vxrot
        vy = vyrot
        
        if filt=='GN2002W2003DV':
            vx = ff.GN2002W2003_DV(vx)
            vy = ff.GN2002W2003_DV(vy)
            vz1 = ff.GN2002W2003_DV(vz1)
            vz2 = ff.GN2002W2003_DV(vz2)
        elif filt=='ROC':
            [vx,_,] = ff.ROCfilt(vx)
            [vy,_,] = ff.ROCfilt(vy)
            [vz1,_,] = ff.ROCfilt(vz1)
            [vz2,_,] = ff.ROCfilt(vz2)
            
        # Acceptance rate
        acc_threshold = 60 # at least, 60 % of data should have been accepted 
                           # if we want a meaningful estimation.
        [Vx_m[i], Vy_m[i], Vz_m[i],
         vxvy_s[i], vxvz_s[i], vyvz_s[i],
         vxvx_n[i], vyvy_n[i], vzvz_n[i],
         Tx[i], Ty[i], Tz[i],
         acc_rate[i]] = get_mean_and_turb(t,vx,vy,vz1,acc_threshold)
        
    return (Z,Vx_m,Vy_m,Vz_m,
            vxvx_n,vyvy_n,vzvz_n,
            vxvy_s,vxvz_s,vyvz_s,
            Tx, Ty, Tz,
            acc_rate,rot)


def q_estimate(ushear, delta, ufs, Hmax, ks, PI):
    """
    This function integrates the discharge supported by the fitted log-wake law
    profile and considers the potential flow in the free stream region.
    """
    N = 300
    z = np.linspace(0.0001, Hmax, N)
    U = np.zeros(N)
    U[z<delta] = logwakelaw(z[z<delta], delta, ks, PI, ushear)
    U[z>delta] = ufs
    
    q = sp.integrate.trapz(U,z)
    
    return q


def get_BL_thickness(Z,VX_m): # - function verified 17/05/2021
    """
    This function inspects a mean velocity profile and returns the boundary
    layer thickness (99 % of free stream velocity) and the free stream
    velocity.
    """
    # Obtain related flow variables:        
    N = len(Z)
    ufs = np.nanmax(VX_m) # max non-nan value
    
    dZ_dVXm = np.gradient(Z)/np.gradient(VX_m, edge_order=2)
    dZ2dVXm2 = np.gradient(dZ_dVXm, edge_order=2)
    
    delta = Z[-1] # If no clear Z is detected, assume it is at free surface
    for i in range(N):
        if dZ2dVXm2[i] < 0 and VX_m[i] > 0.95*ufs: # 0.95 clause, in case of
                                                  # zig-zag err in lower measurements
            delta = Z[i]
            break
    # Quick plot -------------------
    # plt.figure()
    # plt.title('Data mean vel gradient')
    # plt.scatter(dZ_dVXm,Z/delta)
    # plt.xlabel('$\Delta Z / \Delta \\overline{ V_x}$')
    # plt.show()
    # -------------------------------
    return ufs,delta


def get_ushear_dzdu(Z, Vxm, delta, method='robust'): # - function verified 17/05/2021 
    """
    This function uses the slope of the mean velocity profile to obtain the
    shear velocity. The approach is to compute the inverse of the mean
    velocity gradient, which holds a linear slope. That linear slope is
    obtained through the correlation.
    
    delta here is a based on the dz/du gradient, so it is consistent for the
    purpose.
    
    method = 'robust': uses robust estimators, with breakthrough 50 %, to 
    estimate the slope.
    
    method = 'lstsq': uses np.polyfit, and polynominal order 1, to estimate 
    the slope.
    
    """
    if type('robust') != str:
        print('method is not a string. Using \'robust\' method instead.')
        method='robust' 
    elif not(method in ['Pearson', 'lstsq', 'robust']):
        print('Chosen method is not in the available options. Using \'robust\' method instead.')
        method='robust'
    
    kappa = 0.41
    # Cut the region inside the BL
    z = Z[Z<delta]
    u = Vxm[Z<delta]
    # Let's delete the last point, to avoid points too close to BL edge
    #z = z[:-1]
    #u = u[:-1]
    
    # plt.scatter(Vxm,Z,c='b') # profile
    # plt.scatter(u,z,c='r') # profile up to BL thickness
    # print(len(u))
    dz_du = np.gradient(z)/np.gradient(u, edge_order=2)
    if method=='Pearson':
    # 1. Through Pearson correlation:
        [r, p] = pearsonr(dz_du,z)
        m = r*np.std(dz_du)/np.std(z) # As suggested in Valero and Bung (2018, JHE)
        ushear = kappa/m
    
    if method=='lstsq':
        # 2. Through lsqm:
        m, c = np.polyfit(z[:-1],dz_du[:-1], 1) # I do not fully trust the end point with edge order = 1, by default
        ushear = kappa/m
    
    if method=='robust':
        # Alternative, robust-based estimation (50 % breadown point):
            [medslope, medintercept] = stats.siegelslopes(dz_du[:-1], z[:-1])
            ushear = kappa/medslope
    
    return ushear


def get_ushear_u1u2(Z, Vx_m, delta, ks):
    """
    This function returns the shear velocity based on a two points
    approximation of the log law with displacement thickness.
    The displacement thickness is approximated as 1/2 of ks.
    
    By default, "get_ushear_dzdu" should be more accurate.
    """
    kappa = 0.41
    # dr = 1/2*ks # displacement thickness, approximation.
    # dr was already introduced in Z as virtual origin of the data within the
    # function: ADV_profile()
    
    
    ZZ = (Z-delta)**2 # squared so min is at Z == delta
    index = np.argmin(ZZ)-1
    
    aux1 = kappa*(Vx_m[index] - Vx_m[0])
    a3 = Z[index] # - dr
    a4 = Z[0] # - dr
    aux2 = np.log(a3/a4)
    ushear = aux1/aux2
    
    return ushear


def logwakelaw(z, delta, ks, PI, ushear): # - function verified 02/05/2021
    """
    This function returns the mean velocity at a given depth based on the
    log-law mean velocity profile, as described in Eq. 1 of:
        
        Castro-Orgaz, O. (2010). Velocity Profile and Flow Resistance Models
        for Developing Chute Flow. J. Hydraul. Eng., 136, 447-452.
    
    It takes as an input:
        z: level position at which the log-law is evaluated.
        delta: boundary layer thickness.
        ks: sand roughness.
        ufs: free stream velocity.
        ushear: shear velocity.
        
    This function returns:
        Um: streamwise mean velocity.
        
    If needed, a Re dependant expression for B is provided in Eq. 3:39 of:
        
        Dey, S. (2014). Fluvial hydrodynamics. Berlin: Springer.
    
    I think this wake function is of Granville (1976):
        
        Granville PS (1976) A modified law of the wake for turbulent shear
        layers. J Fluids Eng 98:578–580
        
    -------------------------------------------------------------------------
    Function was verified against Fig. 2 of Castro-Orgaz (2010), using
    parameters:   
        ufs = 0.30
        delta=0.20
        ks = 1.5E-3
        z = np.linspace(0.001, delta, 1000)
        Cf = get_Cf(delta, ks)
        ushear = get_ushear(ufs, Cf)
        Um = logwakelaw(z, delta, ks, ufs, ushear)
        plt.semilogy((ufs-Um)/ushear, z/delta)
    Which cross-verifies:
        get_Cf
        get_ushear
    
    """
    # Parameters
    ka = 0.41 # von karman coefficient
    B = 8.5 #8.5 # log-law constant for rough flow.
    # PI = 0.0 #0.2 wake constant, calibrated by Castro-Orgaz 2010 in dev BL.

    # terms
    term1 = (1./ka)*np.log(z/ks) # naperian log (LN).
    term2 = (1./ka)*(1+6*PI)*(z/delta)**2
    term3 = -(1./ka)*(1+4*PI)*(z/delta)**3

    # All together:
    Um = ushear*(term1 + B + term2 + term3)

    return Um

def loglaw(z, delta, ks, ushear): # - function verified 02/05/2021
    """
    This function returns the mean velocity at a given depth based on the
    log-law mean velocity profile, as described in Eq. 1 of:
        
        Castro-Orgaz, O. (2010). Velocity Profile and Flow Resistance Models
        for Developing Chute Flow. J. Hydraul. Eng., 136, 447-452.
    
    It takes as an input:
        z: level position at which the log-law is evaluated.
        delta: boundary layer thickness.
        ks: sand roughness.
        ufs: free stream velocity.
        ushear: shear velocity.
        
    This function returns:
        Um: streamwise mean velocity.
        
    If needed, a Re dependant expression for B is provided in Eq. 3:39 of:
        
        Dey, S. (2014). Fluvial hydrodynamics. Berlin: Springer.       

    """
    # Parameters
    ka = 0.41 # von karman coefficient
    B = 8.5 #8.5 # log-law constant for rough flow.
    

    # terms
    term1 = (1./ka)*np.log(z/ks) # neperian log (LN).


    # All together:
    Um = ushear*(term1 + B)

    return Um

def veldefect_logwakelaw(delta, PI, z): # - function verified 19/05/2021
    """
    This function implements Eq. 3 of:
        
        Castro-Orgaz, O. (2010). Velocity Profile and Flow Resistance Models
        for Developing Chute Flow. J. Hydraul. Eng., 136, 447-452.
    
    to compute the velocity defect.
    """
    # Parameters
    ka = 0.41 # von karman coefficient
    #B = 8.5 #8.5 # log-law constant for rough flow.
    #PI = 0.2 #0.2 # Wake constant, calibrated by Castro-Orgaz 2010 in dev BL.

    # terms
    term1 = -(1./ka)*np.log(z/delta) # neperian log (LN).
    term2 = (1./ka)*2*PI
    term3 = -(1./ka)*(1+6*PI)*(z/delta)**2 # negative here!
    term4 = (1./ka)*(1+4*PI)*(z/delta)**3
    
    defect = term1 + term2 + term3 + term4
    
    defect[z>delta] = np.nan
    # All together:
    #Um = ushear*(term1 + B + term2 + term3)

    return defect

def find_delta_logwakelaw(z,Vxm,ufs,ushear, PI): # - function verified 19/05/2021
    """
    This function minimizes the error between ADV data and velocity defect law
    (log-wake law) to obtain delta. Note that a PI value needs to be pressumed
    and that delta does not satisfy that u(z=delta) = ufs.
    
    In order to focus better the outer layer, the relative position inside the
    BL (z/delta) is used as a weight in the error function.
    
    Inluded management of nan.
    """

    N = 500
    # nan data can be a pain here. Let's create aux vectors without nan
    z = z[~np.isnan(Vxm)]
    Vxm = Vxm[~np.isnan(Vxm)]
    defect_data = (ufs-Vxm)/ushear
    
    delta_i = np.linspace(z.min(),z.max(),N)
    error = 100*np.ones(N)
    for i in range(1,N):
        defect_model = veldefect_logwakelaw(delta_i[i], PI, z)
        aa = defect_model[z<=delta_i[i]]
        bb = defect_data[z<=delta_i[i]]
        
        # w = z[z<=delta_i[i]]/delta_i[i] # we could weight data for the least-squares
        # error[i] = np.sum(w*(aa-bb)**2)/len(aa)
        error[i] = np.sum((aa-bb)**2)/len(aa)
    
    i_ = np.argmin(error)
    delta = delta_i[i_]
    
    print("error2 obtaining delta: ", error[i_])
    return delta


def find_ushear_delta_logwakelaw(z,Vxm,ufs,PI,ks): # - function verified 19/05/2021
    """
    This function minimizes the error between ADV data and velocity log-wake law
    to obtain ushear and delta. Note that a PI value needs to be pressumed
    and that delta does not satisfy that u(z=delta) = ufs.
    
    In order to focus better the outer layer, the relative position inside the
    BL (z/delta) is used as a weight in the error function.
    
    Inluded management of nan.
    """

    kappa = 0.41
    B = 8.5
    
    # N = 500
    # nan data can be a pain here. Let's create aux vectors without nan
    z = z[~np.isnan(Vxm)]
    Vxm = Vxm[~np.isnan(Vxm)]
    
    if z[0] < ks: # The first point may fall within the roughness viscous layer, so we can delete.
        z = z[1:]
        Vxm = Vxm[1:]
        
    
    """
    Based on the stack overflow solution:
    f = lambda x: np.sum((np.sqrt((x[0]-xi)**2+(x[1]-yi)**2)-d)**2)
    res = optimize.minimize(f, (initial_x, initial_y))
    see: https://stackoverflow.com/questions/27682547/how-to-write-the-scipy-optimize-minimizes-parameter
    
    and Scipy documentation: https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.optimize.minimize.html
    
    """
    
    N = len(Vxm)
    # x[0] = ushear
    # x[1] = delta
    f = lambda x: 1/N*np.sum((Vxm - x[0]*(1/kappa* np.log(z/ks) + B + 1/kappa*(1+6*PI)*(z/x[1])**2 - 1/kappa*(1+4*PI)*(z/x[1])**3))**2)
        # np.sum((np.sqrt((x[0]-xi)**2+(x[1]-yi)**2)-d)**2)
    res = minimize(f, (0.05*ufs, np.nanmax(z)/2)) # f, (initial_x, initial_y)
    
    [ushear, delta] = res.x
   
    
    return ushear, delta, res


def get_ushear_uxuy(Z, vxvz, delta): # - function verified 17/05/2021 
    """
    This function uses the slope of the Reynolds shear stresses x-y to
    obtain the shear velocity. The slope intersection at Z = 0 give us
    the value of ushear.
    We use a robust slope detector (median-based, breakdown point = 50 %)
    
    Z: vertical coordinate
    vxvz: Reynolds shear stresses at x-z
    delta: boundary layer thickness
    """
    [medslope, medintercept] = \
    stats.siegelslopes(-vxvz[Z<delta], Z[Z<delta])    
    ushear = np.sqrt(medintercept) # where the line intersects Z = 0 -> u^*2

    return ushear
