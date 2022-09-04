# -*- coding: cp1252 -*-
"""
Daniel Valero, 04.11.2021

This library brings together several basic despiking-filtering
functions, based on statistical methods.

ROC and DVED are explained in:
    Valero, D., Chanson, H., & Bung, D. B. (2020).
Robust estimators for free surface turbulence characterization:
A stepped spillway application. Flow Measurement and
Instrumentation, 76, 101809.

ROC and DVED are however simplifications of the method of
Goring & Nikora (2002), after the modifications suggested by
Wahl (2003):

    Goring, D. G., & Nikora, V. I. (2002). Despiking acoustic Doppler
velocimeter data. Journal of hydraulic engineering, 128(1), 117-126.

    Wahl, T. L. (2003). Discussion of "Despiking acoustic doppler
velocimeter data" by Derek G. Goring and Vladimir I. Nikora.
Journal of Hydraulic Engineering, 129(6), 484-487.

These two papers merit reference when using ROC or DVED, but
I suggest to clarify that ROC and DVED are simplified versions.
Note that ROC is the simplest of the filters.

The main reasoning behind these statistical methods is that "good
data tend to cluster together". There are a couple of strong
assumptions:
- Normality (when using the universal threshold), and thus symmetry.
- Stationarity of the data (the mean and variance does not change
during the data recording).

"""

import numpy as np
import scipy as sp
from scipy import signal as sp_signal
# from statsmodels import robust

from matplotlib import pyplot as plt
# from matplotlib import cm
# from mpl_toolkits.mplot3d import Axes3D


def ROCfilt(raw):
    """
    This function performs: Robust Outlier Cutoff (ROC) filtering.
    Also, it returns the percentaje of outliers.
    """
    signal = np.copy(raw)
    mur = np.nanmedian(signal)
    kn = 1.483    # for normal distributions!!

    # robust.mad() is not prepared for np.nan!
    # Thus, we need to: robust.mad(a[~np.isnan(a)])
    
    # sigmar = robust.mad(signal[~np.isnan(signal)]) # care, in stats.robust they included kn as part of MAD
    sigmar = kn*np.nanmedian(abs(signal - np.nanmedian(signal)))
    # Universal threshold
    # Ns = len(signal)
    Nu = np.count_nonzero(~np.isnan(signal))
    lambu = np.sqrt(2*np.log(Nu))

    Ns = len(signal)

    Vmax = lambu*sigmar # the maximum possible for this number of outliers...
    
    outliers = 0
    
    for i in range(0, Ns):

        # nan is neither < Lb nor > Ub, so it will not
        # be accounted again as an outlier.

        if np.abs(signal[i] - mur) > Vmax:

            signal[i] = np.nan
            outliers += 1
    

    return signal, float(outliers)/Ns


def DVEDfilt(raw):
    """
    This function performs: Depth-Velocity Elliptical Despiking (D-VED) filtering.
    Also, it returns the percentaje of outliers.
    """
    signal = np.copy(raw)
    mur = np.nanmedian(signal)
    eta = signal - mur
    vgrad = np.gradient(eta)
    
    kn = 1.483    # for normal distributions!!
    # eta_sigmar = robust.mad(eta[~np.isnan(eta)]) # Care, this already includes kn
    eta_sigmar = kn*np.nanmedian(abs(eta - np.nanmedian(eta)))
    
    # v_sigmar = robust.mad(vgrad[~np.isnan(vgrad)])
    v_sigmar = kn*np.nanmedian(abs(vgrad - np.nanmedian(vgrad)))
    
    # Universal threshold
    # Ns = len(signal)
    Nu = np.count_nonzero(~np.isnan(signal))
    lambu = np.sqrt(2*np.log(Nu))

    Ns = len(signal)

    eta_max = lambu*eta_sigmar # the maximum possible for this number of outliers...
    v_max = lambu*v_sigmar

    rx2 = eta_max**2.0
    ry2 = v_max**2.0
    outliers = 0
    for i in range(0, Ns):

        # nan is neither < Lb nor > Ub, so it will not
        # be accounted again as an outlier.
        
        a2 = eta[i]**2.0
        b2 = vgrad[i]**2.0

        if (a2 / rx2 + b2 / ry2) > 1.0: # falls outside of the ellipse

            signal[i] = np.nan
            outliers += 1           

    return signal, float(outliers)/Ns


def GN2002W2003_DV(U):
    """
    This function implements the Goring and Nikora (2002) method
    as modified by Wahl (2003), with a few additions from my end (Daniel
    Valero).
    It performs no outliers replacement, but fills in nan instead.
    
    Visually verified.
    """
    # Initialize
    Ufilt = np.copy(U)
    Umu = np.nanmedian(Ufilt) # the data could already come with nan values
    kn = 1.483    # for normal distributions  
    
    #robust.mad(Ufilt[~np.isnan(Ufilt)])
    
    # robust.mad() is not prepared for np.nan! Thus, we need to:
    # robust.mad(a[~np.isnan(a)])    
    # robust.mad includes the kn coefficient!! beware!!

    # Despiking threshold ---------------------------------------------------
    NU = np.count_nonzero(~np.isnan(Ufilt)) # NU as independent number of
        # samples measured.
    # Universal threshold:
    # lambU = np.sqrt(2*np.log(Nu)) # The Universal Threshold, suggested
    # by Goring and Nikora, yet Wahl 2003 suggested using erf() as per
    # Chauvenet´s criterion:
    lambU = np.sqrt(2)*sp.special.erfinv(1 - 1./(2*NU))

    # Create vectors for vel difs and difdifs -------------------------------
    Ns = len(U)
    U_ = np.copy(U) - Umu # U fluctuations
    Usigma = kn*np.nanmedian(abs(U_ - np.nanmedian(U_)))
    # Before we start.... what if we have np.nan() values.
    # We should consider them.
    
    dU_ = np.gradient(U_,edge_order=2)
    ddU_ = np.gradient(dU_,edge_order=2)
    dUsigma = kn*np.nanmedian(abs(dU_ - np.nanmedian(dU_)))
    # robust.mad(dU_[~np.isnan(dU_)])
    ddUsigma = kn*np.nanmedian(abs(ddU_ - np.nanmedian(ddU_)))
    
    # robust.mad(ddU_[~np.isnan(ddU_)])

    # visz3D(U_,dU_,ddU_)
    """
    sns.histplot(x=vx, y=ddvx,bins=50,cbar=True)
    plt.show()
    sns.histplot(x=dvx, y=ddvx,bins=50,cbar=True)
    plt.show()
    """
    new = U_+ddU_
    UddUcor = np.corrcoef(U_[~np.isnan(new)],ddU_[~np.isnan(new)])[1,0]
    m = UddUcor*ddUsigma/Usigma
    """
    Goring and Nikora recommended 
    m2 =np.nansum(U_*ddU_)/np.nansum(U_*U_)
    for the slope. However, when doing U_*U_ or U_*ddU_ outliers have a big 
    role and in the U_*U_ amplify further than in U_*ddU_, which means that
    we tend to detect a smaller slope of the ellipsoid.
    When doing: UddUcor*ddUsigma/Usigma, sigma variables have been obtained
    with robust estimators and so are less sensitive.
    
    """
    theta = np.arctan2(m,1) # rotation angle between U_ and ddU_
    # print("Angle U-ddU", np.rad2deg(theta))
    # Now, it is a matter of bringing down (rotate) the points U_, ddU_
    # bz -theta (theta is the angle from axis to points, we want to move points
    # back to axis).

    U_rot = U_*np.cos(-theta) - ddU_*np.sin(-theta)
    ddU_rot = U_*np.sin(-theta) + ddU_*np.cos(-theta)
    # visz3D(U_rot,dU_,ddU_rot)
    # Goring and Nikora + Wahl use the Usigma, dUsigma, ddUsigma times lambU
    # as axes of the ellipsoid. We have now rotated the points, and should
    # obtain the new axes instead. Otherwise, our ellipsoid uses axes based on
    # unrotated data.
    Urotsigma = kn*np.nanmedian(abs(U_rot - np.nanmedian(U_rot)))
    # robust.mad(U_rot[~np.isnan(U_rot)]) # This already includes kn!! check documentation robust.MAD, it is 1/c
    ddUrotsigma = kn*np.nanmedian(abs(ddU_rot - np.nanmedian(ddU_rot)))
    # robust.mad(ddU_rot[~np.isnan(ddU_rot)])
    

    
    """Ellipsoid axes are clear now: Urotsigma, dUsigma, ddUrotsigma * lambU"""
    a = (U_rot/Urotsigma/lambU)**2   
    b = (dU_/dUsigma/lambU)**2
    c = (ddU_rot/ddUrotsigma/lambU)**2    
    INSIDE = np.zeros(Ns)
        
    # visz3D(a,b,c)
 
    # print("a-axis max", np.nanmax(a))
    # print("b-axis max", np.nanmax(b))
    # print("c-axis max", np.nanmax(c))
    # print("a+b+c max", np.nanmax(a+b+c))
    for i in range(0,Ns): # 0, to number of samples
        
        if a[i]+b[i]+c[i] > 1:
            Ufilt[i] = np.nan  # this point is out of the ellipsoid.
            INSIDE[i] = False
        else:
            INSIDE[i] = True 
    
    # Number of outliers detected
    print("GN2002W2003_DV, Number of outliers detected: ", np.count_nonzero((INSIDE==False)))
    
    # Quick visualization
    # [xe,ye,ze] = generate_ellipsoid(Urotsigma*lambU,
    #                             dUsigma*lambU,
    #                             ddUrotsigma*lambU)
    # xeunrot = xe*np.cos(theta) - ze*np.sin(theta)
    # zeunrot = xe*np.sin(theta) + ze*np.cos(theta)
    """
    This was for rotating:
    U_rot = U_*np.cos(-theta) - ddU_*np.sin(-theta)
    ddU_rot = U_*np.sin(-theta) + ddU_*np.cos(-theta)
    But we do now unrotation (so, -theta).
    """
    # visz3D_theta_INSIDE_ellipsoid(U_,dU_,ddU_,theta,INSIDE,xeunrot,ye,zeunrot)
    
    return Ufilt


def visz3D(U_,dU_,ddU_):
    """Quick visualization of the velocity points and their derivatives."""
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(U_,dU_,ddU_, c='k', marker='o',alpha=0.5)
    ax.set_xlabel('$U$ (m/s)');ax.set_ylabel('$\Delta U$ (m/s)')
    ax.set_zlabel('$\Delta^2 U$ (m/s)')
    plt.show()
    
    return

def visz3D_INSIDE(U_,dU_,ddU_,I):
    """Quick visualization of the velocity points and their derivatives.
    I stands for the vector INSIDE"""
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(U_[I==True],dU_[I==True],ddU_[I==True], c='k', marker='o',
               alpha=0.05,label = 'Points accepted')
    ax.scatter(U_[I==False],dU_[I==False],ddU_[I==False], c='r', marker='o',
               alpha=0.5,label = 'Points filtered out')
    ax.set_xlabel('$U$ (m/s)');ax.set_ylabel('$\Delta U$ (m/s)')
    ax.set_zlabel('$\Delta^2 U$ (m/s)')
    plt.legend(loc='best')
    plt.show()
    
    return


def visz3D_theta(U_,dU_,ddU_,theta):
    """Quick visualization of the velocity points and their derivatives. 
    Includes the axis of the ellipsoid."""
    # Axis of the ellipsoid
    xx = np.linspace(np.nanmin(U_),np.nanmax(U_),100)
    zz = xx*np.sin(theta)
    yy = np.zeros(100)
    
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(U_,dU_,ddU_, c='k', marker='o',alpha=0.5)
    ax.plot3D(xx,yy,zz, 'k')
        
    ax.set_xlabel('$U$ (m/s)');ax.set_ylabel('$\Delta U$ (m/s)')
    ax.set_zlabel('$\Delta^2 U$ (m/s)')
    plt.show()
    
    return

def visz3D_theta_INSIDE_ellipsoid(U_,dU_,ddU_,theta,I,xeunrot,ye,zeunrot):
    """Quick visualization of the velocity points and their derivatives. 
    Includes the axis of the ellipsoid."""
    # Axis of the ellipsoid
    xx = np.linspace(np.nanmin(U_),np.nanmax(U_),100)
    zz = xx*np.sin(theta)
    yy = np.zeros(100)
    
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(U_[I==True],dU_[I==True],ddU_[I==True], c='k', marker='o',
               alpha=0.15,label = 'Points accepted')
    ax.scatter(U_[I==False],dU_[I==False],ddU_[I==False], c='r', marker='o',
               alpha=0.5,label = 'Points filtered out')
    # Ellipsoid
    ax.plot3D(xx,yy,zz, c='b', label='axis-1 ellipsoid')
    surf = ax.plot_surface(xeunrot, ye, zeunrot, rstride=1, cstride=1, color='k',
                    alpha=0.10,  label='Ellipsoid') # cmap=cm.viridis
    surf._facecolors2d = surf._facecolor3d
    surf._edgecolors2d = surf._edgecolor3d
    """"
    Past two lines are required as per a matplotlib bug:
        https://stackoverflow.com/questions/54994600/pyplot-legend-poly3dcollection-object-has-no-attribute-edgecolors2d
    """
        
    # Labels
    ax.set_xlabel('$U$ (m/s)')
    ax.set_ylabel('$\Delta U$ (m/s)')
    ax.set_zlabel('$\Delta^2 U$ (m/s)')
    
    plt.legend(loc='best')
    plt.show()
    
    return


def generate_ellipsoid(a1,b1,c1):
    """
    This is taken from:
        https://stackoverflow.com/questions/7819498/plotting-ellipsoid-with-matplotlib
    """
    
    coefs = (a1**2, b1**2, c1**2)  # Coefficients: (a1*x)**2 + (b1*y)**2 + (c1*z)**2 = 1 
    # Radii corresponding to the coefficients:
    rx, ry, rz = np.sqrt(coefs)#1/np.sqrt(coefs)
    
    # Set of all spherical angles:
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    
    # Cartesian coordinates that correspond to the spherical angles:
    # (this is the equation of an ellipsoid):
    x = rx * np.outer(np.cos(u), np.sin(v))
    y = ry * np.outer(np.sin(u), np.sin(v))
    z = rz * np.outer(np.ones_like(u), np.cos(v))

    return x,y,z


# Interpolation ----------------------------------------------


def fill_linear_nan(signal):
    """
    This functions finds np.nan values and replaces them by the
    linear interpolation using the neighbouring values.

    However, it may happen that all the values that enter are
    non valid (np.nan) and so, it should be robust to this situation.

    It has been here decided that, when more than 50 % of the data
    is nan, no interpolation is done and thus the estimators later
    compute will also return nan.
    """

    Nu = np.count_nonzero(~np.isnan(signal))
    N = len(signal)

    if float(Nu)/N > 0.50:
        

        xi = np.arange(len(signal))
        mask = np.isfinite(signal)
    
        signal_filt = np.interp(xi, xi[mask], signal[mask])

        return signal_filt

    else:

        return signal


# Other useful time series tools -------------------------------------------------

def xcor_half(ya, yb):
    """
    Dimensionless cross-correlation. Half is returned to match
    the dimensions of ya, yb.
    """

    yax = (ya - np.mean(ya))/(np.std(ya)*len(ya))
    ybx = (yb - np.mean(yb))/np.std(yb)

    corr = sp_signal.correlate(ybx, yax, mode='full')

    return corr[int(corr.size/2):]

def xcov_half(ya, yb):
    """
    Dimensional cross-correlation. Half is returned to match
    the dimensions of ya, yb.
    """

    #yax = (ya - np.mean(ya))/(np.std(ya)*len(ya))
    #ybx = (yb - np.mean(yb))/np.std(yb)

    cov = sp_signal.correlate(yb, ya, mode='full')

    return cov[int(cov.size/2):]

def acor_func(t, ya):
    """
    Autocorrelation function, provided:
        - t: time vector
        - ya: a time series over t
    """
    
    ya_interp = fill_linear_nan(ya) # nan values are not most friendly here
    acor = xcor_half(ya_interp, ya_interp)
    
    dt = np.median( np.gradient(t)) # may there be some shorter times due to bottom check or similar, or even at startup.
    s = dt*np.linspace(0, len(acor), len(acor))
    
    return s, acor


