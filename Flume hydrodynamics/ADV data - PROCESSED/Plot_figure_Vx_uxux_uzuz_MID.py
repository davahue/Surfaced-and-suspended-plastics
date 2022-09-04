# -*- coding: utf-8 -*-
"""
Created on Sun May 16, 2021

@author: Daniel Valero (daniel.valero@kit.edu)
---------------------------
This script produces a figure with mean velocity and normal turbulent fluctiations.

Cross-section: middle

"""
# Import libraries ----------------
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

# Physical parameters
Nu = 1E-6


def logwakelaw(z, delta, ks, PI, ushear): # - function verified 02/05/2021
    """
    This function is originally included in faux_input.py, copied here to 
    reduce dependencies.
    --
    
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


# Import data V1 -------------------------------------------------------------
# Mean vel data
Z1 = np.loadtxt('Middle//Middle_V1_Z.txt')
Vxm1 = np.loadtxt('Middle//Middle_V1_Vxm.txt')
[ushear_logwakelaw1,ushear_uxuz1] = np.loadtxt('Middle//Middle_V1_ushear.txt')
[delta_logwakelaw1,delta_uxuz1] = np.loadtxt('Middle//Middle_V1_delta.txt')
[q1,Hmax1,ks1,ufs1,PI1] = np.loadtxt('Middle//Middle_V1_q_Hmax_ks_ufs_PI.txt')
# uxuxn
vxvxn1 = np.loadtxt('Middle//Middle_V1_vxvxn.txt')
vzvzn1 = np.loadtxt('Middle//Middle_V1_vzvzn.txt')

# Import data V2 -------------------------------------------------------------
# Mean vel data
Z2 = np.loadtxt('Middle//Middle_V2_Z.txt')
Vxm2 = np.loadtxt('Middle//Middle_V2_Vxm.txt')
[ushear_logwakelaw2,ushear_uxuz2] = np.loadtxt('Middle//Middle_V2_ushear.txt')
[delta_logwakelaw2,delta_uxuz2] = np.loadtxt('Middle//Middle_V2_delta.txt')
[q2,Hmax2,ks2,ufs2,PI2] = np.loadtxt('Middle//Middle_V2_q_Hmax_ks_ufs_PI.txt')
# uxuxn
vxvxn2 = np.loadtxt('Middle//Middle_V2_vxvxn.txt')
vzvzn2 = np.loadtxt('Middle//Middle_V2_vzvzn.txt')

# Import data V3 -------------------------------------------------------------
# Mean vel data
Z3 = np.loadtxt('Middle//Middle_V3_Z.txt')
Vxm3 = np.loadtxt('Middle//Middle_V3_Vxm.txt')
[ushear_logwakelaw3,ushear_uxuz3] = np.loadtxt('Middle//Middle_V3_ushear.txt')
[delta_logwakelaw3,delta_uxuz3] = np.loadtxt('Middle//Middle_V3_delta.txt')
[q3,Hmax3,ks3,ufs3,PI3] = np.loadtxt('Middle//Middle_V3_q_Hmax_ks_ufs_PI.txt')
# uxuxn
vxvxn3 = np.loadtxt('Middle//Middle_V3_vxvxn.txt')
vzvzn3 = np.loadtxt('Middle//Middle_V3_vzvzn.txt')

# Import data V4 -------------------------------------------------------------
# Mean vel data
Z4 = np.loadtxt('Middle//Middle_V4_Z.txt')
Vxm4 = np.loadtxt('Middle//Middle_V4_Vxm.txt')
[ushear_logwakelaw4,ushear_uxuz4] = np.loadtxt('Middle//Middle_V4_ushear.txt')
[delta_logwakelaw4,delta_uxuz4] = np.loadtxt('Middle//Middle_V4_delta.txt')
[q4,Hmax4,ks4,ufs4,PI4] = np.loadtxt('Middle//Middle_V4_q_Hmax_ks_ufs_PI.txt')
# uxuxn
vxvxn4 = np.loadtxt('Middle//Middle_V4_vxvxn.txt')
vzvzn4 = np.loadtxt('Middle//Middle_V4_vzvzn.txt')

# Import data V5 -------------------------------------------------------------
# Mean vel data
Z5 = np.loadtxt('Middle//Middle_V5_Z.txt')
Vxm5 = np.loadtxt('Middle//Middle_V5_Vxm.txt')
[ushear_logwakelaw5,ushear_uxuz5] = np.loadtxt('Middle//Middle_V5_ushear.txt')
[delta_logwakelaw5,delta_uxuz5] = np.loadtxt('Middle//Middle_V5_delta.txt')
[q5,Hmax5,ks5,ufs5,PI5] = np.loadtxt('Middle//Middle_V5_q_Hmax_ks_ufs_PI.txt')
# uxuxn
vxvxn5 = np.loadtxt('Middle//Middle_V5_vxvxn.txt')
vzvzn5 = np.loadtxt('Middle//Middle_V5_vzvzn.txt')

# ----------------------------------------------------------------------------
# Parsing data (2): Nezu and Rodi (1986)
# black circles
filename_N_bO = '_liter-data stream fluctuations\\1986 Nezu ALL.xlsx'
xl = pd.ExcelFile(filename_N_bO)
df = xl.parse(xl.sheet_names[0])
data_N_bO = np.asarray(df)
# ----------------------------------------------------------------------------
# Cameron et al (2017)-JFM
# Parsing data (4): Cameron et al (2017)
# circles
filename_Cam_O = '_liter-data stream fluctuations\\2017 Cameron uu circles.xlsx'
xl = pd.ExcelFile(filename_Cam_O)
df = xl.parse(xl.sheet_names[0])
data_Cam_O = np.asarray(df)
# diamonds
filename_Cam_D = '_liter-data stream fluctuations\\2017 Cameron uu rombos.xlsx'
xl = pd.ExcelFile(filename_Cam_D)
df = xl.parse(xl.sheet_names[0])
data_Cam_D = np.asarray(df)
# squares
filename_Cam_s = '_liter-data stream fluctuations\\2017 Cameron uu squares.xlsx'
xl = pd.ExcelFile(filename_Cam_s)
df = xl.parse(xl.sheet_names[0])
data_Cam_s = np.asarray(df)
# triangles1
filename_Cam_t1 = '_liter-data stream fluctuations\\2017 Cameron uu triangles1.xlsx'
xl = pd.ExcelFile(filename_Cam_t1)
df = xl.parse(xl.sheet_names[0])
data_Cam_t1 = np.asarray(df)
# triangles2
filename_Cam_t2 = '_liter-data stream fluctuations\\2017 Cameron uu triangles2.xlsx'
xl = pd.ExcelFile(filename_Cam_t2)
df = xl.parse(xl.sheet_names[0])
data_Cam_t2 = np.asarray(df)


# ----------------------------------------------------------------------------
label1= 'V1' #'$u_{fs} = $'+str(round(100*ufs1)/100)+' m/s'
label2= 'V2' #'$u_{fs} = $'+str(round(100*ufs2)/100)+' m/s'
label3= 'V3' #'$u_{fs} = $'+str(round(100*ufs3)/100)+' m/s'
label4= 'V4' #'$u_{fs} = $'+str(round(100*ufs4)/100)+' m/s'
label5= 'V5' #'$u_{fs} = $'+str(round(100*ufs5)/100)+' m/s'


# Figure ---------------------------------------------------------------------

fig, axs = plt.subplots(1, 3, figsize=(6.75, 3), sharey=False)
# Mean velocities and log-wake law -------------------------------------------
# V1
axs[0].scatter(Vxm1/ushear_logwakelaw1,Z1*ushear_logwakelaw1/Nu, marker='<',
               facecolors='1', edgecolors='k',s=20, label=label1)
z1 = np.linspace(0.001,delta_logwakelaw1,1000)
Ulog1 = logwakelaw(z1, delta_logwakelaw1, ks1, PI1, ushear_logwakelaw1)
axs[0].plot(Ulog1/ushear_logwakelaw1,z1*ushear_logwakelaw1/Nu,c='g',lw=1, label='Log-wake law')
# V2
axs[0].scatter(5+Vxm2/ushear_logwakelaw2,Z2*ushear_logwakelaw2/Nu, marker='s',
               facecolors='0.75', edgecolors='k',s=20, label=label2)
z2 = np.linspace(0.0025,delta_logwakelaw2,1000)
Ulog2 = logwakelaw(z2, delta_logwakelaw2, ks2, PI2, ushear_logwakelaw2)
axs[0].plot(5+Ulog2/ushear_logwakelaw2,z2*ushear_logwakelaw2/Nu,c='g',lw=1)
# V3
axs[0].scatter(10+Vxm3/ushear_logwakelaw3,Z3*ushear_logwakelaw3/Nu, marker='P',
               facecolors='w', edgecolors='k',s=20, label=label3)
z3 = np.linspace(0.005,delta_logwakelaw3,1000)
Ulog3 = logwakelaw(z3, delta_logwakelaw3, ks3, PI3, ushear_logwakelaw3)
axs[0].plot(10+Ulog3/ushear_logwakelaw3,z3*ushear_logwakelaw3/Nu,c='g',lw=1)
# V4
axs[0].scatter(15+Vxm4/ushear_logwakelaw4,Z4*ushear_logwakelaw4/Nu, marker='o',
               facecolors='0.25', edgecolors='k',s=20, label=label4)
z4 = np.linspace(0.0075,delta_logwakelaw4,1000)
Ulog4 = logwakelaw(z4, delta_logwakelaw4, ks4, PI4, ushear_logwakelaw4)
axs[0].plot(15+Ulog4/ushear_logwakelaw4,z4*ushear_logwakelaw4/Nu,c='g',lw=1)
# V5
axs[0].scatter(20+Vxm5/ushear_logwakelaw5,Z5*ushear_logwakelaw5/Nu, marker='>',
               facecolors='0', edgecolors='k',s=20, label=label5)
z5 = np.linspace(0.010,delta_logwakelaw5,1000)
Ulog5 = logwakelaw(z5, delta_logwakelaw5, ks5, PI5, ushear_logwakelaw5)
axs[0].plot(20+Ulog5/ushear_logwakelaw5,z5*ushear_logwakelaw5/Nu,c='g',lw=1)

# Formatting figure
axs[0].set_xlabel('$u^+$ (-)')
axs[0].set_ylabel('$z^+$ (-)')
axs[0].text(4, 10000,"A",verticalalignment='center',weight='bold')

axs[0].legend(loc='best',bbox_to_anchor=(1.025, 1.45), fontsize='small', ncol=2)

axs[0].semilogy()
# ----------------------------------------------------------------------------
# vxvx -----------------------------------------------------------------------
# Nezu and Rodi (1986)
y_del = np.linspace(0.10, 1.5, 200)
Ku = 2.26 # According to Nezu and Rodi (1986)
Lambdu = 0.88
u_NR =Ku*np.exp(-Lambdu*y_del)

axs[1].scatter(data_N_bO[:,1]**2, data_N_bO[:,0], facecolor='orange', edgecolor='k',
               alpha=0.40, marker='.', label='Nezu and Rodi (1986) data')
axs[1].plot(u_NR**2, y_del, c='darkorange', ls='-.', lw=2, label='Modelled data of Nezu and Rodi (1986)')

# For Cameron, Stuart and Nikora 2017-JFM data

CC = 'green'
axs[1].scatter(data_Cam_O[:,1], data_Cam_O[:,0], facecolor=CC, edgecolor='k',
               alpha=0.40, marker='x', label='Cameron et al. (2017) data', zorder=4)
axs[1].scatter(data_Cam_D[:,1], data_Cam_D[:,0], facecolor=CC, edgecolor='k',
               alpha=0.40, marker='x', zorder=4)
axs[1].scatter(data_Cam_s[:,1], data_Cam_s[:,0], facecolor=CC, edgecolor='k',
               alpha=0.40, marker='x', zorder=4)
axs[1].scatter(data_Cam_t1[:,1], data_Cam_t1[:,0], facecolor=CC, edgecolor='k',
               alpha=0.40, marker='x', zorder=4)
axs[1].scatter(data_Cam_t2[:,1], data_Cam_t2[:,0], facecolor=CC, edgecolor='k',
               alpha=0.40, marker='x', zorder=4)



u_NS14 = 2.222*np.exp(-0.8367*y_del) # Formula from thesis of Nezu  (1977)
u_NS_min = 2.181*np.exp(-0.8773*y_del) # min and max, from Valero 2018 reanalysis.
u_NS_max = 2.264*np.exp(-0.7961*y_del)
axs[1].plot(u_NS14**2, y_del, c='b', ls='--', lw=2, label='Modelled data of Cameron et al. (2017)', zorder=9)
axs[1].plot(u_NS_min**2, y_del, c='b', ls=':', lw=2, label='5$-$95% uncertainty - Cameron et al. (2017)', zorder=9)
axs[1].plot(u_NS_max**2, y_del, c='b', ls=':', lw=2, zorder=9) # , label='95% uncertainty - Cameron et al. (2017)'
#axs[1].fill_between(y=y_del, x1=u_NS_min**2, x2=u_NS_max**2, color='g', alpha=0.25, label='95 % confidence, Cameron et al. (2017) data')

# V1
axs[1].scatter(vxvxn1/ushear_logwakelaw1**2,Z1/delta_logwakelaw1, marker='<',
               facecolors='1', edgecolors='k',s=20)
# V2
axs[1].scatter(vxvxn2/ushear_logwakelaw2**2,Z2/delta_logwakelaw2, marker='s',
               facecolors='0.75', edgecolors='k',s=20)
# V3
axs[1].scatter(vxvxn3/ushear_logwakelaw3**2,Z3/delta_logwakelaw3, marker='P',
               facecolors='w', edgecolors='k',s=20)
# V4
axs[1].scatter(vxvxn4/ushear_logwakelaw4**2,Z4/delta_logwakelaw4, marker='o',
               facecolors='0.25', edgecolors='k',s=20)
# V5
axs[1].scatter(vxvxn5/ushear_logwakelaw5**2,Z5/delta_logwakelaw5, marker='>',
               facecolors='0', edgecolors='k',s=20)


# Formatting figure
axs[1].set_xlabel('$\overline{u_x\' u_x\'} / u_* ^2$ (-)')
axs[1].set_ylabel('$z/\\delta$ (-)')
axs[1].set_xlim(0,10)
#axs[1].set_xlim(0.1,10)
#axs[1].semilogx()
axs[1].text(2.5, 2.35,"B",verticalalignment='center',weight='bold')
axs[1].legend(loc='best',bbox_to_anchor=(2.15, 1.45), fontsize='small', ncol=2)
# ----------------------------------------------------------------------------
# vzvz -----------------------------------------------------------------------

Kv = 1.23 # According to Nezu and Rodi (1986)
Lambdv = 0.67
v_NR =Kv*np.exp(-Lambdv*y_del)
# Parsing Nezu
filename_N_bO = '_liter-data normal fluctuations\\1986 Nezu black circles.xlsx'
xl = pd.ExcelFile(filename_N_bO)
df = xl.parse(xl.sheet_names[0])
data_N_bO = np.asarray(df)
# black crosses
filename_N_bc = '_liter-data normal fluctuations\\1986 Nezu crosses.xlsx'
xl = pd.ExcelFile(filename_N_bc)
df = xl.parse(xl.sheet_names[0])
data_N_bc = np.asarray(df)
# black triangles
filename_N_t = '_liter-data normal fluctuations\\1986 Nezu triangles.xlsx'
xl = pd.ExcelFile(filename_N_t)
df = xl.parse(xl.sheet_names[0])
data_N_t = np.asarray(df)
# white circles
filename_N_wc = '_liter-data normal fluctuations\\1986 Nezu white circles.xlsx'
xl = pd.ExcelFile(filename_N_wc)
df = xl.parse(xl.sheet_names[0])
data_N_wc = np.asarray(df)

# Plot Nezu and Rodi (1986)

axs[2].scatter(data_N_bO[:,1]**2, data_N_bO[:,0], facecolor='orange', edgecolor='k',
               alpha=0.40, marker='.', label='Nezu and Rodi (1986)')
axs[2].scatter(data_N_bc[:,1]**2, data_N_bc[:,0], facecolor='orange', edgecolor='k',
               alpha=0.40, marker='.')
axs[2].scatter(data_N_t[:,1]**2, data_N_t[:,0], facecolor='orange', edgecolor='k',
               alpha=0.40, marker='.')
axs[2].scatter(data_N_wc[:,1]**2, data_N_wc[:,0], facecolor='orange', edgecolor='k',
               alpha=0.40, marker='.')
axs[2].plot(v_NR**2, y_del, c='darkorange', ls='-.', lw=2, label='Fit to Nezu and Rodi (1986)', zorder=9)





# For Cameron, Stuart and Nikora 2017-JFM data
# Parse data Cameron et al. (2017)
# circles
filename_Cam_O = '_liter-data normal fluctuations\\2017 Cameron ww circles.xlsx'
xl = pd.ExcelFile(filename_Cam_O)
df = xl.parse(xl.sheet_names[0])
data_Cam_O = np.asarray(df)
# diamonds
filename_Cam_D = '_liter-data normal fluctuations\\2017 Cameron ww rombos.xlsx'
xl = pd.ExcelFile(filename_Cam_D)
df = xl.parse(xl.sheet_names[0])
data_Cam_D = np.asarray(df)
# squares
filename_Cam_s = '_liter-data normal fluctuations\\2017 Cameron ww squares.xlsx'
xl = pd.ExcelFile(filename_Cam_s)
df = xl.parse(xl.sheet_names[0])
data_Cam_s = np.asarray(df)
# triangles1
filename_Cam_t1 = '_liter-data normal fluctuations\\2017 Cameron ww triangles1.xlsx'
xl = pd.ExcelFile(filename_Cam_t1)
df = xl.parse(xl.sheet_names[0])
data_Cam_t1 = np.asarray(df)
# triangles1
filename_Cam_t2 = '_liter-data normal fluctuations\\2017 Cameron ww triangles2.xlsx'
xl = pd.ExcelFile(filename_Cam_t2)
df = xl.parse(xl.sheet_names[0])
data_Cam_t2 = np.asarray(df)

# Now fits
v_NS14 = 1.108*np.exp(-0.6626*y_del) # Formula from thesis of Nezu (1977)
v_NS_min = 1.063*np.exp(-0.7469*y_del)
v_NS_max = 1.153*np.exp(-0.5783*y_del)

CC = 'green'

axs[2].scatter(data_Cam_O[:,1], data_Cam_O[:,0], facecolor=CC, edgecolor='k',
               alpha=0.40, marker='x', label='Cameron et al. (2017)', zorder=4)
axs[2].scatter(data_Cam_D[:,1], data_Cam_D[:,0], facecolor=CC, edgecolor='k',
               alpha=0.40, marker='x', zorder=4)
axs[2].scatter(data_Cam_s[:,1], data_Cam_s[:,0], facecolor=CC, edgecolor='k',
               alpha=0.40, marker='x', zorder=4)
axs[2].scatter(data_Cam_t1[:,1], data_Cam_t1[:,0], facecolor=CC, edgecolor='k',
               alpha=0.40, marker='x', zorder=4)
axs[2].scatter(data_Cam_t2[:,1], data_Cam_t2[:,0], facecolor=CC, edgecolor='k',
               alpha=0.40, marker='x', zorder=4)

axs[2].plot(v_NS14**2, y_del, c='b', ls='--',lw=2, zorder=9)
axs[2].plot(v_NS_min**2, y_del, c='b', ls=':', lw=2, zorder=9)
axs[2].plot(v_NS_max**2, y_del, c='b', ls=':', lw=2, zorder=9)

# Our data here...

# V1
axs[2].scatter(vzvzn1/ushear_logwakelaw1**2,Z1/delta_logwakelaw1, marker='<',
               facecolors='1', edgecolors='k',s=20)
# V2
axs[2].scatter(vzvzn2/ushear_logwakelaw2**2,Z2/delta_logwakelaw2, marker='s',
               facecolors='0.75', edgecolors='k',s=20)
# V3
axs[2].scatter(vzvzn3/ushear_logwakelaw3**2,Z3/delta_logwakelaw3, marker='P',
               facecolors='w', edgecolors='k',s=20)
# V4
axs[2].scatter(vzvzn4/ushear_logwakelaw4**2,Z4/delta_logwakelaw4, marker='o',
               facecolors='0.25', edgecolors='k',s=20)
# V5
axs[2].scatter(vzvzn5/ushear_logwakelaw5**2,Z5/delta_logwakelaw5, marker='>',
               facecolors='0', edgecolors='k',s=20)


# Formatting figure
axs[2].set_xlabel('$\overline{u_z\' u_z\'} / u_* ^2$ (-)')
axs[2].set_ylabel('$z/\\delta$ (-)')
axs[2].set_xlim(0,2)
#axs[2].set_xlim(0.1,10)
#axs[2].semilogx()

#axs[2].legend(loc='best',bbox_to_anchor=(1.25, 1.45), fontsize='small', ncol=2)
axs[2].text(0.55, 2.35,"C",verticalalignment='center',weight='bold')

fig.subplots_adjust(left = -0.4,wspace = 0.275)
#fig.tight_layout()

# ----------------------------------------------------------------------------
#plt.subplots_adjust(bottom=0.15)
plt.savefig('Vx_uxux_uzuz_middle.pdf', dpi=600, format='pdf',  bbox_inches='tight')
plt.savefig('Vx_uxux_uzuz_middle.svg', dpi=600, format='svg',  bbox_inches='tight')
plt.savefig('Vx_uxux_uzuz_middle.png', dpi=600, format='png',  bbox_inches='tight')
plt.show()