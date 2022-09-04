# -*- coding: utf-8 -*-
"""
Created on Mon Jan 10 19:21:55 2022

This code produces a figures showing differences between suspended transport profiles

@author: D. Valero
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors
#import matplotlib.gridspec as gridspec

def get_label(filename):
    """Remove d5 from the label."""
    if "_d5" in filename:
        label_plot = filename.replace('_d5','')
    else:
        label_plot = filename
        
    return label_plot 


# FROM PREVIOUS CODES....
kmax_list = np.loadtxt('kmax_list11330.txt')
st_list = np.loadtxt('st_list11330.txt')
pv_list = np.loadtxt('pv_list11330.txt')


# Main plot ------------------------------------------------------------------

fig = plt.figure(constrained_layout=True, figsize=(7.5,8.5))


gs = fig.add_gridspec(2, 3)
f_ax_st = fig.add_subplot(gs[0, :]) # TOP
f_ax_pv = fig.add_subplot(gs[1, :]) # BOTTOM

Ntotal,NT = np.shape(kmax_list)

y_space = 10**np.linspace(-2,1,100)
CB_st = f_ax_st.hist2d(kmax_list.flatten(),st_list.flatten(),
                       bins=(kmax_list[0,:], y_space), cmap=plt.cm.binary,
                       norm=colors.LogNorm())    
fig.colorbar(CB_st[3], shrink=0.2, location='top',anchor=(1,1), label='count(-)',ax=f_ax_st) # anchor=(1,1.1), orientation='horizontal',

f_ax_st.set_rasterized(True)

f_ax_st.loglog()
f_ax_st.plot(kmax_list[0,:], np.median(st_list,axis=0), c='r', ls='--', label='Median (true suspended)')
f_ax_st.plot(kmax_list[0,:], np.percentile(st_list,q=5,axis=0), c='r', ls=':', lw=2,label='5$^{th}$/95$^{th}$ perc. (true suspended)')
f_ax_st.plot(kmax_list[0,:], np.percentile(st_list,q=95,axis=0), c='r', ls=':',lw=2)

f_ax_st.axes.xaxis.set_ticklabels([])

f_ax_st.text(6.6,9,'A',weight="bold")

f_ax_st.set_ylabel('Kolmogorov-Smirnov statistic\n [KS stat] (-)')    
f_ax_st.legend(loc='upper left', ncol=2)
# --------------
# BOTTOM

CB_pv = f_ax_pv.hist2d(kmax_list.flatten(),pv_list.flatten(),
                       bins=(kmax_list[0,:], y_space), norm=colors.LogNorm(),
                       cmap=plt.cm.binary) # cmap=plt.cm.Blues)
fig.colorbar(CB_pv[3], shrink=0.2, location='bottom',anchor=(1,0), label='count(-)',ax=f_ax_pv) # anchor=(1,1.1), orientation='horizontal',

f_ax_pv.set_rasterized(True)

f_ax_pv.loglog()
f_ax_pv.plot(kmax_list[0,:], np.median(pv_list,axis=0), c='r', ls='--', label='Median (true suspended)')
f_ax_pv.plot(kmax_list[0,:], np.percentile(pv_list,q=5,axis=0), c='r', ls=':', lw=2,label='5$^{th}$/95$^{th}$ perc. (true suspended)')
f_ax_pv.plot(kmax_list[0,:], np.percentile(pv_list,q=95,axis=0), c='r', ls=':', lw=2)

f_ax_pv.set_xlabel('Number of samples')
f_ax_pv.set_ylabel('p-value (-)')
f_ax_pv.set_ylim([1E-4, 1E1])


# --------------

# Import data ----------------------------------------------------------------


filelist = ['Cup_PP_100_V1_d5',  'Cup_PP_100_V2_d5', 'Cup_PP_100_V3_d5',
            'Cup_PP_100_V4_d5', 'Cup_PP_100_V5_d5', # end
            'Cup_PP_98_def_V2_d5', 'Cup_PP_98_def_V3_d5', 'Cup_PP_98_def_V4_d5',
            'Cup_PP_98_def_V5_d5',
            'Cup_PP_50_V1_d5', 'Cup_PP_50_V2_d5', 
            'Cup_PP_50_V4_d5', 'Cup_PP_50_V5_d5', # end
            'Cup_PP_05_V1_d5', 'Cup_PP_05_V2_d5', 'Cup_PP_05_V3_d5', 
            'Cup_PP_05_V4_d5', 'Cup_PP_05_V5_d5',
            'Film_HDPE_100_V1_d5', 'Film_HDPE_100_V2_d5', 
            'Film_HDPE_100_V4_d5', 'Film_HDPE_100_V5_d5',
            'Film_HDPE_15_V1_d5', 'Film_HDPE_15_V2_d5', 'Film_HDPE_15_V3_d5', 
            'Film_HDPE_15_V4_d5', 'Film_HDPE_15_V5_d5',
            'Mask_V1_d5', 'Mask_V2_d5', 'Mask_V3_d5', 'Mask_V4_d5', 'Mask_V5_d5'
            ]


# Ink the lines
i = 0
colors = plt.cm.PuOr(np.linspace(0,1,len(filelist))) # tab20
#colors = plt.cm.cividis_r(np.linspace(0,1,len(filelist))) # tab20
colors2 = plt.cm.Blues_r(np.linspace(0,1,len(filelist))) # tab20

for filename in filelist:
    
    path = filename
    Nsamples = np.loadtxt(path+'_Nsamples'+'.txt')
    pv = np.loadtxt(path+'_pv'+'.txt')
    st = np.loadtxt(path+'_st'+'.txt')
    zp = np.loadtxt(path+'_zp'+'.txt')
    surfp = np.loadtxt(path+'_surfp'+'.txt')
    [plane, beta, binsize, Ca, a_surfaced_video, a_surfaced_KS] = np.loadtxt(path+'_scalars'+'.txt')
    
    # Into the plot:
    if np.sum(1-surfp) > 15: # over 50 samples in suspension
        label_plot = get_label(filename)
        if "HDPE" in filename:
            f_ax_st.plot(Nsamples, st, c=colors2[i], ls='-', alpha=1) #, label='KS stat') 'k'# marker='+',
            f_ax_pv.plot(Nsamples, pv, c=colors2[i], ls='-', alpha=1, label=label_plot) #, label='p-value') # 'deeppink'
            i = i+1
        else:
            f_ax_st.plot(Nsamples, st, c=colors[i], ls='-', alpha=1) #, label='KS stat') 'k'# marker='+',
            f_ax_pv.plot(Nsamples, pv, c=colors[i], ls='-', alpha=1, label=label_plot) #, label='p-value') # 'deeppink'
            i = i+1

# ----------------------------------------------------------------------------

#f_ax_pv.legend(loc='lower right', ncol=2) # Legend me.

# Plot settings
plt.tight_layout(h_pad=2.25)

f_ax_pv.text(6.6,4,'B',weight="bold")

plotname = 'KS_surface_transport_detection'
plt.savefig(plotname+'.pdf', dpi=600, format='pdf',  bbox_inches='tight')
plt.savefig(plotname+'.png', dpi=600, format='png',  bbox_inches='tight')
plt.savefig(plotname+'.svg', dpi=600, format='svg',  bbox_inches='tight')
