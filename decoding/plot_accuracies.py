#!/gpfs01/bartels/user/vplikat/anaconda3/bin/python
# -*- coding: utf-8 -*-

import os
import glob
from pathlib import Path
# import/export data
import h5py
import numpy as np   # most important numerical calculations
# plotting
import matplotlib.pyplot as plt

HOME            = str(Path.home())
PROJ_DIR        = os.path.join(HOME, 'Documents/Master_Thesis/DATA/MRI')
DERIVATIVES_DIR = os.path.join(PROJ_DIR, 'derivatives')
ANALYSIS = 'ROI-analysis'
DECODER = 'LDA'
RESULTS_DIR     = os.path.join(DERIVATIVES_DIR, 'decoding', 'object_decoding',
                               'train75_test25_notmagic', ANALYSIS, DECODER)
SUB = 'sub-01'
ROIS = [
        'V1', 'V2', 'V3', 'hV4', 
        'V3A', 'V3B', 
        'LO1', 'LO2', 
        'VO1', 'VO2', 
        'TO1', 'TO2', 
        'FEF', 'SPL1',
        'PHC1', 'PHC2',
        'IPS1' ,'IPS2', 'IPS3', 'IPS4','IPS5',
        'ventricle'
      ]

COLORS = ['r','b','y']
MARKERSTYLE = ['s','D','^', '>','*','<','o','v','+','.','h','H']

SMOOTH_SIZES = ['nativespace', '2mm-smoothed-nativespace','3mm-smoothed-nativespace']# 
SCALES = ['scale_None', 'scale_z']
CUTOFFS = ['cutoff_inf', 'cutoff_2.0']
TRANSFORMATIONS = ['feat_trans_None', 'feat_trans_PCA']

fig, (ax0, ax1) = plt.subplots(2,1,sharex=True)
fig.set_figwidth(30)
fig.set_figheight(30)

for s,smooth in enumerate(SMOOTH_SIZES):
    color = COLORS[s]
    idx=0
    for root, directories, files in os.walk(os.path.join(RESULTS_DIR, smooth)):
        if 'sub-01' in directories:
            linestyle = color+MARKERSTYLE[idx]
            idx += 1
            path_as_list = root.split(os.sep)
            index = path_as_list.index(smooth)+1
            
            p_values = []
            accuracies = []
            for roi in ROIS:
                roi_file = os.path.join(root,'sub-01', roi + '.hdf5')
                res = h5py.File(roi_file, 'r')
                accuracies.append(res['accuracy'][()]-1/3)
                p_values.append(res['p_value'][()])
            ax0.plot(accuracies,linestyle,markersize=12, label='_'.join(path_as_list[index:]))#, label=scale+cutoff+transform)
            ax0.plot(accuracies, ls=':')
            ax0.tick_params(axis='both',labelsize=20)
            
            ax1.plot(p_values, linestyle,markersize=12, label='_'.join(path_as_list[index:]))#,label=scale+cutoff+transform)
            ax1.plot(p_values, ls=':')
            ax1.tick_params(axis='both',labelsize=20)
    ax1.plot(np.NaN,np.NaN,c=COLORS[s] , label = smooth)
    
    # for sc, scale in enumerate(SCALES):
        
    #     for cu, cutoff in enumerate(CUTOFFS):
            
    #         for t, transform in enumerate(TRANSFORMATIONS):
            
    #             linestyle = color+MARKERSTYLE[t+cu*len(TRANSFORMATIONS)+sc*len(CUTOFFS)+sc*len(TRANSFORMATIONS)]
    #             datapath = os.path.join(RESULTS_DIR, smooth,scale,cutoff,transform, SUB)
    #             p_values = []
    #             accuracies = []
                
    #             if not os.path.exists(datapath):
    #                 continue
    #             for roi in ROIS:
    #                 roi_file = os.path.join(datapath, roi + '.hdf5')
    #                 res = h5py.File(roi_file, 'r')
    #                 accuracies.append(res['accuracy'][()]-1/3)
    #                 p_values.append(res['p_value'][()])
                    
    #             ax0.plot(accuracies,linestyle, label=scale+cutoff+transform,markersize=12)
    #             ax0.plot(accuracies, ls=':')
    #             ax0.tick_params(axis='both',labelsize=20)
                
    #             ax1.plot(p_values, linestyle,label=scale+cutoff+transform,markersize=12)
    #             ax1.plot(p_values, ls=':')
    #             ax1.tick_params(axis='both',labelsize=20)
    # ax1.plot(np.NaN,np.NaN,c=COLORS[s] , label = smooth)

ax1.invert_yaxis()
ax0.set_ylim(0.0)
ax1.set_ylim(0.2)
fig.suptitle('Comparing performance of LDA with different parameters',fontsize=25)
ax0.set_title('Accuracy minus chance',fontsize=25)
#ax0.set_ylabel('Accuracy minus chance in %')
ax1.set_title('p-value',fontsize=25)
ax1.set_xlabel('ROIs',fontsize=25)
                
plt.xticks(np.arange(len(ROIS)),ROIS,rotation=45,fontsize=20)
plt.yticks(fontsize=20)
handles, labels = plt.gca().get_legend_handles_labels()
labels, ids = np.unique(labels, return_index=True)
handles = [handles[i] for i in ids]
fig.legend(handles, labels, loc='upper left',prop={'size': 18})

fig.savefig(os.path.join(RESULTS_DIR,'accuracies.png'))