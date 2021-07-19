#!/gpfs01/bartels/user/vplikat/anaconda3/bin/python
# -*- coding: utf-8 -*-

import os
import sys
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
DECODER = 'SVM'
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

COLORS = ['r','b','g','c','m','y','tab:orange','tab:pink','tab:gray','indianred','brown']
SUB_MEAN_COLORS = ['c','m','y']
MARKERSTYLE = ['s','D','^', '>','*','<','o','v','+','.','h','H']

SMOOTH_SIZES = ['nativespace']#,'3mm-smoothed-nativespace']# 
SCALES = ['scale_None', 'scale_z']
CUTOFFS = ['cutoff_inf', 'cutoff_2.0']
TRANSFORMATIONS = ['feat_trans_None', 'feat_trans_PCA']
overall_accuracy = []
overall_p_values = []
try:

    fig, (ax0, ax1) = plt.subplots(2,1,sharex=True)
    fig.set_figwidth(30)
    fig.set_figheight(30)
    
    for s,smooth in enumerate(SMOOTH_SIZES):
        color = COLORS[s]
        overall_color = SUB_MEAN_COLORS[s]
        marker_idx=0
        for root, directories, files in os.walk(os.path.join(RESULTS_DIR, smooth)):
            if 'sub-01' in directories:
                marker_idx += 1
            if len(files)>0:
                markerstyle = color+MARKERSTYLE[marker_idx]
                path_as_list = root.split(os.sep)
                index = path_as_list.index(smooth)+1
                try:
                    p_values = []
                    accuracies = []
                    for roi in ROIS:
                        roi_file = os.path.join(root, roi + '.hdf5')
                        res = h5py.File(roi_file, 'r')
                        accuracies.append(res['accuracy'][()]-1/3)
                        p_values.append(res['p_value'][()])
                    overall_accuracy.append(accuracies)
                    overall_p_values.append(p_values)
                    ax0.plot(accuracies,markerstyle,markersize=12, alpha=.3, label='_'.join(path_as_list[index:-1]))#, label=scale+cutoff+transform)
                    #ax0.plot(accuracies, color=color, ls=':')
                    
                    #ax1.plot(p_values, markerstyle,markersize=12, alpha=.1, label='_'.join(path_as_list[index:-1]))#,label=scale+cutoff+transform)
                    #ax1.plot(p_values, color=color, ls=':')
                    
                except Exception as e:
                    exc_type, exc_obj, exc_tb = sys.exc_info()
                    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                    print(exc_type, fname, exc_tb.tb_lineno,e.args)
            if len(overall_accuracy)==24:
                overall_accuracy_array = np.asarray(overall_accuracy)
                overall_p_values_array = np.asarray(overall_p_values)
                sig_ps = overall_p_values_array<0.05
                ax0.plot(overall_accuracy_array.mean(axis=0), overall_color+MARKERSTYLE[marker_idx],markersize=18)
                ax0.plot(overall_accuracy_array.mean(axis=0),'--', color=COLORS[marker_idx], linewidth=12)
                ax1.plot(np.NaN,np.NaN,'--', linewidth=12, color=COLORS[marker_idx], label=smooth+'_'.join(path_as_list[index:-1]))
                # ax1.plot(overall_p_values_array.mean(axis=0), overall_color+MARKERSTYLE[marker_idx],markersize=25,label='mean over subjects')
                # ax1.plot(overall_p_values_array.mean(axis=0),':')
                #ax1.plot(sig_ps.sum(axis=0), color+MARKERSTYLE[marker_idx],markersize=12, label='_'.join(path_as_list[index:-1]))
                ax1.plot(sig_ps.sum(axis=0), MARKERSTYLE[marker_idx],markersize=12, label='_'.join(path_as_list[index:-1]))
                overall_accuracy = []
                overall_p_values = []
        ax1.plot(np.NaN,np.NaN,c=COLORS[s] , label = smooth)
    
    ax0.tick_params(axis='both',labelsize=20)
    ax1.tick_params(axis='both',labelsize=20)
    #ax1.invert_yaxis()
    ax0.set_ylim(0.0)
    #ax1.set_ylim(.2,-0.0005)
    fig.suptitle('Comparing performance of '+ DECODER +' with different parameters',fontsize=25)
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
    
except Exception as e:
    exc_type, exc_obj, exc_tb = sys.exc_info()
    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
    print(exc_type, fname, exc_tb.tb_lineno,e.args)