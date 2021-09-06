#!/gpfs01/bartels/user/vplikat/anaconda3/bin/python

##########
# HEADER #
##########
# Befor you run this script you have to have previously run the script 
# 'univariate_ROI-analysis.py'.
# The data for this script is derived from an fMRI experiment in which
# subjects viewed different videos. Either magic, control or surprise videos. 
# The experiment was devided into 3 blocks. Each block consisted of 4 
# experimental runs. In each run subjects viewed 24 videos (each video is 
# considered a trial).
# The videos in each block were associated with one object (Balls, Cards and 
# Sticks) and there were 3 magic effects (Appear, Change and Vanish). For each
# magic effect and object there are two trick versions (i.e. Appear1, Appear2,
# Change1,...). This resulted in 6 magic videos per object = 18 magic videos 
# and for every magic video there was a corresponding control video showing
# the same movements without the magical effect. Additionally per object there
# were 3 surprise videos showing unusual surprising actions performed with the
# objects (e.g. eating a playing card).
# After the second run in each block the underlying method behind each magic 
# trick was presented.

#                       TIME
#   ---------------------------------->
#   OBJECT1 R1  R2  Revelation  R3  R4  |
#   OBJECT2 R1  R2  Revelation  R3  R4  |   TIME
#   OBJECT3 R1  R2  Revelation  R3  R4  v

# RUNS: 2*Appear1 Magic 2*Appear2 Magic 2*Appear1 Control 2*Appear2 Control
#       2*Vanish1 Magic 2*Vanish2 Magic 2*Vanish1 Control 2*Vanish2 Control
#       2*Change1 Magic 2*Change2 Magic 2*Change1 Control 2*Change2 Control
#       2*Surprise1     2*Surprise2     2*Surprise13
#       = 24 Videos

# The aim of the experiment was to find neural correlates of surprise and in 
# particular of surprising events we consider "impossible". 
# The data used are beta estimate NIfTI images derived from a GLM using SPM12 
# in MATLAB. 

##########################
# PURPOSE OF THIS SCRIPT #
##########################
# The purpose of this script is to read in the repeated measures ANOVA results
# derived from 'univariate_ROI-analysis.py' and create easy to understand plots

################################
# FUNCTIONALITY OF THIS SCRIPT #
################################
# FIRST STEP
# After getting all important path information and names of ROIs (here we
# devide the ROIs into two subsets, the ROIs generated by the Benson toolbox 
# and the rest generated with the Glasser atlas) and the magic effects (Appear
# Change Vanish), the data frame generated by 'univariate_ROI-analysis.py' is
# read in
# SECOND STEP
# For each subset of ROIs a figure with 3 subplots is generated (one for data 
# pre revelation, post revelation and all data together). Then the script 
# iterates over the magic effects and ROIs and saves the data for one effect
# but all ROIs in a list. With this list a boxplot is generated for every ROI
# based on the data of all subjects. The boxplots for the three effects are 
# all plotted in one subplot.
# THIRD STEP
# For each of the two subsets of ROIs a figure is created with one subplot for
# each ROI. In every subplot we plot the mean over all subjects for each 
# effect (Appear Change Vansih) once for pre revelation data and once for post
# revelation data

######################
# COMMAND LINE FLAGS #
######################
# --smooth, const=0, default=0 
# if smoothed data is used and if so what smoothing kernel 
# --analyzed, const='moment',  default='moment'
# The GLM was either calculated using all scans during one trial or depending
# on a specific moment during each video. This flag decides which of the GLMs
# is used

#############
# LIBRARIES #
#############

# interact with the operating system 
import os
import argparse
from pathlib import Path
import git
# data structuration and calculations
import pandas as pd  # to create data frames
import numpy as np   # most important numerical calculations
import scipy.stats as st
# optimize time performance
import time
# plotting
import matplotlib.pyplot as plt

def draw_plot(data, ax, offset,edge_color='blue', fill_color='white'):
    """
    draw_plot: A function that plots data in boxplots with offset, so that 
    multiple boxplots for one variable but different conditions can be plot
    next to each other.
    data: an array like object, with N x M dimensions, resulting in M boxplots
    ad: an axis object in which the data should be plotted
    offset: how much the boxplot should be shifted away from the default 
    position.
    edge_color: color of boxplot
    fill_color: with to fill the boxplot
    
    returns
    bp: a dictionary of artists
    """
    pos = np.arange(data.shape[1])+offset 
    bp  = ax.boxplot(data, 
                     positions=pos, 
                     widths=0.2, 
                     patch_artist=True, 
                     manage_ticks=False)
    for element in ['boxes', 'whiskers', 'fliers', 'medians', 'caps']:
        plt.setp(bp[element], color=edge_color)
    for patch in bp['boxes']:
        patch.set(facecolor=fill_color)
        
    return bp


# get start time
T_START = time.time()

# create a parser object which handles the input variables
parser = argparse.ArgumentParser()

# add all the input arguments
parser.add_argument("--smooth",     nargs='?', const=0,         
                    default=0,          type=int)    
parser.add_argument("--analyzed",   nargs='?', const='moment',  
                    default='moment',   type=str)    
# parse the arguments to a parse-list(???)
ARGS = parser.parse_args()
# assign values 
SMOOTHING_SIZE  = ARGS.smooth
ANALYZED        = ARGS.analyzed

if ANALYZED == 'moment':
    data_analyzed = 'SpecialMoment'
elif ANALYZED == 'video':
    data_analyzed = 'WholeVideo'

# variables for path selection and data access
HOME            = str(Path.home())
PROJ_DIR        = os.path.join(HOME, 'Documents/Master_Thesis/DATA/MRI')
RAWDATA_DIR     = os.path.join(PROJ_DIR, 'rawdata')
DERIVATIVES_DIR = os.path.join(PROJ_DIR, 'derivatives')
DATA_DIR        = os.path.join(DERIVATIVES_DIR, 'univariate-ROI',
                               data_analyzed)
DATA_DF = pd.read_hdf(os.path.join(DATA_DIR,'prepared_dataframe.hdf5'),
                      key='df')
###############
# define ROIs #
###############
# ROIs are devided into the hypothesis driven visual areas provided by the 
# Benson toolbox (which uses the Wang 2015 atlas)
BENSON_ROIS = [
    'V1', 'V2', 'V3', 'hV4',
    'V3A', 'V3B',  
    'LO', 'VO', 
    'IPS','FEF'
    ]

# and additional ROIs extracted from univariate results by Danek et al (2014)
GLASSER_ROIS = [
    'PHT', 'PF', 
    'ACC', 'PCC', 
    'IFG', 'aINSULA', 'IFJ'
    ]

# Combine the two ROI sets and give them names
ROIS = [BENSON_ROIS, GLASSER_ROIS]
ROI_SET_NAMES = ['Benson ROIs', 'Glasser ROIs']

# the magical effects used in the magic study. For each effect there will be 
# one boxplot
EFFECTS = [
    'Appear',
    'Change',
    'Vanish'
    ]

# colors for the boxplots and therefore magic effects (apper: red, change:green ...)
COLORS = ['red','green','blue']

# Markers for pre post comparision
MARKERSTYLE = ['^', 'v']

# in the experiment the same videos were shown to subjects before (pre) and 
# after (post) presenting them the underlying method behind a magic trick.
# These two conditions shall be compared
PRE_POST = ['pre','post']
##############
# BOXPLOTS ! #
##############
# how much each boxplot should be shifted.
OFFSETS = [-0.3,0,0.3]
# an empty list, that will be filled with dictionaries of artists (returned by
# draw_plot), to later create a legend in the figure
bps = []

# iterate over ROI subsets, since we want a figure for each subset
for rs,roi_set in enumerate(ROIS):
    # figure settings
    fig, ax = plt.subplots(3,1,sharex=True)
    fig.set_figwidth(30)
    fig.set_figheight(30)
    
    # iterate over effect (Appear, Change, Vanish)
    for e,eff in enumerate(EFFECTS):
        # and over conditions (pre, post) to get data from only one effect and
        # one condition (e.g. Vanish pre) to draw a boxplot from
        for c,cond in enumerate(PRE_POST):
            # empty list to be filled with data from one condition and effect
            # but all subjects. For each ROI in the subset the data will be
            # appended
            sub_data = []
            
            # iterate over ROIs in subset
            for roi in roi_set:
                # select all data from the current roi (DATA_DF.loc[:,roi]) 
                # but only from the current effect (DATA_DF.Effect==eff) and
                # the current condition (DATA_DF.pre_post==cond). When using 
                # more than one 'filter' parentesis are necessary
                effect_data = DATA_DF.loc[:,roi][(DATA_DF.Effect==eff)&
                                                  (DATA_DF.pre_post==cond)]
                sub_data.append(effect_data)
            
            # turn sub_data list into an array an transverse it for draw_plot
            # otherwise you get one boxplot per subject using data from all 
            # rois
            sub_data = np.array(sub_data)
            draw_plot(data=sub_data.T,
                      ax=ax[c],
                      offset=OFFSETS[e],
                      edge_color=COLORS[e],
                      fill_color='white')
        # Do the same for pre AND post data
        sub_data = []
        # iterate over ROIs in subset
        for roi in roi_set:
            # select all data from the current roi (DATA_DF.loc[:,roi]) 
            # but only from the current effect (DATA_DF.Effect==eff) 
            effect_data = DATA_DF.loc[:,roi][(DATA_DF.Effect==eff)]
            sub_data.append(effect_data)
        
        # turn sub_data list into an array an transverse it for draw_plot
        # otherwise you get one boxplot per subject using data from all 
        # rois
        sub_data = np.array(sub_data)
        bp = draw_plot(data=sub_data.T,
                  ax=ax[2],
                  offset=OFFSETS[e],
                  edge_color=COLORS[e],
                  fill_color='white')
        # add the dictionary of artists to bps list. We need only one per 
        # effect
        bps.append(bp)
    
    ##################
    # FIGURE MAKE UP #
    ##################
    # create a legend in the each subplot
    for a in ax: a.legend([b["boxes"][0] for b in bps], EFFECTS)
    # set the x ticks to the ROIs in subset
    plt.xticks(np.arange(len(roi_set)),roi_set,rotation=45,fontsize=12)
    # set title and fontsizes
    ax[0].set_title('Before revelation', fontsize=25)
    ax[1].set_title('After revelation', fontsize=25)
    ax[2].set_title('All data together', fontsize=25)
    for a in ax: a.tick_params(axis='both',labelsize=20)
    fig.savefig(os.path.join(DATA_DIR,ROI_SET_NAMES[rs]+'_boxplots.png'))

################
# EFFECT MEANS #
################
# figure settings
fig = plt.figure(figsize=(30,30))

# iterate over Benson ROIs 
for r,roi in enumerate(BENSON_ROIS):
    # For each ROI we create a single subplot. 
    # We have selected 10 from the Benson ROIs therefore we go for 4 rows and
    # 3 columns
    num_rows = 4
    num_cols = 3
    # The row is ROI index divided by number of rows the column is ROI index 
    # modulo number of rows
    # The subplot index is column + row*number of rows +1
    row = r//num_rows
    col = r%num_rows
    idx = col+row*num_rows+1
    # create subplot and set title, ticks etc.
    ax = fig.add_subplot(num_rows,num_cols,idx)
    ax.set_title(roi, fontsize=25)
    ax.set_xticks([0,1,2])
    ax.set_xticklabels(EFFECTS)
    ax.tick_params(axis='both',labelsize=20)
    # iterate over conditions (pre, post)
    for c,cond in enumerate(PRE_POST):
        # empty list to be filled with the mean of data from one condition 
        means = []
        stds = []
        # iterate over effect (Appear, Change, Vanish)
        for e,eff in enumerate(EFFECTS):
            # select all data from the current roi (DATA_DF.loc[:,roi]) 
            # but only from the current effect (DATA_DF.Effect==eff) and
            # the current condition (DATA_DF.pre_post==cond). When using 
            # more than one 'filter' parentesis are necessary
            effect_data = DATA_DF.loc[:,roi][(DATA_DF.Effect==eff)&
                                              (DATA_DF.pre_post==cond)]
            means.append(effect_data.mean())
            stds.append(effect_data.std())
        
        # plot the mean value of the effects once with markers and once with
        # a dashed line
        ax.plot(means, color = COLORS[c], marker=MARKERSTYLE[c], label=cond)
        ax.plot(means, color = COLORS[c], linestyle='--')
        
# create a nice legend
legend = plt.legend(title='Legend', 
                    bbox_to_anchor=(1.3,0), 
                    loc='lower right', 
                    fontsize=20)
# fontsize of legend's title is not affected by argument 'fontsize' in function
# plt.legend. Therefore we change it here
plt.setp(legend.get_title(),fontsize=20)
fig.savefig(os.path.join(DATA_DIR,'BensonROIs_means.png'))

# figure settings
fig = plt.figure(figsize=(30,30))
# iterate over Benson ROIs 
for r,roi in enumerate(GLASSER_ROIS):
    # For each ROI we create a single subplot. 
    # We have selected 7 from the Glasser ROIs therefore we go for 3 rows and
    # 3 columns
    num_rows = 3
    num_cols = 3
    # The row is ROI index divided by number of rows the column is ROI index 
    # modulo number of rows
    # The subplot index is column + row*number of rows +1
    row = r//num_rows
    col = r%num_rows
    idx = col+row*num_rows+1
    # create subplot and set title, ticks etc.
    ax = fig.add_subplot(num_rows,num_cols,idx)
    ax.set_title(roi, fontsize=25)
    ax.set_xticks([0,1,2])
    ax.set_xticklabels(EFFECTS)
    ax.tick_params(axis='both',labelsize=20)
    # iterate over conditions (pre, post)
    for c,cond in enumerate(PRE_POST):
        # empty list to be filled with data from one condition 
        means = []
        stds = []
        # iterate over effect (Appear, Change, Vanish)
        for e,eff in enumerate(EFFECTS):
            # select all data from the current roi (DATA_DF.loc[:,roi]) 
            # but only from the current effect (DATA_DF.Effect==eff) and
            # the current condition (DATA_DF.pre_post==cond). When using 
            # more than one 'filter' parentesis are necessary
            effect_data = DATA_DF.loc[:,roi][(DATA_DF.Effect==eff)&
                                              (DATA_DF.pre_post==cond)]
            means.append(effect_data.mean())
            stds.append(effect_data.std())
        # plot the mean value of the effects once with markers and once with
        # a dashed line
        ax.plot(means, color = COLORS[c], marker=MARKERSTYLE[c], label=cond)
        ax.plot(means, color = COLORS[c], linestyle='--')
# create a nice legend  
legend = plt.legend(title='Legend', 
                    bbox_to_anchor=(0.3,-.4), 
                    loc='lower right', 
                    fontsize=20)
# fontsize of legend's title is not affected by argument 'fontsize' in function
# plt.legend. Therefore we change it here
plt.setp(legend.get_title(),fontsize=20)
fig.savefig(os.path.join(DATA_DIR,'GlasserROIs_means.png'))