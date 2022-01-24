#!/gpfs01/bartels/user/vplikat/anaconda3/bin/python

##########
# HEADER #
##########
# This script is the first one to run before doing any ROI based univariate
# analysis!
# The data for this script is derived from an fMRI experiment in which subjects
# viewed different videos. Either magic, control or surprise videos. 
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
# The purpose of this script is to read in the data frame containing beta 
# values of a set of ROIs and surprise ratings (created by 'create_Rating_Betas_DF.py').
# and perform a set of analyses on these data:
# Test if there are main effects of the different experimental conditions
# --> Pre-Post, Magic-Control-Surprise, Object (or Run), Effect
# Test for interaction effects (how to interprete)

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
from statsmodels.stats.anova import AnovaRM
import pingouin as pg
import matplotlib.pyplot as plt
# optimize time performance
import time

# get start time
T_START = time.time()

################################
# Command line input arguments #
################################

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
else:
    raise

# variables for path selection and data access
HOME            = str(Path.home())
PROJ_DIR        = os.path.join(HOME, 'Documents/Magic_fMRI/DATA/MRI')
RAWDATA_DIR     = os.path.join(PROJ_DIR, 'rawdata')
DERIVATIVES_DIR = os.path.join(PROJ_DIR, 'derivatives')
DATA_DIR        = os.path.join(DERIVATIVES_DIR, 'univariate-ROI',
                               data_analyzed,'EveryVideo')

data_frame = pd.read_csv(filepath_or_buffer=os.path.join(DATA_DIR,'ratings_df.csv'))

VIDEO_TYPE = [
    'Magic',
    'Control',
    'Surprise'
    ]

EFFECTS = [
    'Appear',
    'Change',
    'Vanish'
    ]

############
# ANALYSES #
############
# correlation between behavioral surprise ratings and beta images

# first remove missing values and average over conditions
# since pingouin 0.5.0 -> (not remove_rm_na anymore)
na_removed = data_frame.pivot_table(index='ids', 
                                    columns=['PrePost', 'Type'],
                                    values='Ratings',
                                    aggfunc='mean',
                                    dropna=True,
                                    fill_value=None)
avg_data = na_removed.melt(ignore_index=False, value_name='Ratings').reset_index()

# repeated measures ANOVA of pre-post x Video type

prepost_type_spher = pg.sphericity(data=avg_data,
                      dv='Ratings', 
                      subject='ids',
                      within=['PrePost', 'Type'])
pre_post_norm = pg.normality(data=avg_data,
                      dv='Ratings',
                      group= 'PrePost')

type_norm = pg.normality(data=avg_data,
                      dv='Ratings',
                      group= 'Type')

if all(pre_post_norm.normal.values) and all(type_norm.normal.values):
    
    # print if sphericity is given for the data
    print('Spericity for rating data is given: {}\n'.format(prepost_type_spher.spher))
    aov_res = pg.rm_anova(data=avg_data,
                          dv='Ratings', 
                          subject='ids',
                          within=['PrePost', 'Type'],
                          correction = not prepost_type_spher.spher,
                          detailed=True)
    
    aov_res.to_csv(path_or_buf=os.path.join(DATA_DIR,'aov_PrePost_Type_results.csv'),
                   index=False, float_format='%.3f')
    
    # report results of anov
    if prepost_type_spher.spher:
        p_to_report = aov_res['p-unc']
    else:
        p_to_report = aov_res['p-GG-corr']
    print('{} main effect: F = {} (p={})\n'.format(aov_res.Source[0], aov_res.F[0],
                                           p_to_report[0]))
    print('{} main effect: F = {} (p={})\n'.format(aov_res.Source[1], aov_res.F[1],
                                           p_to_report[1]))
    print('{} interaction effect: F = {} (p={})\n'.format(aov_res.Source[2], aov_res.F[2],
                                           p_to_report[2]))
    
    fig = plt.figure()
    for vid_type in VIDEO_TYPE:
        rating_pre = avg_data.loc[:,'Ratings'][(avg_data.Type==vid_type)&
                                                 (avg_data.PrePost==0)]
        rating_post = avg_data.loc[:,'Ratings'][(avg_data.Type==vid_type)&
                                                  (avg_data.PrePost==1)]
        pre_post_means  = [rating_pre.mean(),rating_post.mean()]
        pre_post_means  = np.array(pre_post_means)
        pre_post_sems   = [rating_pre.sem(),rating_post.sem()]
        pre_post_sems   = np.array(pre_post_sems)
        upper = pre_post_means+pre_post_sems
        lower = pre_post_means-pre_post_sems
        plt.plot(pre_post_means)
        plt.fill_between([0,1], upper, lower, alpha=0.2)
    
    # # if the interaction effect is significant
    if p_to_report[2]<0.05:
        # pairwise t-tests
        
        post_hoc_res_PrePost = pg.pairwise_ttests(data=avg_data, 
                                                  dv='Ratings', 
                                                  within= ['PrePost','Type'],
                                                  subject='ids',
                                                  alternative='two-sided',
                                                  effsize='cohen',
                                                  padjust='bonf')
        
        post_hoc_res_Type = pg.pairwise_ttests(data=avg_data, 
                                               dv='Ratings', 
                                               within= ['Type','PrePost'],
                                               subject='ids',
                                               alternative='two-sided',
                                               effsize='cohen',
                                               padjust='bonf')
        
        post_hoc_res = pd.concat ([post_hoc_res_PrePost, post_hoc_res_Type])
        post_hoc_res.to_csv(path_or_buf=os.path.join(DATA_DIR,'post_hoc_PrePost_Type_results.csv'),
                            index=False, float_format='%.3f')
        
else:
    # do non-parametric statistics
    # for now just do a friedman test for the video types, pooled and one for 
    # pre and post each
    # and a wilkoxon test for pre post and for every video type pre and post
    # wilcoxon test for pre post pooled
    pre_data = avg_data[avg_data.PrePost==0].Ratings
    post_data = avg_data[avg_data.PrePost==1].Ratings
    
    pre_post_pooled_wilcox = pg.wilcoxon(x=pre_data,
                                         y=post_data,
                                         alternative='greater')
    pre_post_pooled_wilcox.to_csv(os.path.join(DATA_DIR,'pre_post_wilcoxon.csv'),
                                  index=False, float_format='%.3f')
    
    x_var   = []
    y_var   = []
    data    = []
    
    # wilcoxon test for pre post Magic
    pre_data_magic = avg_data[(avg_data.PrePost==0) & 
                                (avg_data.Type=='Magic')].Ratings
    post_data_magic = avg_data[(avg_data.PrePost==1) & 
                                 (avg_data.Type=='Magic')].Ratings
    
    pre_post_magic_wilcox = pg.wilcoxon(x=pre_data_magic,
                                        y=post_data_magic,
                                        alternative='greater')
    x_var.append('pre')
    y_var.append('post')
    data.append('Magic')
    
    # wilcoxon test for pre post Surprise
    pre_data_surprise = avg_data[(avg_data.PrePost==0) & 
                                   (avg_data.Type=='Surprise')].Ratings
    post_data_surprise = avg_data[(avg_data.PrePost==1) & 
                                    (avg_data.Type=='Surprise')].Ratings
    
    pre_post_surprise_wilcox = pg.wilcoxon(x=pre_data_surprise,
                                           y=post_data_surprise,
                                           alternative='greater')
    x_var.append('pre')
    y_var.append('post')
    data.append('Surprise')
    
    # wilcoxon test for pre post Control
    pre_data_control = avg_data[(avg_data.PrePost==0) & 
                                  (avg_data.Type=='Control')].Ratings
    post_data_control = avg_data[(avg_data.PrePost==1) & 
                                   (avg_data.Type=='Control')].Ratings
    
    pre_post_control_wilcox = pg.wilcoxon(x=pre_data_control,
                                          y=post_data_control,
                                          alternative='greater')
    x_var.append('pre')
    y_var.append('post')
    data.append('Control')
    
    # friedman for video types all data:
    friedman_type_pooled = pg.friedman(data=avg_data,
                                       dv='Ratings', 
                                       subject='ids',
                                       within='Type')

    if friedman_type_pooled['p-unc'].values[0]<0.05:
        # do post hoc wilcoxon tests
        pooled_data_magic = avg_data[avg_data.Type=='Magic'].Ratings
        pooled_data_surprise = avg_data[avg_data.Type=='Surprise'].Ratings
        pooled_data_control = avg_data[avg_data.Type=='Control'].Ratings
        
        magic_surprise_pooled_wilcox = pg.wilcoxon(x=pooled_data_magic,
                                                   y=pooled_data_surprise,
                                                   alternative='greater')
        x_var.append('Magic')
        y_var.append('Surprise')
        data.append('pooled')
        magic_control_pooled_wilcox = pg.wilcoxon(x=pooled_data_magic,
                                                  y=pooled_data_control,
                                                  alternative='greater')
        x_var.append('Magic')
        y_var.append('Control')
        data.append('pooled')
        surprise_control_pooled_wilcox = pg.wilcoxon(x=pooled_data_surprise,
                                                     y=pooled_data_control,
                                                     alternative='greater')
        x_var.append('Surprise')
        y_var.append('Control')
        data.append('pooled')
    
    # friedman for video types pre revelation:
    pre_data = avg_data[avg_data.PrePost==0]
    friedman_type_pre = pg.friedman(data=pre_data,
                                       dv='Ratings', 
                                       subject='ids',
                                       within='Type')
    
    if friedman_type_pre['p-unc'].values[0]<0.05:
        # do post hoc wilcoxon tests
        magic_surprise_pre_wilcox = pg.wilcoxon(x=pre_data_magic,
                                                y=pre_data_surprise,
                                                alternative='greater')
        x_var.append('Magic')
        y_var.append('Surprise')
        data.append('pre')
        magic_control_pre_wilcox = pg.wilcoxon(x=pre_data_magic,
                                               y=pre_data_control,
                                               alternative='greater')
        x_var.append('Magic')
        y_var.append('Control')
        data.append('pre')
        surprise_control_pre_wilcox = pg.wilcoxon(x=pre_data_surprise,
                                                  y=pre_data_control,
                                                  alternative='greater')
        x_var.append('Surprise')
        y_var.append('Control')
        data.append('pre')
    
    # friedman for video types pre revelation:
    post_data = avg_data[avg_data.PrePost==1]
    friedman_type_post = pg.friedman(data=post_data,
                                          dv='Ratings', 
                                          subject='ids',
                                          within='Type')

    if friedman_type_post['p-unc'].values[0]<0.05:
        # do post hoc wilcoxon tests
        magic_surprise_post_wilcox = pg.wilcoxon(x=post_data_magic,
                                                 y=post_data_surprise,
                                                 alternative='greater')
        x_var.append('Magic')
        y_var.append('Surprise')
        data.append('post')
        magic_control_post_wilcox = pg.wilcoxon(x=post_data_magic,
                                                y=post_data_control,
                                                alternative='greater')
        x_var.append('Magic')
        y_var.append('Control')
        data.append('post')
        surprise_control_post_wilcox = pg.wilcoxon(x=post_data_surprise,
                                                   y=post_data_control,
                                                   alternative='greater')
        x_var.append('Surprise')
        y_var.append('Control')
        data.append('post')
    
    
    post_hoc_df = pd.concat([pre_post_magic_wilcox,         # pre post
                             pre_post_surprise_wilcox, 
                             pre_post_control_wilcox,       
                             magic_surprise_pooled_wilcox,  # type pooled
                             magic_control_pooled_wilcox, 
                             surprise_control_pooled_wilcox,
                             magic_surprise_pre_wilcox,     # type pre
                             magic_control_pre_wilcox, 
                             surprise_control_pre_wilcox, 
                             magic_surprise_post_wilcox,    # type post
                             magic_control_post_wilcox, 
                             surprise_control_post_wilcox])
    
    post_hoc_df['X']=x_var
    post_hoc_df['Y']=y_var
    post_hoc_df['data']=data
    post_hoc_df['p-corrected']=post_hoc_df['p-val']*len(post_hoc_df)
    
    post_hoc_df.to_csv(os.path.join(DATA_DIR,'prepost_type_posthoc_wilcoxon.csv'),
                       index=False, float_format='%.3f')
    
    friedman_type_tests = pd.concat([friedman_type_pooled, 
                                     friedman_type_pre, 
                                     friedman_type_post])
    friedman_type_tests['p-corrected'] = friedman_type_tests['p-unc']*len(friedman_type_tests)
    friedman_type_tests['data'] = ['pooled', 'pre', 'post']
    
    friedman_type_tests.to_csv(os.path.join(DATA_DIR,'friedman_type.csv'),
                       index=False, float_format='%.3f')
    
magic_df = data_frame[data_frame.Type=='Magic']
na_removed = magic_df.pivot_table(index='ids', 
                                    columns=['Objects', 'Effect'],
                                    values='Ratings',
                                    aggfunc='mean',
                                    dropna=True,
                                    fill_value=None)
avg_data = na_removed.melt(ignore_index=False, value_name='Ratings').reset_index()

effect_norm = pg.normality(data=avg_data,
                           dv='Ratings',
                           group= 'Effect')

object_norm = pg.normality(data=avg_data,
                           dv='Ratings',
                           group= 'Objects')

if all(effect_norm.normal.values) and all(object_norm.normal.values):
    # repeated measures ANOVA for effects x objects (only magic videos)

    rm_aov = AnovaRM(data=avg_data, 
                     depvar='Ratings', 
                     subject='ids',
                     within=['Objects', 'Effect'])
    res = rm_aov.fit()
    res.anova_table.to_csv(os.path.join(DATA_DIR,'effect_object_rmANOVA.csv'),
                       index=False, float_format='%.3f')
    # report results of anov
    p_to_report = res.anova_table['Pr > F']
    
    # if spher.spher:
    #     p_to_report = res.anova_table['Pr > F']
    # else:
    #     p_to_report = aov_res['p-GG-corr']
    print('{} main effect: F = {} (p={})\n'.format(res.anova_table.index[0], 
                                                   res.anova_table['F Value'][0],
                                                   p_to_report[0]))
    print('{} main effect: F = {} (p={})\n'.format(res.anova_table.index[1], 
                                                   res.anova_table['F Value'][1],
                                                   p_to_report[1]))
    print('{} interaction effect: F = {} (p={})\n'.format(res.anova_table.index[2], 
                                                          res.anova_table['F Value'][2],
                                                          p_to_report[2]))
    
    ####################################################################
    # check for significant effects and do post hoc tests if there are #
    ####################################################################
    # first check for an interaction. If there is one, do all post hoc
    if res.anova_table['Pr > F'][2]<0.05:
        post_hoc_res_Obj = pg.pairwise_ttests(data=magic_df,
                                              dv='Ratings',
                                              within=['Objects', 'Effect'],
                                              subject='ids',
                                              effsize='cohen',
                                              padjust='bonf')
        post_hoc_res_Eff = pg.pairwise_ttests(data=magic_df,
                                              dv='Ratings',
                                              within=['Effect', 'Objects'],
                                              subject='ids',
                                              effsize='cohen',
                                              padjust='bonf')
        post_hoc_res = pd.concat ([post_hoc_res_Obj, post_hoc_res_Eff])
        post_hoc_res.to_csv(path_or_buf=os.path.join(DATA_DIR,'post_hoc_Obj_Eff_results.csv'),
                            index=False, float_format='%.3f')
        
    elif res.anova_table['Pr > F'][1]<0.05:
        post_hoc_res_Eff = pg.pairwise_ttests(data=magic_df,
                                              dv='Ratings',
                                              within=['Effect'],
                                              subject='ids',
                                              effsize='cohen',
                                              padjust='bonf')
        post_hoc_res_Eff.to_csv(path_or_buf=os.path.join(DATA_DIR,'post_hoc_Eff_results.csv'),
                            index=False, float_format='%.3f')
    elif res.anova_table['Pr > F'][0]<0.05:
        post_hoc_res_Obj = pg.pairwise_ttests(data=magic_df,
                                              dv='Ratings',
                                              within=['Objects'],
                                              subject='ids',
                                              effsize='cohen',
                                              padjust='bonf')
        post_hoc_res_Obj.to_csv(path_or_buf=os.path.join(DATA_DIR,'post_hoc_Obj_results.csv'),
                            index=False, float_format='%.3f')
        
else:
    # do non-parametric statistics
    pass

fig = plt.figure()
eff_values  = []
eff_sems    = []
for eff in EFFECTS:
    eff_ratings = avg_data.loc[:,'Ratings'][(avg_data.Effect==eff)]
    eff_values.append(eff_ratings.mean())
    eff_sems.append(eff_ratings.sem())
    
eff_values  = np.array(eff_values)
eff_sems    = np.array(eff_sems)
upper = eff_values+eff_sems
lower = eff_values-eff_sems
plt.plot(eff_values)
plt.fill_between([0,1,2], upper, lower, alpha=0.2)

# Since appear magic tricks are significantly less surprising than change and
# vanish we compare appear against surprise
# first remove missing values and average over conditions
effect_df = data_frame[data_frame.Type=='Magic']
eff_na_removed = effect_df.pivot_table(index='ids', 
                                    columns='Effect',
                                    values='Ratings',
                                    aggfunc='mean',
                                    dropna=True,
                                    fill_value=None)
eff_avg_data = eff_na_removed.melt(ignore_index=False, value_name='Ratings').reset_index()

type_na_removed = data_frame.pivot_table(index='ids', 
                                    columns='Type',
                                    values='Ratings',
                                    aggfunc='mean',
                                    dropna=True,
                                    fill_value=None)
type_avg_data = type_na_removed.melt(ignore_index=False, value_name='Ratings').reset_index()

xvar = eff_avg_data.Ratings[eff_avg_data.Effect=='Appear']
yvar = type_avg_data.Ratings[type_avg_data.Type=='Surprise']

ttest_res_all = pg.ttest(x=xvar, y=yvar,paired=True,alternative='greater')

eff_na_removed = effect_df.pivot_table(index='ids', 
                                    columns=['Effect','PrePost'],
                                    values='Ratings',
                                    aggfunc='mean',
                                    dropna=True,
                                    fill_value=None)
eff_avg_data = eff_na_removed.melt(ignore_index=False, value_name='Ratings').reset_index()

type_na_removed = data_frame.pivot_table(index='ids', 
                                    columns=['Type','PrePost'],
                                    values='Ratings',
                                    aggfunc='mean',
                                    dropna=True,
                                    fill_value=None)
type_avg_data = type_na_removed.melt(ignore_index=False, value_name='Ratings').reset_index()

xvar = eff_avg_data.Ratings[(eff_avg_data.Effect=='Appear') & 
                                 (eff_avg_data.PrePost==0)]
yvar = type_avg_data.Ratings[(type_avg_data.Type=='Surprise') &
                               (type_avg_data.PrePost==0)]
ttest_res_pre = pg.ttest(x=xvar, y=yvar,paired=True,alternative='greater')

xvar = eff_avg_data.Ratings[(eff_avg_data.Effect=='Appear') & 
                                 (eff_avg_data.PrePost==1)]
yvar = type_avg_data.Ratings[(type_avg_data.Type=='Surprise') &
                               (type_avg_data.PrePost==1)]
ttest_res_post = pg.ttest(x=xvar, y=yvar,paired=True,alternative='greater')

appear_ttest_res = pd.concat ([ttest_res_all, ttest_res_pre, ttest_res_post])
appear_ttest_res['PrePost'] = ['all', 'pre', 'post']

appear_ttest_res.to_csv(path_or_buf=os.path.join(DATA_DIR,
                                                 'appear_surprise_ttest_res.csv'),
                        index=False, float_format='%.3f')

##################
# WRITE LOG FILE #
##################
# We want to save all important information of the script execution
# To get the git hash we have to check if the script was run locally or on the
# cluster. If it is run on the cluster we want to get the $PBS_O_WORKDIR 
# variable, which preserves the location from which the job was started. 
# If it is run locally we want to get the current working directory.

try:
    script_file_directory = os.environ["PBS_O_WORKDIR"]
except KeyError:
    script_file_directory = os.getcwd()
    
try:
    rep = git.Repo(script_file_directory, search_parent_directories=True)
    git_hash = rep.head.object.hexsha
except git.InvalidGitRepositoryError:
    git_hash = 'not-found'

# create a log file, that saves some information about the run script
with open(os.path.join(DATA_DIR,'rating_analysis-logfile.txt'), 'w+') as writer:
    writer.write('Codeversion: {} \n'.format(git_hash))
    writer.write('Time for computation: {}h'.format(str((time.time() - T_START)/3600)))
