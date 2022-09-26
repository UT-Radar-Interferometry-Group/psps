from psps import sario, scr, similarity
import glob
import os
import numpy as np
from matplotlib import pyplot as plt

######################################
# Set-up
######################################
plt.rcParams['image.interpolation'] = "none"
nrow,ncol = 600,600 # number of rows and columns in each interferogram
amp_disp_th = 0.4   # threshold of amplitude dispersion (used in the MLE method)
scr_th = 2          # threshold of Signal-to-Clutter-Ratio (used for
                    # preliminary PS selection before phase similarity
                    # calculation)
med_sim_th = 0.3    # threshold of median phase similarity to remove false
                    # positive PS detections from the initial PS set
alpha = 0.01        # maximum acceptable false positive rate (to determine
                    # the treshold used in maximum phase similarity
                    # calculation)

# list of SAR amplitude files
amplist = glob.glob('data/amplitude/*.amp')
# list of interferograms
ifglist = glob.glob('data/igrams/*.int')
# where to save the scr/similarity files
output_dir = 'ps_data'
if not os.path.exists(output_dir):
    os.mkdir(output_dir)

######################################
# Step 1, select initial PS candidates
######################################

# Calculate amplitude dispersion from amplitude files
print('Calculate amplitude dispersion from amplitude files')
amp_disp = scr.amp_dispersion(amplist,nrow,ncol)
amp_disp.tofile(os.path.join(output_dir,'amp_disp'))
print(f'Amplitude dispersion values have been saved to {output_dir}/amp_disp')

# calculate SCR
print('Calculate amplitude dispersion from amplitude files')
scr_val = scr.cal_scr(ifglist,
                      nrow,
                      ncol,
                      ntilerow=2,
                      ntilecol=2,
                      windowsize = 11,
                      candidate_mask=(amp_disp<amp_disp_th),
                      model='constant',
                      output_dir=output_dir)
scr_val.tofile(os.path.join(output_dir,'scr'))
#scr_val = np.fromfile('ps_data/scr',dtype=np.float32).reshape((nrow,ncol))
print(f'scr values have been saved to {output_dir}/scr')
ps0 = scr_val > scr_th

#################################################
# Step 2, refine PS set based on phase similarity
#################################################
print('Calculate median phase similarity')
med_sim = similarity.median_similarity(
          ifglist,nrow,ncol,ps0=ps0,N=20,rdmin=3,rdmax=50)
ps_refined = med_sim > med_sim_th

####################################################################
# Step 3, recover false negative PS pixels based on phase similarity
####################################################################
# This is an example of how to obtain a set of calibration non-PS pixels 
# to estimate the phase similarity threshold in maximum phase similarity
# calculation. The users do not need to follow exactly the same procedure.
# The calibration pixel set can also be derived from other data sets (e.g.,
# a land-cover map of the area of interest).

# file to read the average InSAR phase correlation
avg_corr_file = 'data/correlation/avg_correlation'
# read in average InSAR phase correlation
avg_corr = np.fromfile(avg_corr_file,dtype=np.float32).reshape((nrow,ncol))
# calculate the 1st percentile of 'avg_corr'
avg_corr_th = np.percentile(avg_corr,1)
# calculate a mask of calibration (non-PS) pixels 
nonps = avg_corr < avg_corr_th

print('Estimate the threshold used for maximum phase similarity calculation')
nonps_sim = similarity.nonps_similarity(ifglist,nrow,ncol,nonps)
sim_th = np.percentile(nonps_sim[nonps],(1-alpha)*100)
print(f'The estimated threshold is {sim_th}')

print('Calculate maximum phase similarity')
max_sim = similarity.max_similarity(
          ifglist,nrow,ncol,ps0=ps_refined,threshold=sim_th,N=20,rdmin=3,rdmax=50)
max_sim.tofile(os.path.join(output_dir,'max_sim'))
#max_sim = np.fromfile(os.path.join(output_dir,'max_sim'),dtype=np.float32).reshape(nrow,ncol)
print('Maximum phase similarity values have been saved to '+os.path.join(output_dir,'max_sim'))

