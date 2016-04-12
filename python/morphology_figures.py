# Generate the figures for Sutton et al., Quantifying STG morphology

"""
Find a way to import 
1. neuron_getProperties.py
2. pretty_plot.py
3. Data dictionaries (first part of each sub-section below)
"""

# All statistics were performed in R. Dictionaries were saved as csv
#   files and analyzed in R. Descriptive statistics were done with 
#   apply or aggregate and comparisons were done with TukeyHSD.

#########################################################################
# Use this to re-assemble all dicts from the geofiles
import os
fdir = '~/data/morphology/morphology-hoc-files/hocs/all/' # Change as needed
fils = os.listdir(fdir)
fils = [f for f in fils if f.split('.')[-1] == 'hoc'
fils = [fdir + f for f in fils]
geofiles = [demoReadsilent(f) for f in fils]

# Paths
paths_lengths, paths_torts = [], [] 
for g in geofiles:
  p, t = path_lengths2(g)
  paths_lengths.append(p)
  paths_torts.append(t)
# Tip stuffs
tip_to_tip_dists = [tip_to_tip(g) for g in geofiles]
tip_dists = [where_tips(g, False) for g in geofiles]
# Branch angles
angles = [branch_angles(g) for g in geofiles]
# Torques
torqs = [get_torques(g) for g in geofiles]
# Asymmetry
path_asym = [path_asymmetry(g) for g in geofiles]
# Sholl
sholl = [hooser_sholl(g)[0] for g in geofiles]
# Neuropil fitting
npil_fit = [neuropil_fit(g) for g in geofiles]
# Fractal dimension
frac_dim = [fractal_dimension(g)[0][0] for g in geofiles]

# Note: Soma position and subtrees require feedback from the user;
# For these, please see the IPython notebook "analysis" at
# github.com/marderlab/QuantifyingMorphology

figdir = '/home/alex/data/morphology/848/848_117b/images/'
#########################################################################
# Figure 3 - path lengths, tips stuff
# Path lengths
length_labs, length_vals = condition_by_name(props['cellTypes'], props['path_lengths'])
plot_hori_hist(length_vals, length_labs, axes=['','Length (um)'], switch=True)
# Tip to tip path distances
tipp_labs, tipp_vals = condition_by_name(props['cellTypes'], props['tip_tip_path'])
plot_hori_hist(tipp_vals, tipp_labs, axes=['','Distance (um)'], switch=True)
# Tip to tip euclidean distances
tipe_labs, tipe_vals = condition_by_name(props['cellTypes'], props['tip_tip_euclid'])
plot_hori_hist(tipe_vals, tipe_labs, axes=['','Distance (um)'], switch=True)
# Tip to center of neuropil distance
tipc_labs, tipc_vals = condition_by_name(props['cellTypes'], props['tip_center_dists'])
plot_hori_hist(tipc_vals, tipc_labs, axes=['','Distance from center (um)'], switch=True)

#########################################################################
# Figure 4 - path tortuosity, etc
tort_labs, tort_vals = condition_by_name(props['cellTypes'], props['path_torts'])
plot_hori_hist(tort_vals, tort_labs, axes=['','Tortuosity'], switch=True, rrange=[1,7])



#########################################################################
# Figure 5 - branch angles
ang_labs, ang_vals = condition_by_name(props['cellTypes'], props['branch_angles'])
fnames = ['GM_branch-angles.png', 'LG_branch-angles.png', 
          'PD_branch-angles.png', 'LP_branch-angles.png']
fnames = [figdir+f for f in fnames]
for u in range(4):
  circular_hist(ang_vals[u*4:(u+1)*4], ang_labs[u*4:(u+1)*4], same=u)
# For an overview
plot_hori_hist(ang_vals, ang_labs, axes=['','Branch angle (degrees)'], switch=True)
# Gouped
mean_scatter(props['branch_angles'], props['cellTypes'], axes=['','Angle (deg)'])

#########################################################################
# Figure 6 - torques
torq_labs, torq_vals = condition_by_name(props['cellTypes'], props['torques'])
for u in range(4):
  circular_hist(torq_vals[u*4:(u+1)*4], torq_labs[u*4:(u+1)*4], same=u)
# Grouped
props['torques'] = [rm_nan(i) for i in props['torques']]
mean_scatter(props['torques'], props['cellTypes'], axes=['','Angle (deg)'])

#########################################################################
# Figure 7 - length symmetry




#########################################################################
# Figure 8 - sholl, soma positions




#########################################################################
# Figure 9 - neuropil fit, fractal dimension


frac_dim = [fractal_dimension(g)[0] for g in geofiles]



#########################################################################
# Figure 10 - Subtrees
# data = subtrees_dict; subtrees_dict.keys():
# dict_keys(['files', 'cellTypes', 'wires', 'sum_paths', 'main_paths', 'total'])

subtrees_dict = json.load(open('~/data/morphology/morpho-paper/datasets/subtrees.txt', 'r'))

# Plot the lengths of each subtree
wires_labs, wires_vals = condition_by_name(subtrees_dict['cellTypes'], subtrees_dict['wires'])
hori_scatter(wires_vals, wires_labs, axes=['','Subtree length (um)'], switch=True, llog=True, counts=True)

# Plot the total wiring 
total_labs, total_vals = condition_by_name(subtrees_dict['cellTypes'], subtrees_dict['total'])
pretty_scatter(total_vals, total_labs, axes=['','Total wiring (um)'])

# Plot the main path statistics
main_labs, main_percent = condition_by_name(subtrees_dict['cellTypes'], [subtrees_dict['sum_paths'][u]/subtrees_dict['total'][u] for u in range(16)])
main_l_labs, main_length = condition_by_name(subtrees_dict['cellTypes'], subtrees_dict['sum_paths'])
if main_l_labs == main_labs:
  pretty_scatter(main_length, main_labs, axes=['','Length (um)','% total wiring'], moreD=main_percent)

# Get summary statistics of a desired cell type
def summarize_type(subtrees_dict, cellType='LP'):
  idx = [i for i in range(len(subtrees_dict['cellTypes'])) if 
         subtrees_dict['cellTypes'][i] == cellType]
  for i in idx:
    print('Found %i %s cells' %(len(idx, cellType)))
    print('No. of subtrees: %i' %len(subtrees_dict['wires'][i]))
    print('Max: %.5f, Min: %.5f' %(max(subtrees_dict['wires'][i]), 
                                   min(subtrees_dict['wires'][i])))
    print('Mean: %.5f, Med: %.5f' %(np.mean(subtrees_dict['wires'][i]), 
                                    np.median(subtrees_dict['wires'][i])))
    print('Path wiring: %.5f, Path wiring percent: %.5f'
          %(subtrees_dict['sum_paths'][i], subtrees_dict['sum_paths'][i]/
                                           subtrees_dict['total'][i]))
    print('Total tree wiring: %.5f' %subtrees_dict['total'][i])
    print('No. of axons: %i ' %len(subtrees_dict['main_paths'][i]))



#########################################################################

# Figure 11 - Radius properties
# data = radius; radius.keys():
# dict_keys(['dist_hand', 'handyes', 'soma2_combined', 'cellTypes', 'tips_single', 
# 'soma1_combined', 'soma2_single', 'initial_expand', 'tips_combined', 'soma1_single'])

radius = json.load(open('~/data/morphology/848_848_117b/radius_full.txt', 'r'))

# Plot initial distance
init_lab, init_vals = condition_by_name(radius['cellTypes'], radius['initial_expand'])
_d, d_hand = condition_by_name(radius['cellTypes'], radius['dist_hand'])
if _d == init_lab:
  pretty_scatter(init_vals, init_lab, axes=['','Primary neurite expansion', 'Distance to "hand"'], moreD=d_hand, showleg='upper center')

# Plot all the other things
# Soma 1 (primary)
soma1_s_labs, soma1_s_vals = condition_by_name(radius['cellTypes'], radius['soma1_single'])            
hori_scatter(soma1_s_vals, soma1_s_labs, switch=True, title='soma1', axes=['','Daughter-Parent Ratio'])
soma1_c_labs, soma1_c_vals = condition_by_name(radius['cellTypes'], radius['soma1_combined'])                          
hori_scatter(soma1_c_vals, soma1_c_labs, switch=True, title='soma1', axes=['','sum(Daughter)-Parent Ratio'])           
# Soma 2 (secondary)
soma2_s_labs, soma2_s_vals = condition_by_name(radius['cellTypes'], radius['soma2_single'])                     
hori_scatter(soma2_s_vals, soma2_s_labs, switch=True, title='soma2', axes=['','Daughter-Parent Ratio'])          
soma2_c_labs, soma2_c_vals = condition_by_name(radius['cellTypes'], radius['soma2_combined'])            
hori_scatter(soma2_c_vals, soma2_c_labs, switch=True, title='soma2', axes=['','sum(Daughter)-Parent Ratio']) 
# Tips
tips_s_labs, tips_s_vals = condition_by_name(radius['cellTypes'], radius['tips_single'])
hori_scatter(tips_s_vals, tips_s_labs, switch=True, axes=['','Daughter-Parent Ratio'])
tips_c_labs, tips_c_vals = condition_by_name(radius['cellTypes'], radius['tips_combined'])
hori_scatter(tips_c_vals, tips_c_labs, switch=True, axes=['','sum(Daughter)-Parent Ratio'])

