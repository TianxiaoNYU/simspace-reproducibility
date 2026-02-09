import warnings
warnings.filterwarnings("ignore")

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns
import colorcet as cc
import os

import simspace as ss

########## Panel A:
shape = (100, 100)
script_dir = os.path.dirname(os.path.abspath(__file__))
param_file = os.path.join(script_dir, 'fitted_params.json')
param_1 = ss.util.load_params(param_file)

n_group = param_1['n_group']
diag_length = len(param_1['theta_list'][0])
n_state = int((1 + np.sqrt(1 + 8 * diag_length)) / 2)
if n_state * (n_state - 1) != 2 * diag_length:
    raise ValueError("Invalid theta matrix size.")

niche_theta = np.zeros((n_group, n_group))
niche_theta[np.triu_indices(n_group, 1)] = param_1['niche_theta']
niche_theta = niche_theta + niche_theta.T - np.diag(niche_theta.diagonal())
np.fill_diagonal(niche_theta, 1)

theta_list = []
for i in range(n_group):
    theta_tmp = np.zeros((n_state, n_state))
    theta_tmp[np.triu_indices(n_state, 1)] = param_1['theta_list'][i]
    theta_tmp = theta_tmp + theta_tmp.T - np.diag(theta_tmp.diagonal())
    np.fill_diagonal(theta_tmp, 1)
    theta_list.append(theta_tmp)

density_replicates = np.array(param_1['density_replicates'])
density_replicates[density_replicates < 0] = 0
phi_replicates = param_1['phi_replicates']

sim1 = ss.SimSpace(
    shape=shape,
    num_states=n_state,
    num_iterations=4,
    theta=theta_list,
    phi=phi_replicates,
    neighborhood=ss.spatial.generate_offsets(3, 'manhattan'), 
    random_seed=0,
    )
sim1.initialize()
sim1.create_niche(num_niches=n_group, n_iter=6, theta_niche=niche_theta)
sim1.gibbs_sampler()    # Gibbs sampling
sim1.plot_niche(figsize=(5, 5), dpi=100)

# Save the figure instead of showing it
plt.savefig(f'{script_dir}/Fig1_panel_A.png', bbox_inches='tight', dpi=300)
plt.close()  # Close the figure to prevent it from showing


########## Panel B:
n_group = param_1['n_group']
diag_length = len(param_1['theta_list'][0])
n_state = int((1 + np.sqrt(1 + 8 * diag_length)) / 2)
if n_state * (n_state - 1) != 2 * diag_length:
    raise ValueError("Invalid theta matrix size.")

niche_theta = np.zeros((n_group, n_group))
niche_theta[np.triu_indices(n_group, 1)] = param_1['niche_theta']
niche_theta = niche_theta + niche_theta.T - np.diag(niche_theta.diagonal())
np.fill_diagonal(niche_theta, 1)

theta_list = []
for i in range(n_group):
    theta_tmp = np.zeros((n_state, n_state))
    theta_tmp[np.triu_indices(n_state, 1)] = param_1['theta_list'][i]
    theta_tmp = theta_tmp + theta_tmp.T - np.diag(theta_tmp.diagonal())
    np.fill_diagonal(theta_tmp, 1)
    theta_list.append(theta_tmp)

density_replicates = np.array(param_1['density_replicates'])
density_replicates[density_replicates < 0] = 0
phi_replicates = param_1['phi_replicates']

sim1 = ss.SimSpace(
    shape = shape,
    num_states = n_state,
    num_iterations= 4,
    theta=theta_list,
    phi=phi_replicates,
    neighborhood=ss.spatial.generate_offsets(3, 'manhattan'), 
    random_seed=42,
    )
sim1.initialize()
sim1.create_niche(num_niches=n_group, n_iter=6, theta_niche=niche_theta)
sim1.plot_niche(figsize=(5, 5))
plt.savefig(f'{script_dir}/Fig1_panel_B1.png', bbox_inches='tight', dpi=300)
plt.close()  # Close the figure to prevent it from showing

sim1.gibbs_sampler()    # Gibbs sampling
sim1.plot_grid(figsize=(5, 5), dpi=150)
# Save the figure instead of showing it
plt.savefig(f'{script_dir}/Fig1_panel_B2.png', bbox_inches='tight', dpi=300)
plt.close()  # Close the figure to prevent it from showing
sim1.density_sampler(density_replicates)  # Cell density of each niche
sim1.plot_grid(figsize=(5, 5), dpi=150)
plt.savefig(f'{script_dir}/Fig1_panel_B3.png', bbox_inches='tight', dpi=300)
plt.close()  # Close the figure to prevent it from showing
sim1.perturbation(step = 0.2)

cmap = sns.color_palette(cc.glasbey, sim1.meta['state'].nunique())

fig = plt.figure(figsize=(5,5), dpi=300)
ax = fig.add_subplot()
ax.set_aspect('equal')
ax = sns.scatterplot(data=sim1.meta, x='col', y='row', hue='state', s=20, palette=cmap)
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title='state')
ax.set_xlabel('X')
ax.set_ylabel('Y')
plt.tight_layout()
plt.savefig(f'{script_dir}/Fig1_panel_B4.png', bbox_inches='tight')
plt.close()  # Close the figure to prevent it from showing
