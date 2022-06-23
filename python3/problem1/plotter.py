import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
from typing import List
import numpy as np
import os

import constants

sns.set_theme()

def plot_compare(arrs: List[np.ndarray],
                 colors: List[str],
                 linestyles: List[str],
                 labels: List[str],
                 title: str,
                 ylabel: str,
                 filename: str,
                 add_area: bool = True,
                 legend_loc: str = 'lower right'):
    w, h = matplotlib.figure.figaspect(9/16)
    fig = plt.figure(figsize=(w,h))
    
    plt.grid()
    for arr, color, label, linestyle in zip(arrs, colors, labels, linestyles):
        values_mean = np.nanmean(arr, axis=0)
        
        plt.plot(range(len(values_mean)), values_mean, label=label, c=color, lw=2, linestyle=linestyle)
        if add_area:
            values_std = np.std(arr, axis=0)
            
            plt.fill_between(range(len(values_std)), (values_mean - values_std), (values_mean + values_std), color=color, alpha=0.1)
    
    plt.legend(loc=legend_loc)
    plt.title(title)
    plt.ylabel(ylabel)
    plt.xlabel('Generations')
    plt.autoscale(enable=True, axis='x', tight=True)
    plt.savefig(f'./{filename}.png', transparent=True, bbox_inches='tight')
    plt.show()

# TODO
# load npz files and store them in a list

arrs = []
algs = ['simple_ea', 'shd_1_plus_lambda', 'shd_mu_plus_lambda']
colors = ['darkgreen', 'red', 'blue']
linestyles = ['solid']
title = 'Algorithms comparison'
ylabel = 'Fitness'

for alg in algs:
    arrs.append(np.load(os.path.join(constants.results_folder, f'{alg}_res.npz'))['arr_0'])

# call plot_compare with appropriate params
plot_compare(arrs=arrs, colors=colors, linestyles=linestyles, labels=algs, title=title, ylabel=ylabel, filename=os.path.join(constants.results_folder, 'comparison_of_methods'), add_area=True)