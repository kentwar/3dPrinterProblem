import matplotlib.pyplot as plt
import matplotlib
from typing import List
import numpy as np

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
# call plot_compare with appropriate params