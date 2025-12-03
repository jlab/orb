from glob import glob
from os.path import join, basename
import yaml
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from compile_data import get_environments, getdata_recovery


def plot_recovery(fp_orb_basedir, settings, num_columns:int=3):
    """Plots contig recovery.
    
    Parameters
    ----------
    fp_orb_basedir : str
        The filepath to Orb's base data dir
    settings : yaml-dict
        General settings for plotting Orb graphs. Here, we need
        a) the order of the environments
        b) pretty labels for environments and assembler
        c) category information about contigs
    num_columns : int
        Maximal number of panels in a row.

    Returns
    -------
    plt figure of the multi-pabel plot.
    """
    # which environments to plot and in which order
    environments = get_environments(fp_orb_basedir, settings)
    # load data
    data_recovery = getdata_recovery(fp_orb_basedir, settings)
    
    fig, axes = plt.subplots(
        int(np.ceil(len(environments) / num_columns)), 
        num_columns * 2, figsize=(2 * num_columns * 4, np.ceil(len(environments) / num_columns) * 5),
        gridspec_kw={"wspace": 0.31, "hspace": 0.3})
    
    for i, environment in enumerate(environments):
        orb = data_recovery[data_recovery['environment'] == environment]

        # bad contigs
        ax_bad = axes[i // num_columns, (i % num_columns) * 2]
        ax_bad.invert_xaxis()
        orb.loc[:, [c['label'] for _, c in settings['contig_classes'].items() if c['class'] == 'bad']].plot(kind='barh', stacked=True, ax=ax_bad, color={c['label']: c['color'] for _, c in settings['contig_classes'].items()})
        ax_bad.set_xlabel("number contigs")
        ax_top_bad = ax_bad.twiny()
        ax_top_bad.xaxis.set_label_position('top')
        ax_top_bad.set_xticks([])
        ax_top_bad.set_xlabel("bad: artifacts")
        ax_bad.xaxis.set_label_coords(1, -0.08)
        ax_bad.text(-0.5, 1.05, chr(97+i), transform=ax_bad.transAxes, fontsize=16, fontweight='bold',)

        # good contigs
        ax_good = axes[i // num_columns, (i % num_columns) * 2  + 1]
        orb.loc[:, [c['label'] for _, c in settings['contig_classes'].items() if c['class'] != 'bad']].plot(kind='barh', stacked=True, ax=ax_good, color={c['label']: c['color'] for _, c in settings['contig_classes'].items()})
        ax_good.set_yticks([])
        ax_good.set_xlabel("good: recovered")
        ax_good.xaxis.set_label_position('top')
    
        # concat right (=good) axis directly adjacent to left (=bad) axis
        ax_good.set_position([
            ax_bad.get_position().x1,
            ax_good.get_position().y0,
            ax_good.get_position().width,
            ax_good.get_position().height
            ])
        ax_bad.set_title(settings['labels']['environments'].get(environment, environment), loc='right', horizontalalignment='center')
    
        # one legend for all panels
        if i+1 == len(environments):
            handles = ax_good.get_legend_handles_labels()[0] + list(reversed(ax_bad.get_legend_handles_labels()[0]))
            labels = ax_good.get_legend_handles_labels()[1] + list(reversed(ax_bad.get_legend_handles_labels()[1]))
            ax_good.legend(handles, labels, ncol=8, bbox_to_anchor=(-0.1, -0.20))
        else:
            ax_good.legend().remove()
        ax_bad.legend().remove()

    return fig


# # above function should be called as following:
# # 1) load plotting settings
# with open("/homes/sjanssen/Git/jlab/orb/plotting/style.yaml", "r") as f:
#    settings = yaml.safe_load(f)
# # 2) then call plot_recovery with the filepath to Orb's base dir and the loaded settings
# _ = plot_recovery('/vol/jlab/tlin/all_project/nf_results/orb', settings)