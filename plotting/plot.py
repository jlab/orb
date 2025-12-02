from glob import glob
from os.path import join, basename
import yaml
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def get_environments(fp_orb_basedir:str, settings) -> [str]:
    """Iterated the base dir of Orb to find environments for which Orb has calculated metrics.
    
    Parameters
    ----------
    fp_orb_basedir : str
        The filepath to Orb's base data dir
    settings : yaml-dict
        General settings for plotting Orb graphs. Here, we need
        a) the list of environments to skip
        b) a potential pre-defined order of the found environments
    
    Returns
    -------
    A list of environment names, i.e. sub-directories in Orb's base dir.

    """
    # collect environments from orb directory
    environments = [basename(fp)
     for fp in glob(join(fp_orb_basedir, '*'))
     if basename(fp) not in settings['skip_environments']]
    # sort environments according to YAML occurrence
    environments = sorted(environments, key=lambda x: list(settings['labels']['environments'].keys()).index(x) if x in settings['labels']['environments'] else float('inf'))

    return environments

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
    environments = get_environments(fp_orb_basedir)
    
    fig, axes = plt.subplots(
        int(np.ceil(len(environments) / num_columns)), 
        num_columns * 2, figsize=(2 * num_columns * 4, np.ceil(len(environments) / num_columns) * 5),
        gridspec_kw={"wspace": 0.31, "hspace": 0.3})
    
    for i, environment in enumerate(environments):
        # load orb data
        orb = pd.read_csv(join(fp_orb_basedir, environment, 'mergedresults', '%s_all_scores.tsv' % environment), sep="\t", index_col=0)
        # rename metrics to names used in publication
        orb = orb.rename(index={k: c['label'] for k, c in settings['contig_classes'].items()})
        # add information about missed blocks
        trueBlocks = pd.read_csv(join(fp_orb_basedir, environment, 'mergedataframessummaries', '%s_all_contigs.tsv' % environment), sep="\t", index_col=0).loc['count', ['Blocks', 'Chimeric Blocks']].sum()
        orb.loc['missed blocks', :] = trueBlocks - orb.loc[[c['label'] for _, c in settings['contig_classes'].items() if c['class'] == 'good'], :].sum()
        # select metrics for "recovery analysis"
        orb = orb.loc[[c['label'] for _, c in settings['contig_classes'].items()], :]
        # sort assembler by amount of good contigs
        orb = orb[orb.loc[[c['label'] for _, c in settings['contig_classes'].items() if c['class'] == 'good']].sum().sort_values(ascending=True).index]
        # use pretty label for assembler
        orb = orb.rename(columns=settings['labels']['assemblers'])
    
        # drawing
        # bad contigs
        ax_bad = axes[i // num_columns, (i % num_columns) * 2]
        ax_bad.invert_xaxis()
        orb.loc[[c['label'] for _, c in settings['contig_classes'].items() if c['class'] == 'bad']].T.plot(kind='barh', stacked=True, ax=ax_bad, color={c['label']: c['color'] for _, c in settings['contig_classes'].items()})
        ax_bad.set_xlabel("number contigs")
        ax_top_bad = ax_bad.twiny()
        ax_top_bad.xaxis.set_label_position('top')
        ax_top_bad.set_xticks([])
        ax_top_bad.set_xlabel("bad: artifacts")
        ax_bad.xaxis.set_label_coords(1, -0.08)
        ax_bad.text(-0.5, 1.05, chr(97+i), transform=ax_bad.transAxes, fontsize=16, fontweight='bold',)
    
        # good contigs
        ax_good = axes[i // num_columns, (i % num_columns) * 2  + 1]
        orb.loc[[c['label'] for _, c in settings['contig_classes'].items() if c['class'] != 'bad']].T.plot(kind='barh', stacked=True, ax=ax_good, color={c['label']: c['color'] for _, c in settings['contig_classes'].items()})
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