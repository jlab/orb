from glob import glob
from os.path import join, basename
import yaml
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.colors as mcolors
from matplotlib import patches as mpatches
from matplotlib.lines import Line2D
from matplotlib_venn import venn2
import colorsys
import seaborn as sns
from compile_data import get_environments, getdata_recovery, getdata_gene_recovery, getdata_runtime_memory, getdata_DEgenes, getdata_DEvennOrtho, getdata_DEorthogroups
from tqdm import tqdm


def plot_recovery(fp_orb_basedir, settings, num_columns:int=3, verbose=True):
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
    verbose : boolean
        Report progress on sys.stderr

    Returns
    -------
    plt figure of the multi-pabel plot.
    """
    # which environments to plot and in which order
    environments = get_environments(fp_orb_basedir, settings)
    # load data
    data_recovery = getdata_recovery(fp_orb_basedir, settings, verbose)
    
    fig, axes = plt.subplots(
        int(np.ceil(len(environments) / num_columns)), 
        num_columns * 2, figsize=(2 * num_columns * 4, np.ceil(len(environments) / num_columns) * 5),
        gridspec_kw={"wspace": 0.31, "hspace": 0.3})
    
    for i, environment in tqdm(enumerate(environments), disable=not verbose, desc='Drawing panels for contig recovery plot'):
        orb = data_recovery[data_recovery['environment'] == environment]

        # bad contigs
        ax_bad = axes[i // num_columns, (i % num_columns) * 2]
        ax_bad.invert_xaxis()
        orb.loc[:, [c['label'] for _, c in settings['contig_classes'].items() if c['class'] == 'bad']].plot(kind='barh', stacked=True, ax=ax_bad, color={c['label']: c['color'] for _, c in settings['contig_classes'].items()})
        ax_bad.set_xlabel("number contigs")
        ax_top_bad = ax_bad.twiny()
        ax_top_bad.xaxis.set_label_position('top')
        ax_top_bad.set_xticks([])
        ax_top_bad.set_xlabel("weak")
        ax_bad.xaxis.set_label_coords(1, -0.08)
        ax_bad.text(-0.5, 1.05, chr(97+i), transform=ax_bad.transAxes, fontsize=16, fontweight='bold',)

        # good contigs
        ax_good = axes[i // num_columns, (i % num_columns) * 2  + 1]
        orb.loc[:, [c['label'] for _, c in settings['contig_classes'].items() if c['class'] != 'bad']].plot(kind='barh', stacked=True, ax=ax_good, color={c['label']: c['color'] for _, c in settings['contig_classes'].items()})
        ax_good.set_yticks([])
        ax_good.set_xlabel("robust")
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

def _color_desaturate(color, saturation=0.75):
    h, s, v = colorsys.rgb_to_hsv(*mcolors.to_rgb(color))
    return colorsys.hsv_to_rgb(h, s, v*saturation)

def plot_gene_recovery(fp_orb_basedir, settings, num_columns:int=3, verbose=True):
    # which environments to plot and in which order
    environments = get_environments(fp_orb_basedir, settings)
    # load data
    recovered_genes = getdata_gene_recovery(fp_orb_basedir, settings)

    fig, axes = plt.subplots(
        int(np.ceil(len(environments) / num_columns)), 
        num_columns, figsize=(num_columns * 5, np.ceil(len(environments) / num_columns) * 3),
        gridspec_kw={"wspace": 0.6, "hspace": 0.4}
    )

    palette = {'core': 'gold', 'shared': 'cyan', 'total': 'lightgray'}
    assemblers = [x for x in recovered_genes.index if x not in palette.keys()]
    palette.update({assembler: sns.color_palette()[0] for assembler in assemblers})

    for i, environment in tqdm(enumerate(environments), disable=not verbose, desc='Drawing panels for gene recovery plot'):
        ax = axes[i // num_columns, i % num_columns]

        order = list(recovered_genes.loc[assemblers, environment].sort_values(ascending=False).index) + ['core', 'shared']
        sns.barplot(data=recovered_genes[environment].to_frame().reset_index(), orient='h', y='assembler', hue='assembler', x=environment, ax=ax, palette=palette, order=order)
        ax.axvline(x=recovered_genes.loc['core', environment], color=palette['core'])
        ax.set_ylabel("")
        ax.set_xlabel("number recovered genes")
        ax.set_title(settings['labels']['environments'].get(environment, environment))
        ax.set_xscale('log')
        
        if i+1 == len(environments):
            ax.legend(handles=[
                mpatches.Patch(color=palette[assemblers[0]], label='exclusively recovered'),
                mpatches.Patch(color=_color_desaturate(palette['core'], 0.75), label='recovered by all'),
                mpatches.Patch(color=_color_desaturate(palette['shared'], 0.75), label='recovered by some')],
                    bbox_to_anchor=(-0.1, -0.25), ncols=3)
        
        # panel labels
        ax.text(-0.43, 1.05, chr(97+i), transform=ax.transAxes, fontsize=16, fontweight='bold',)

    return fig


def plotTimeMemory(fp_caviar_basedir:str, settings, verbose=True):
    timemem = getdata_runtime_memory(fp_caviar_basedir, settings)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), gridspec_kw={"wspace": 0.5})
    
    for ax, (pType, factor, label, title) in zip(axes, [
            ('CPU time (seconds)', 3600, 'CPU hours', 'Runtime'), 
            ('Maximum resident set size (kbytes): ', 1024**2, 'RAM in GB', 'Memory footprint')]):
        plotdata = timemem[timemem['type'] == pType].copy()
        plotdata['unit'] = plotdata['value'] / factor
        plotdata = plotdata.sort_values(by='unit')
        sns.boxplot(data=plotdata, x='unit', y='assembler', orient='h', ax=ax, color='lightgray')
        sns.stripplot(data=plotdata, x='unit', y='assembler', orient='h', ax=ax, hue='environment', hue_order=settings['labels']['environments'])
        ax.set_xlabel(label)
        ax.set_ylabel("")
        ax.set_title(title)
        if ax != axes[-1]:
            ax.legend().remove()
        else:
            ax.axvline(x=64, linestyle='-.', color='gray', zorder=-1, label="64 GB laptop")
            ax.legend(bbox_to_anchor=(1.1, -0.15), ncols=7)

    # panel labels
    for i, ax in enumerate(axes):
        ax.text(-0.43, 1.05, chr(97+i), transform=ax.transAxes, fontsize=16, fontweight='bold',)

    return fig


def get_rank_shifts(ranksA: pd.Series, ranksB: pd.Series, mergeAssembler={'idba': ['IDBA-MT', 'IDBA-tran']}):
    """Given two rankings for the same set of assembler, compute the difference in rank positions.
    
    Note: optionally treat different assemblers as identical, e.g. IDBA-MT and IDBA-tran
    """
    def _merge_assembler(ranks, mergeAssembler):
        # map different assemblers to same label, e.g. for IDBA-MT and IDBA-tran, because we don't see differences and want to keep them as one
        ranks.name = 'old_ranks'
        ranks = ranks.to_frame()
        ranks['merged'] = list(map(lambda x: {label: mergelabel for mergelabel, assemblers in mergeAssembler.items() for label in assemblers}.get(x, x), ranks.index))
        return ranks

    def _rerank_nomergedassemblers(ranks):
        # assign new ranks according to the given ones, but assign same rank to same assembler name
        rank = 0
        currAss = None
        newranks = []
        ranks = ranks.sort_values('old_ranks')
        for assembler in ranks['merged'].values:
            if currAss != assembler:
                rank +=1
            currAss = assembler
            newranks.append(rank)
        ranks['merged_ranks'] = newranks
        return ranks['merged_ranks']

    cmp = pd.concat([_rerank_nomergedassemblers(_merge_assembler(ranksA, mergeAssembler)).rename('rank_reference'),
                     _rerank_nomergedassemblers(_merge_assembler(ranksB, mergeAssembler)).rename('rank_other')], axis=1)
    # should an assembler be missing in one of the ranks, assign worst+1 rank to it
    for col in cmp.columns:
        cmp[col] = cmp[col].fillna(cmp[col].max() + 1)
    cmp['shift'] = cmp.iloc[:, 0] - cmp.iloc[:, 1]
    
    def _create_label(row):
        if row['shift'] > 0:
            return '%s ⬆+%i' % (row.name, row['shift'])
        elif row['shift'] < 0:
            return '%s ⬇%i' % (row.name, row['shift'])
        else:
            return '%s %i' % (row.name, row['shift'])
    cmp['label'] = cmp.apply(_create_label, axis=1)

    def _create_color(row):
        if row['shift'] > 0:
            return 'darkgreen'
        elif row['shift'] < 0:
            return 'darkorange'
        else:
            return 'black'
    cmp['color'] = cmp.apply(_create_color, axis=1)
    
    return  cmp

# def plot_DEgenes(fp_orb_basedir, settings, forOrthogroups=False, fp_marbel_basedir:str=None, fp_ogtruth_basedir:str=None, num_columns:int=3, verbose=True):
#     # which environments to plot and in which order
#     environments = get_environments(fp_orb_basedir, settings)
    
#     # load data
#     if forOrthogroups is False:
#         degenes, _ = getdata_DEgenes(fp_orb_basedir, settings, verbose)
#     else:
#         degenes, _ = getdata_DEorthogroups(fp_orb_basedir, fp_marbel_basedir, fp_ogtruth_basedir, settings)
#     recovery = getdata_recovery(fp_orb_basedir, settings, verbose)

#     fig, axes = plt.subplots(
#         int(np.ceil(len(environments) / num_columns)), 
#         num_columns, figsize=(num_columns * 7, np.ceil(len(environments) / num_columns) * 4),
#         gridspec_kw={"hspace": 0.31, "wspace": 0.5})
#     BARWIDTH=0.8

#     palette = {'True Positive': '#238cc3',
#             'False Positive': '#5d5e60',
#             'False Negative': '#c06364'}
#     for i, environment in tqdm(enumerate(environments), disable=not verbose, desc='Drawing panels for DE plot'):
#         ax = axes[i // num_columns, i % num_columns]
        
#         order = list(reversed(
#             pd.pivot_table(data=degenes.loc[environment, :], index='assembler', columns='class', values='num_genes', aggfunc="sum").sort_values(
#                 by=       ['True Positive', 'False Positive', 'False Negative'],
#                 ascending=[False,           True,             True]).index))
#         for y, assembler in enumerate(order):
#             cls = 'True Positive'
#             pos_TP = degenes.loc[environment, assembler, cls]['num_genes']
#             ax.add_patch(plt.Rectangle((0, y), pos_TP, BARWIDTH, facecolor=palette[cls], edgecolor='white', linewidth=1, label=cls if y == 0 else None))
        
#             cls = 'False Positive'
#             pos_FP = degenes.loc[environment, assembler, cls]['num_genes']
#             ax.add_patch(plt.Rectangle((pos_TP, y), pos_FP, BARWIDTH, facecolor=palette[cls], edgecolor='white', linewidth=1, label=cls if y == 0 else None))
        
#             cls = 'False Negative'
#             pos_FN = degenes.loc[environment, assembler, cls]['num_genes']
#             ax.add_patch(plt.Rectangle((0, y), -1 * pos_FN, BARWIDTH, facecolor=palette[cls], edgecolor='white', linewidth=1, label=cls if y == 0 else None))
#         ax.set_ylim((-1 * (1 - BARWIDTH), len(order)))
        
#         ranks = get_rank_shifts(recovery[recovery['environment'] == environment]['recovery_rank'],
#                                 pd.Series(index=order, data=reversed(range(1, len(order) + 1))))
#         ax.set_yticks(list(map(lambda x: x + BARWIDTH/2, range(len(order)))), ranks.loc[order, 'label'].values)
#         for tick_label, color in zip(ax.get_yticklabels(), ranks.loc[order, 'color'].values):
#             tick_label.set_color(color)
                
#         ax.set_title(settings['labels']['environments'].get(environment, environment))
#         ax.set_xlabel('number genes')
        
#         num_positives = degenes.loc[environment, :].reset_index().set_index('truth').loc[True].groupby('assembler')['num_genes'].sum().iloc[0]
#         ax.axvline(x=num_positives, color=palette['True Positive'], label='Positive')
#         maxX = (degenes.loc[environment, :, 'True Positive']['num_genes'] + degenes.loc[environment, :, 'False Positive']['num_genes']).max()
#         for (ex_env, ex_ass) in [('seawater', 'Trinity'),
#                                 ('freshwater', 'Trinity'),
#                                 ('healthy_gut', 'dbg')]:
#             if environment == ex_env:
#                 pdata = degenes.loc[environment, [ass for ass in degenes.loc[environment, :].index.levels[0] if ass != ex_ass], :]
#                 maxX = (pdata.loc[environment, :, 'True Positive']['num_genes'] + pdata.loc[environment, :, 'False Positive']['num_genes']).max()

#                 ax.text(max(1, num_positives, maxX), list(order).index(ex_ass) + 0.35, '%i' % degenes.loc[environment, ex_ass, 'False Positive']['num_genes'],
#                         verticalalignment='center', horizontalalignment='right', color='white')

#         # panel labels
#         ax.text(-0.43, 1.05, chr(97+i), transform=ax.transAxes, fontsize=16, fontweight='bold',)

#         ax.set_xlim((
#             -1.1 * max(1, degenes.loc[environment, :, 'False Negative']['num_genes'].max()),
#             1.1 * max(1, num_positives, maxX)))

#         if i+1 == len(environments):
#             ax.legend(handles=[mpatches.Patch(color=palette[cls], label=cls) for cls in palette.keys()] + \
#                               [Line2D([0], [0], color=palette['True Positive'], lw=2, label='Positive')],
#                       #title='Category',
#                       bbox_to_anchor=(-0.4, -0.25),
#                       ncols=4)
#     return fig

def plot_DEgenes(fp_orb_basedir, settings, forOrthogroups=False, fp_marbel_basedir:str=None, fp_ogtruth_basedir:str=None, num_columns:int=3, verbose=True):
    # which environments to plot and in which order
    environments = get_environments(fp_orb_basedir, settings)
    
    # load data
    DEfeatures, _ = getdata_DEgenes(fp_orb_basedir, settings, verbose)
    if forOrthogroups is False:
        rank_data = getdata_recovery(fp_orb_basedir, settings, verbose)
    else:
        rank_data = DEfeatures.loc[:, :, 'True Positive']['rank'].reset_index().rename(columns={'rank': 'recovery_rank'}).set_index('assembler').copy()
        DEfeatures, _ = getdata_DEorthogroups(fp_orb_basedir, fp_marbel_basedir, fp_ogtruth_basedir, settings)

    fig, axes = plt.subplots(
        int(np.ceil(len(environments) / num_columns)), 
        num_columns, figsize=(num_columns * 7, np.ceil(len(environments) / num_columns) * 4),
        gridspec_kw={"hspace": 0.31, "wspace": 0.5})
    BARWIDTH=0.8

    palette = {'True Positive': '#238cc3',
               'False Positive': '#5d5e60',
               'False Negative': '#c06364'}
    for i, environment in tqdm(enumerate(environments), disable=not verbose, desc='Drawing panels for DE plot'):
        ax = axes[i // num_columns, i % num_columns]
        
        order = list(DEfeatures.loc[environment, :].reset_index().groupby('assembler').head(1).set_index('assembler')['rank'].sort_values(ascending=False).index)
        for y, assembler in enumerate(order):
            cls = 'True Positive'
            pos_TP = DEfeatures.loc[environment, assembler, cls]['num_genes']
            ax.add_patch(plt.Rectangle((0, y), pos_TP, BARWIDTH, facecolor=palette[cls], edgecolor='white', linewidth=1, label=cls if y == 0 else None))
        
            cls = 'False Positive'
            pos_FP = DEfeatures.loc[environment, assembler, cls]['num_genes']
            ax.add_patch(plt.Rectangle((pos_TP, y), pos_FP, BARWIDTH, facecolor=palette[cls], edgecolor='white', linewidth=1, label=cls if y == 0 else None))
        
            cls = 'False Negative'
            pos_FN = DEfeatures.loc[environment, assembler, cls]['num_genes']
            ax.add_patch(plt.Rectangle((0, y), -1 * pos_FN, BARWIDTH, facecolor=palette[cls], edgecolor='white', linewidth=1, label=cls if y == 0 else None))
        ax.set_ylim((-1 * (1 - BARWIDTH), len(order)))
        
        ranks = get_rank_shifts(rank_data[rank_data['environment'] == environment]['recovery_rank'],
                                pd.Series(index=order, data=reversed(range(1, len(order) + 1))))
        ax.set_yticks(list(map(lambda x: x + BARWIDTH/2, range(len(order)))), ranks.loc[order, 'label'].values)
        for tick_label, color in zip(ax.get_yticklabels(), ranks.loc[order, 'color'].values):
            tick_label.set_color(color)
                
        ax.set_title(settings['labels']['environments'].get(environment, environment))
        ax.set_xlabel('number genes')
        
        num_positives = DEfeatures.loc[environment, :].reset_index().set_index('truth').loc[True].groupby('assembler')['num_genes'].sum().iloc[0]
        ax.axvline(x=num_positives, color=palette['True Positive'], label='Positive')
        maxX = (DEfeatures.loc[environment, :, 'True Positive']['num_genes'] + DEfeatures.loc[environment, :, 'False Positive']['num_genes']).max()
        for (ex_forOrtho, ex_env, ex_ass) in [(False, 'seawater', 'Trinity'),
                                              (False, 'freshwater', 'Trinity'),
                                              (False, 'healthy_gut', 'dbg'),
                                             ]:
            if (ex_forOrtho == forOrthogroups) and (environment == ex_env):
                pdata = DEfeatures.loc[environment, [ass for ass in DEfeatures.loc[environment, :].index.levels[0] if ass != ex_ass], :]
                maxX = (pdata.loc[environment, :, 'True Positive']['num_genes'] + pdata.loc[environment, :, 'False Positive']['num_genes']).max()

                ax.text(max(1, num_positives, maxX), list(order).index(ex_ass) + 0.35, '%i' % DEfeatures.loc[environment, ex_ass, 'False Positive']['num_genes'],
                        verticalalignment='center', horizontalalignment='right', color='white')

        # panel labels
        ax.text(-0.43, 1.05, chr(97+i), transform=ax.transAxes, fontsize=16, fontweight='bold',)

        ax.set_xlim((
            -1.1 * max(1, DEfeatures.loc[environment, :, 'False Negative']['num_genes'].max()),
             1.1 * max(1, num_positives, maxX)))

        if i+1 == len(environments):
            ax.legend(handles=[mpatches.Patch(color=palette[cls], label=cls) for cls in palette.keys()] + \
                              [Line2D([0], [0], color=palette['True Positive'], lw=2, label='Positive')],
                      #title='Category',
                      bbox_to_anchor=(-0.4, -0.25),
                      ncols=4)
    return fig


def plot_DEvennOrtho(fp_orb_basedir:str, fp_ogtruth_basedir:str, fp_marbel_basedir:str, settings, num_columns:int=3, verbose=True):
    # which environments to plot and in which order
    environments = get_environments(fp_orb_basedir, settings)

    # get data
    truth = getdata_DEvennOrtho(fp_orb_basedir, fp_ogtruth_basedir, fp_marbel_basedir, settings)

    fig, axes = plt.subplots(
        3,
        len(environments), figsize=(len(environments) * 3, len(environments)/2 * 3),
        #gridspec_kw={"wspace": 0.6, "hspace": 0.4}
    )
    palette = {'DE orthogroups': 'blue',
               'DE genes contained in orthogroups': 'orange'}
    labels = ["", ""]

    reordered_environments = []
    for start in range(num_columns):
        reordered_environments.extend(environments[start::num_columns])

    for i, environment in tqdm(enumerate(reordered_environments), disable=not verbose, desc='Drawing panels for DE Venn diagrams'):
        ax = axes[0][i]
        venn2([set(truth[environment][truth[environment]['DEorthogroup']]['orthogroup'].values),
               set(truth[environment][truth[environment]['DEgene']]['orthogroup'].values)],
              labels,
              ax=ax, set_colors=(palette['DE orthogroups'], palette['DE genes contained in orthogroups']))
        ax.set_title(settings['labels']['environments'].get(environment, environment))
        ax.set_ylabel("all orthogroups")            

        ax = axes[1][i]
        venn2([set(truth[environment][truth[environment]['DEorthogroup'] & (truth[environment]['OGsize'] == 'multi_gene_OG')]['orthogroup'].values),
               set(truth[environment][truth[environment]['DEgene'] & (truth[environment]['OGsize'] == 'multi_gene_OG')]['orthogroup'].values)],
              labels,
              ax=ax, set_colors=(palette['DE orthogroups'], palette['DE genes contained in orthogroups']))
        ax.set_ylabel("only multi gene orthogroups")

        ax = axes[2][i]
        venn2([set(truth[environment][truth[environment]['DEorthogroup'] & (truth[environment]['OGsize'] == 'single_gene_OG')]['orthogroup'].values),
               set(truth[environment][truth[environment]['DEgene'] & (truth[environment]['OGsize'] == 'single_gene_OG')]['orthogroup'].values)],
              labels,
              ax=ax, set_colors=(palette['DE orthogroups'], palette['DE genes contained in orthogroups']))
        ax.set_ylabel("only single gene orthogroups")

        if i == 0:
            for row in range(len(axes)):
                axes[row][i].axison = True
                for p in ['left', 'right', 'bottom', 'top']:
                    axes[row][i].spines[p].set_color('white')
    axes[2][0].legend(handles=[mpatches.Patch(label=grp, facecolor=color, alpha=0.4) for (grp, color) in palette.items()],
                      bbox_to_anchor=(4.8, -0.),
                      ncols=2)
    
    return fig