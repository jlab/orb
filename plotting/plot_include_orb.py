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
import matplotlib.ticker as ticker
import colorsys
import seaborn as sns
from compile_data_orb import (
    get_environments,
    getdata_DEgenes_mod,
    getdata_recovery,
    getdata_gene_recovery,
    getdata_runtime_memory,
    getdata_DEvennOrtho,
    getdata_DEorthogroups_mod,
    get_recovered_contigs,
    getdata_rnaquast,
    getdata_block_recovery,
    get_contigs,
    filter_for_assembler_with_ns,
    getdata_detonate,
)
from tqdm import tqdm
from scipy.cluster.hierarchy import linkage, leaves_list, dendrogram
from scipy.spatial.distance import pdist, squareform
from skbio.stats.distance import DistanceMatrix
from statannotations.Annotator import Annotator


def legendentries_col2rowordering(handles, labels, num_rows, dummiesatpositions=[]):
    entries = list(zip(handles, labels))
    num_columns = int(np.ceil(len(entries) / num_rows))
    for pos in dummiesatpositions:
        entries.insert(pos, (Line2D([], [], linestyle="none"), ""))
    reordered_entries = []
    for i in range(num_columns):
        reordered_entries.extend(entries[i::num_columns])
    return zip(*reordered_entries)


def plot_recovery(fp_orb_basedir, settings, num_columns: int = 3, verbose=True, detail_view=False, report_percent=True):
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
    report_percent : boolean
        If true, legend items also mention the fraction of contigs classified as
        such across all environments / all assemblers. This gives a feeling for
        importance of a class.

    Returns
    -------
    plt figure of the multi-pabel plot.
    """
    # which environments to plot and in which order
    environments = get_environments(fp_orb_basedir, settings)
    # load dataf
    if detail_view:
        data_recovery = getdata_recovery(fp_orb_basedir, settings, detail_view=detail_view, verbose=verbose)
    else:
        data_recovery = getdata_recovery(fp_orb_basedir, settings, verbose=verbose)

    fig, axes = plt.subplots(
        int(np.ceil(len(environments) / num_columns)),
        num_columns * 2,
        figsize=(2 * num_columns * 4, np.ceil(len(environments) / num_columns) * 3),
        gridspec_kw={"wspace": 0.31, "hspace": 0.45},
    )

    if detail_view:
        view_class = "contig_classes"
    else:
        view_class = "simplelabels"

    for i, environment in tqdm(
        enumerate(environments), disable=not verbose, desc="Drawing panels for contig recovery plot"
    ):
        orb = data_recovery[data_recovery["environment"] == environment]

        # bad contigs
        ax_bad = axes[i // num_columns, (i % num_columns) * 2]
        ax_bad.invert_xaxis()
        orb.loc[:, [c["label"] for _, c in settings[view_class].items() if c["class"] == "bad"]].plot(
            kind="barh",
            stacked=True,
            ax=ax_bad,
            color={c["label"]: c["color"] for _, c in settings["contig_classes"].items()},
        )

        ax_bad.set_xlabel("number contigs")
        ax_top_bad = ax_bad.twiny()
        ax_top_bad.xaxis.set_label_position("top")
        ax_top_bad.set_xticks([])
        ax_top_bad.set_xlabel(settings["labels"]["recovery_plot"]["bad"])
        ax_bad.xaxis.set_label_coords(1, -0.23)
        ax_bad.text(
            -0.5,
            1.00,
            chr(97 + i),
            transform=ax_bad.transAxes,
            fontsize=16,
            fontweight="bold",
        )
        ax_bad.tick_params(axis="x", labelrotation=20)

        # good contigs
        ax_good = axes[i // num_columns, (i % num_columns) * 2 + 1]
        orb.loc[
            :,
            [
                c["label"]
                for _, c in settings[view_class].items()
                if (c["class"] == "good") or (c["class"] == "neutral")
            ],
        ].plot(
            kind="barh",
            stacked=True,
            ax=ax_good,
            color={c["label"]: c["color"] for _, c in settings["contig_classes"].items()},
        )
        ax_good.set_yticks([])
        ax_good.set_xlabel(settings["labels"]["recovery_plot"]["good"])
        ax_good.xaxis.set_label_position("top")
        ax_good.tick_params(axis="x", labelrotation=20)

        # concat right (=good) axis directly adjacent to left (=bad) axis
        ax_good.set_position(
            [
                ax_bad.get_position().x1,
                ax_good.get_position().y0,
                ax_good.get_position().width,
                ax_good.get_position().height,
            ]
        )
        # ax_bad.set_title(settings['labels']['environments'].get(environment, environment), loc='right', horizontalalignment='center')
        ax_bad.text(
            0.06,
            0.96,
            settings["labels"]["environments"].get(environment, environment),
            transform=ax_bad.transAxes,
            fontsize=12,
            fontweight="bold",
            va="top",
        )

        # one legend for all panels
        if i + 1 == len(environments):
            handles = ax_good.get_legend_handles_labels()[0] + list(reversed(ax_bad.get_legend_handles_labels()[0]))
            labels = ax_good.get_legend_handles_labels()[1] + list(reversed(ax_bad.get_legend_handles_labels()[1]))
            # if seperate position for label is provided in style yaml use it
            first_elem = list(settings[view_class].keys())[0]
            if "legend_position" in settings[view_class][first_elem]:
                sorted_labels = []
                for i in range(len(settings[view_class])):
                    label = [c["label"] for c in settings[view_class].values() if c["legend_position"] == i]
                    if len(label) != 1:
                        print("Fatal error: legend positions do not match ordered numbered of labels, check style.yaml")
                        return
                    sorted_labels.append(label[0])
                sorted_handels = []
                for label in sorted_labels:
                    sorted_handels.append(handles[labels.index(label)])
                handles = sorted_handels
                labels = sorted_labels

            if report_percent:
                # print percent of contigs in each class in legend
                # compute ratio classes across all environments and assembler (except missed blocks)
                denominators_noneutral = data_recovery[
                    [
                        c["label"]
                        for _, c in settings["contig_classes"].items()
                        if c["class"] != "neutral"
                        if c["label"] in data_recovery.columns
                    ]
                ].sum()
                percent = denominators_noneutral / denominators_noneutral.sum()
                # compute ratio of missed blocks: neutral / (good + neutral) blocks
                denominators_rhs = data_recovery[
                    [
                        c["label"]
                        for _, c in settings["contig_classes"].items()
                        if c["class"] != "bad"
                        if c["label"] in data_recovery.columns
                    ]
                ].sum()
                percent = pd.concat(
                    [
                        percent,
                        (
                            denominators_rhs[
                                [c["label"] for _, c in settings["contig_classes"].items() if c["class"] == "neutral"]
                            ]
                            / denominators_rhs.sum()
                        ),
                    ]
                ).apply(lambda x: "%.2f%%" % (x * 100))
                # updates legend labels
                labels = ["%s (%s)" % (l, percent.loc[l]) if l in percent.index else l for l in labels]

            # reorder entries in legend
            xoffset = 0
            if detail_view is False:
                handles, labels = legendentries_col2rowordering(handles, labels, 2, dummiesatpositions=[4])
            else:
                xoffset = 0.2
            ax_good.legend(handles, labels, ncol=5, bbox_to_anchor=(1 + xoffset, -0.3))
        else:
            ax_good.legend().remove()
        ax_bad.legend().remove()

        # TODO: print DF for check with orb base dir
    return fig


# # above function should be called as following:
# # 1) load plotting settings
# with open("/homes/sjanssen/Git/jlab/orb/plotting/style.yaml", "r") as f:
#    settings = yaml.safe_load(f)
# # 2) then call plot_recovery with the filepath to Orb's base dir and the loaded settings
# _ = plot_recovery('/vol/jlab/tlin/all_project/nf_results/orb', settings)


def _color_desaturate(color, saturation=0.75):
    h, s, v = colorsys.rgb_to_hsv(*mcolors.to_rgb(color))
    return colorsys.hsv_to_rgb(h, s, v * saturation)


def plot_gene_recovery(fp_orb_basedir, settings, num_columns: int = 3, verbose=True, recovered_genes=None):
    # which environments to plot and in which order
    environments = get_environments(fp_orb_basedir, settings)
    # load data
    if recovered_genes is None:
        recovered_genes = getdata_gene_recovery(fp_orb_basedir, settings)

    fig, axes = plt.subplots(
        int(np.ceil(len(environments) / num_columns)),
        num_columns,
        figsize=(num_columns * 5, np.ceil(len(environments) / num_columns) * 2.5),
        gridspec_kw={"wspace": 0.6, "hspace": 0.4},
    )

    palette = {"core": "gold", "shared": "cyan", "total": "lightgray"}
    assemblers = [x for x in recovered_genes.index if x not in palette.keys()]
    palette.update({assembler: sns.color_palette()[0] for assembler in assemblers})

    for i, environment in tqdm(
        enumerate(environments), disable=not verbose, desc="Drawing panels for gene recovery plot"
    ):
        ax = axes[i // num_columns, i % num_columns]

        order = list(recovered_genes.loc[assemblers, environment].sort_values(ascending=False).index) + [
            "core",
            "shared",
        ]
        sns.barplot(
            data=recovered_genes[environment].to_frame().reset_index(),
            orient="h",
            y="assembler",
            hue="assembler",
            x=environment,
            ax=ax,
            palette=palette,
            order=order,
        )
        ax.axvline(x=recovered_genes.loc["core", environment], color=palette["core"])
        ax.set_ylabel("")
        ax.set_xlabel("number recovered genes")
        # ax.set_title(settings['labels']['environments'].get(environment, environment))
        ax.text(
            0.95,
            0.8,
            settings["labels"]["environments"].get(environment, environment),
            transform=ax.transAxes,
            fontsize=12,
            fontweight="bold",
            va="center",
            ha="right",
        )
        ax.set_xscale("log")

        if i + 1 == len(environments):
            ax.legend(
                handles=[
                    mpatches.Patch(color=palette[assemblers[0]], label="exclusively recovered"),
                    mpatches.Patch(color=_color_desaturate(palette["core"], 0.75), label="recovered by all"),
                    mpatches.Patch(color=_color_desaturate(palette["shared"], 0.75), label="recovered by at least two"),
                ],
                bbox_to_anchor=(0.5, -0.35),
                ncols=3,
            )

        # panel labels
        ax.text(
            -0.48,
            0.95,
            chr(97 + i),
            transform=ax.transAxes,
            fontsize=16,
            fontweight="bold",
        )

    return fig


def plotTimeMemory(fp_caviar_basedir: str, settings, verbose=True):
    timemem = getdata_runtime_memory(fp_caviar_basedir, settings)

    fig, axes = plt.subplots(1, 2, figsize=(10, 2.5), gridspec_kw={"wspace": 0.5})

    for ax, (pType, factor, label, title) in zip(
        axes,
        [
            ("CPU time (seconds)", 3600, "CPU hours", "Runtime"),
            ("Maximum resident set size (kbytes): ", 1024**2, "RAM in GB", "Memory footprint"),
        ],
    ):
        plotdata = timemem[timemem["type"] == pType].copy()
        plotdata["unit"] = plotdata["value"] / factor
        plotdata = plotdata.sort_values(by="unit")
        sns.boxplot(data=plotdata, x="unit", y="assembler", orient="h", ax=ax, color="lightgray")
        sns.stripplot(
            data=plotdata,
            x="unit",
            y="assembler",
            orient="h",
            ax=ax,
            hue="environment",
            hue_order=list(settings["labels"]["environments"].values()),
        )
        ax.set_xlabel(label)
        ax.set_ylabel("")
        ax.set_title(title)
        if ax != axes[-1]:
            ax.axvline(x=4 * 24, linestyle=":", color="gray", zorder=-1, label="4 CPU days")
            rt_handles, rt_labels = ax.get_legend_handles_labels()
            ax.legend().remove()
        else:
            ax.axvline(x=64, linestyle="-.", color="gray", zorder=-1, label="64 GB laptop")
            mem_handles, mem_labels = ax.get_legend_handles_labels()
            mem_handles.insert(len(mem_handles) - 1, rt_handles[-1])
            mem_labels.insert(len(mem_labels) - 1, rt_labels[-1])
            ax.legend(mem_handles, mem_labels, bbox_to_anchor=(0.5, -0.25), ncols=4)

    # panel labels
    for i, ax in enumerate(axes):
        ax.text(
            -0.43,
            1.00,
            chr(97 + i),
            transform=ax.transAxes,
            fontsize=16,
            fontweight="bold",
        )

    return fig


def get_rank_shifts(ranksA: pd.Series, ranksB: pd.Series, mergeAssembler={"idba": ["IDBA-MT", "IDBA-tran"]}):
    """Given two rankings for the same set of assembler, compute the difference in rank positions.

    Note: optionally treat different assemblers as identical, e.g. IDBA-MT and IDBA-tran
    """

    def _merge_assembler(ranks, mergeAssembler):
        # map different assemblers to same label, e.g. for IDBA-MT and IDBA-tran, because we don't see differences and want to keep them as one
        ranks.name = "old_ranks"
        ranks = ranks.to_frame()
        ranks["merged"] = list(
            map(
                lambda x: {
                    label: mergelabel for mergelabel, assemblers in mergeAssembler.items() for label in assemblers
                }.get(x, x),
                ranks.index,
            )
        )
        return ranks

    def _rerank_nomergedassemblers(ranks):
        # assign new ranks according to the given ones, but assign same rank to same assembler name
        rank = 0
        currAss = None
        newranks = []
        ranks = ranks.sort_values("old_ranks")
        for assembler in ranks["merged"].values:
            if currAss != assembler:
                rank += 1
            currAss = assembler
            newranks.append(rank)
        ranks["merged_ranks"] = newranks
        return ranks["merged_ranks"]

    cmp = pd.concat(
        [
            _rerank_nomergedassemblers(_merge_assembler(ranksA, mergeAssembler)).rename("rank_reference"),
            _rerank_nomergedassemblers(_merge_assembler(ranksB, mergeAssembler)).rename("rank_other"),
        ],
        axis=1,
    )
    # should an assembler be missing in one of the ranks, assign worst+1 rank to it
    for col in cmp.columns:
        cmp[col] = cmp[col].fillna(cmp[col].max() + 1)
    cmp["shift"] = cmp.iloc[:, 0] - cmp.iloc[:, 1]

    def _create_label(row):
        if row["shift"] > 0:
            return "%s ⬆+%i" % (row.name, row["shift"])
        elif row["shift"] < 0:
            return "%s ⬇%i" % (row.name, row["shift"])
        else:
            return "%s %i" % (row.name, row["shift"])

    cmp["label"] = cmp.apply(_create_label, axis=1)

    def _create_color(row):
        if row["shift"] > 0:
            return "darkgreen"
        elif row["shift"] < 0:
            return "darkorange"
        else:
            return "black"

    cmp["color"] = cmp.apply(_create_color, axis=1)

    return cmp


def plot_DEgenes(fp_orb_basedir, settings, forOrthogroups=False, num_columns: int = 3, verbose=True, sortF1=False):
    # which environments to plot and in which order
    environments = get_environments(fp_orb_basedir, settings)

    # load data
    DEfeatures = getdata_DEgenes_mod(fp_orb_basedir, settings, verbose, sortF1=sortF1)
    if forOrthogroups is False:
        rank_data = getdata_recovery(fp_orb_basedir, settings, verbose=verbose)
    else:
        rank_data = (
            DEfeatures.loc[:, :, "True Positive"]["rank"]
            .reset_index()
            .rename(columns={"rank": "recovery_rank"})
            .set_index("assembler")
            .copy()
        )
        DEfeatures = getdata_DEorthogroups_mod(fp_orb_basedir, settings, sortF1=sortF1)

    fig, axes = plt.subplots(
        int(np.ceil(len(environments) / num_columns)),
        num_columns,
        figsize=(num_columns * 7, np.ceil(len(environments) / num_columns) * 2.5),
        gridspec_kw={"hspace": 0.5, "wspace": 0.5},
    )
    BARWIDTH = 0.8

    palette = {"True Positive": "#238cc3", "False Positive": "#5d5e60", "False Negative": "#c06364"}
    for i, environment in tqdm(enumerate(environments), disable=not verbose, desc="Drawing panels for DE plot"):
        ax = axes[i // num_columns, i % num_columns]

        order = list(
            DEfeatures.loc[environment, :]
            .reset_index()
            .groupby("assembler")
            .head(1)
            .set_index("assembler")["rank"]
            .sort_values(ascending=False)
            .index
        )
        for y, assembler in enumerate(order):
            cls = "True Positive"
            pos_TP = DEfeatures.loc[environment, assembler, cls]["num_genes"]
            ax.add_patch(
                plt.Rectangle(
                    (0, y),
                    pos_TP,
                    BARWIDTH,
                    facecolor=palette[cls],
                    edgecolor="white",
                    linewidth=1,
                    label=cls if y == 0 else None,
                )
            )

            cls = "False Positive"
            pos_FP = DEfeatures.loc[environment, assembler, cls]["num_genes"]
            ax.add_patch(
                plt.Rectangle(
                    (pos_TP, y),
                    pos_FP,
                    BARWIDTH,
                    facecolor=palette[cls],
                    edgecolor="white",
                    linewidth=1,
                    label=cls if y == 0 else None,
                )
            )

            cls = "False Negative"
            pos_FN = DEfeatures.loc[environment, assembler, cls]["num_genes"]
            ax.add_patch(
                plt.Rectangle(
                    (0, y),
                    -1 * pos_FN,
                    BARWIDTH,
                    facecolor=palette[cls],
                    edgecolor="white",
                    linewidth=1,
                    label=cls if y == 0 else None,
                )
            )
        ax.set_ylim((-1 * (1 - BARWIDTH), len(order)))

        ranks = get_rank_shifts(
            rank_data[rank_data["environment"] == environment]["recovery_rank"],
            pd.Series(index=order, data=reversed(range(1, len(order) + 1))),
        )
        ax.set_yticks(list(map(lambda x: x + BARWIDTH / 2, range(len(order)))), ranks.loc[order, "label"].values)
        for tick_label, color in zip(ax.get_yticklabels(), ranks.loc[order, "color"].values):
            tick_label.set_color(color)

        ax.set_title(settings["labels"]["environments"].get(environment, environment))
        if forOrthogroups:
            ax.set_xlabel("number orthologous groups")
        else:
            ax.set_xlabel("number genes")

        num_positives = (
            DEfeatures.loc[environment, :]
            .reset_index()
            .set_index("truth")
            .loc[True]
            .groupby("assembler")["num_genes"]
            .sum()
            .iloc[0]
        )
        ax.axvline(x=num_positives, color=palette["True Positive"], label="Positive")
        maxX = (
            DEfeatures.loc[environment, :, "True Positive"]["num_genes"]
            + DEfeatures.loc[environment, :, "False Positive"]["num_genes"]
        ).max()

        if forOrthogroups:
            shortening_dict = settings["de_image"]["orthogroup_assembler_shortening"]
        else:
            shortening_dict = settings["de_image"]["gene_assembler_shortening"]
        for ex_env, ex_ass in shortening_dict.items():
            if environment == ex_env:
                pdata = DEfeatures.loc[
                    environment, [ass for ass in DEfeatures.loc[environment, :].index.levels[0] if ass != ex_ass], :
                ]
                maxX = (
                    pdata.loc[environment, :, "True Positive"]["num_genes"]
                    + pdata.loc[environment, :, "False Positive"]["num_genes"]
                ).max()
                ax.text(
                    max(1, num_positives, maxX),
                    list(order).index(ex_ass) + 0.35,
                    "%i" % DEfeatures.loc[environment, ex_ass, "False Positive"]["num_genes"],
                    verticalalignment="center",
                    horizontalalignment="right",
                    color="white",
                )

        # panel labels
        ax.text(
            -0.43,
            1.05,
            chr(97 + i),
            transform=ax.transAxes,
            fontsize=16,
            fontweight="bold",
        )

        ax.set_xlim(
            (
                -1.1 * max(1, DEfeatures.loc[environment, :, "False Negative"]["num_genes"].max()),
                1.1 * max(1, num_positives, maxX),
            )
        )

        if i + 1 == len(environments):
            ax.legend(
                handles=[mpatches.Patch(color=palette[cls], label=cls) for cls in palette.keys()]
                + [Line2D([0], [0], color=palette["True Positive"], lw=2, label="Positive")],
                # title='Category',
                bbox_to_anchor=(0.3, -0.35),
                ncols=4,
            )
    return fig


def plot_DEvennOrtho(fp_orb_basedir: str, fp_marbel_basedir: str, settings, verbose=True):
    # which environments to plot and in which order
    environments = get_environments(fp_orb_basedir, settings)

    # get data
    truth = getdata_DEvennOrtho(fp_orb_basedir, fp_marbel_basedir, settings)

    fig, axes = plt.subplots(
        3,
        len(environments),
        figsize=(len(environments) * 3, len(environments) / 2 * 3),
        # gridspec_kw={"wspace": 0.6, "hspace": 0.4}
    )
    palette = {"DE orthologous groups": "blue", "DE genes contained in orthologous groups": "orange"}
    labels = ["", ""]

    for i, environment in tqdm(
        enumerate(environments), disable=not verbose, desc="Drawing panels for DE Venn diagrams"
    ):
        ax = axes[0][i]
        venn2(
            [
                set(truth[environment][truth[environment]["DEorthogroup"]]["orthogroup"].values),
                set(truth[environment][truth[environment]["DEgene"]]["orthogroup"].values),
            ],
            labels,
            ax=ax,
            set_colors=(palette["DE orthologous groups"], palette["DE genes contained in orthologous groups"]),
        )
        ax.set_title(settings["labels"]["environments"].get(environment, environment))
        ax.set_ylabel(r"$\mathbf{all}$" + "\northologous groups")

        ax = axes[1][i]
        venn2(
            [
                set(
                    truth[environment][
                        truth[environment]["DEorthogroup"] & (truth[environment]["OGsize"] == "multi_gene_OG")
                    ]["orthogroup"].values
                ),
                set(
                    truth[environment][
                        truth[environment]["DEgene"] & (truth[environment]["OGsize"] == "multi_gene_OG")
                    ]["orthogroup"].values
                ),
            ],
            labels,
            ax=ax,
            set_colors=(palette["DE orthologous groups"], palette["DE genes contained in orthologous groups"]),
        )
        ax.set_ylabel("only " + r"$\mathbf{multi}$" + " genes\northologous groups")

        ax = axes[2][i]
        venn2(
            [
                set(
                    truth[environment][
                        truth[environment]["DEorthogroup"] & (truth[environment]["OGsize"] == "single_gene_OG")
                    ]["orthogroup"].values
                ),
                set(
                    truth[environment][
                        truth[environment]["DEgene"] & (truth[environment]["OGsize"] == "single_gene_OG")
                    ]["orthogroup"].values
                ),
            ],
            labels,
            ax=ax,
            set_colors=(palette["DE orthologous groups"], palette["DE genes contained in orthologous groups"]),
        )
        ax.set_ylabel("only " + r"$\mathbf{single}$" + " gene\northologous groups")

        if i == 0:
            for row in range(len(axes)):
                axes[row][i].axison = True
                for p in ["left", "right", "bottom", "top"]:
                    axes[row][i].spines[p].set_color("white")
    axes[2][0].legend(
        handles=[mpatches.Patch(label=grp, facecolor=color, alpha=0.4) for (grp, color) in palette.items()],
        bbox_to_anchor=(4.8, -0.0),
        ncols=2,
    )

    return fig


def plot_heatmap(fp_orb_basedir: str, settings, num_columns: int = 3, verbose=True):
    # which environments to plot and in which order
    environments = get_environments(fp_orb_basedir, settings)
    # re-order such that paired-environments are on top of each other
    environments = environments[::2] + environments[1::2]

    # re-order environments such that two additional "environments"
    # are spiked in for the color map and the combined environment
    def _spikein(environments, num_cols=3, spikeelements=["colormap", "all six environments"]):
        chunks = [environments[i : i + num_columns] for i in range(0, len(environments), num_columns)]
        reordered = []
        for i, spike in enumerate(spikeelements):
            reordered.extend(chunks[i])
            reordered.append(spike)
        for chunk in chunks[len(spikeelements) :]:
            reordered.extend(chunk)
        return reordered

    ext_environments = _spikein(environments, num_cols=num_columns)
    num_columns += 1

    # load data
    recovered_contigs = get_recovered_contigs(fp_orb_basedir, settings, verbose)

    fig, axes = plt.subplots(
        int(np.ceil(len(ext_environments) / num_columns)) * 2,
        num_columns,
        figsize=(num_columns * 5, np.ceil((len(ext_environments) + 1) / num_columns) * 4),
        height_ratios=[1, 5] * (len(ext_environments) // num_columns),
        gridspec_kw={"hspace": 0.6, "wspace": 0.6},
    )

    skipped_d = False

    for i, environment in tqdm(enumerate(ext_environments), disable=not verbose, desc="Compute heatmap"):
        col = i % num_columns
        row = i // num_columns
        ax_dendro = axes[row * 2 + 0, col]
        ax_heat = axes[row * 2 + 1, col]

        # skip panel d for jaccard index
        if i == 3:
            pass
        else:
            if i > 3:
                i -= 1
            # panel labels
            ax_dendro.text(
                -0.33,
                1.05,
                chr(97 + i),
                transform=ax_dendro.transAxes,
                fontsize=16,
                fontweight="bold",
            )

        if environment == "colormap":
            ax_dendro.axis("off")
            ax_heat.set_position(
                [
                    ax_heat.get_position().x0,
                    ax_heat.get_position().y0,
                    ax_heat.get_position().width / 5,
                    ax_heat.get_position().height,
                ]
            )
            ax_heat.set_title("Jaccard distance")
            continue
        pd.set_option("future.no_silent_downcasting", True)
        if environment == "all six environments":
            features = pd.concat(
                [
                    pd.concat(
                        [
                            pd.Series(
                                index=list(v["gene_name"].unique()),
                                data=True,
                                name=k,
                            ).rename_axis("gene_name")
                            for k, v in recovered_contigs[env].items()
                        ],
                        axis=1,
                    )
                    .replace(np.nan, False)
                    .astype(bool)
                    for env in environments
                ]
            )
        else:
            features = (
                pd.concat(
                    [
                        pd.Series(
                            index=list(v["gene_name"].unique()),
                            data=True,
                            name=k,
                        ).rename_axis("gene_name")
                        for k, v in recovered_contigs[environment].items()
                    ],
                    axis=1,
                )
                .replace(np.nan, False)
                .astype(bool)
            )

        jaccard_distances = DistanceMatrix(
            squareform(pdist(features.T, metric="jaccard")),
            ids=[settings["labels"]["assemblers"].get(ass, ass) for ass in features.columns],
        ).to_data_frame()

        linkage_matrix = linkage(features.T, method="average", metric="jaccard")  # , optimal_ordering=False)
        col_dendro = dendrogram(linkage_matrix, no_plot=True)
        col_order = col_dendro["leaves"]
        sns.heatmap(
            jaccard_distances.iloc[col_order, col_order],
            ax=ax_heat,
            cbar=row == 0 and col == num_columns - 2,
            vmin=0,
            vmax=1,
            cbar_ax=axes[1, len(axes[0]) - 1],
        )

        dendrogram(linkage_matrix, ax=ax_dendro, color_threshold=0, no_labels=True)
        ax_dendro.axis("off")
        ax_dendro.set_title(settings["labels"]["environments"].get(environment, environment))
        ax_dendro.set_position(
            [
                ax_dendro.get_position().x0,
                ax_heat.get_position().y1,
                ax_dendro.get_position().width,
                ax_dendro.get_position().height,
            ]
        )

    return fig


def plot_rnaquast(fp_orb_basedir: str, fp_quast_basedir: str, settings, num_columns: int = 3, verbose=True):
    # which environments to plot and in which order
    environments = get_environments(fp_orb_basedir, settings)

    # load data
    quast = getdata_rnaquast(fp_orb_basedir, fp_quast_basedir, settings)
    data_recovery = getdata_recovery(fp_orb_basedir, settings, verbose=verbose)

    fig, axes = plt.subplots(
        int(np.ceil(len(environments) / num_columns)),
        num_columns * 2,
        figsize=(2 * num_columns * 4, np.ceil(len(environments) / num_columns) * 2.5),
        gridspec_kw={"wspace": 0.31, "hspace": 0.5},
    )

    for i, environment in tqdm(enumerate(environments), disable=not verbose, desc="Draw RNAquast panels"):
        order = list(data_recovery[data_recovery["environment"] == environment].sort_values(by="recovery_rank").index)

        # bad contigs
        ax_bad = axes[i // num_columns, (i % num_columns) * 2]
        sns.barplot(
            data=quast.loc[environment, :, "Misassemblies"],
            x="score",
            y="assembler",
            ax=ax_bad,
            order=order,
            color=settings["contig_classes"]["multi_mapped_contigs_multi_og"]["color"],
        )
        ax_bad.invert_xaxis()
        ax_bad.set_ylabel("")
        ax_bad.set_xlabel("number contigs")
        # ax_bad.set_title(settings['labels']['environments'].get(environment, environment), loc='right', horizontalalignment='center')
        ax_top_bad = ax_bad.twiny()
        ax_top_bad.xaxis.set_label_position("top")
        ax_top_bad.set_xticks([])
        ax_top_bad.set_xlabel(settings["labels"]["recovery_plot"]["bad"])
        ax_bad.xaxis.set_label_coords(1, -0.18)
        ax_bad.text(
            -0.5,
            1.05,
            chr(97 + i),
            transform=ax_bad.transAxes,
            fontsize=16,
            fontweight="bold",
        )

        ax_good = axes[i // num_columns, (i % num_columns) * 2 + 1]
        sns.barplot(
            data=quast.loc[environment, :, "95%-assembled isoforms"],
            x="score",
            y="assembler",
            ax=ax_good,
            order=order,
            color=settings["contig_classes"]["minimap2_single_recovered"]["color"],
        )
        ax_good.set_yticks([])
        ax_good.set_xlabel(settings["labels"]["recovery_plot"]["good"])
        ax_good.xaxis.set_label_position("top")
        ax_good.set_ylabel("")

        # concat right (=good) axis directly adjacent to left (=bad) axis
        ax_good.set_position(
            [
                ax_bad.get_position().x1,
                ax_good.get_position().y0,
                ax_good.get_position().width,
                ax_good.get_position().height,
            ]
        )

        # one legend for all panels
        if i + 1 == len(environments):
            ax_good.legend(
                handles=[
                    mpatches.Patch(
                        color=settings["contig_classes"]["multi_mapped_contigs_multi_og"]["color"],
                        label="Misassemblies",
                    ),
                    mpatches.Patch(
                        color=settings["contig_classes"]["minimap2_single_recovered"]["color"],
                        label="95%-assembled isoforms",
                    ),
                ],
                ncol=2,
                bbox_to_anchor=(-0.5, -0.30),
            )

        ax_bad.text(
            0.06,
            0.96,
            settings["labels"]["environments"].get(environment, environment),
            transform=ax_bad.transAxes,
            fontsize=12,
            fontweight="bold",
            va="top",
        )

    return fig


def plot_rnaquast_singlescore(
    fp_orb_basedir: str,
    fp_quast_basedir: str,
    settings,
    field: str = "Database coverage",
    num_columns: int = 3,
    verbose=True,
):
    # which environments to plot and in which order
    environments = get_environments(fp_orb_basedir, settings)

    # load data
    quast = getdata_rnaquast(fp_orb_basedir, fp_quast_basedir, settings)
    data_recovery = getdata_recovery(fp_orb_basedir, settings, verbose)

    fig, axes = plt.subplots(
        int(np.ceil(len(environments) / num_columns)),
        num_columns,
        figsize=(num_columns * 5, np.ceil(len(environments) / num_columns) * 2.5),
        gridspec_kw={"wspace": 0.7, "hspace": 0.6},
    )

    for i, environment in tqdm(enumerate(environments), disable=not verbose, desc="Draw RNAquast panels"):
        order = list(data_recovery[data_recovery["environment"] == environment].sort_values(by="recovery_rank").index)
        ax = axes[i // num_columns, i % num_columns]
        sns.barplot(data=quast.loc[environment, :, field], x="score", y="assembler", ax=ax, order=order)
        ax.set_ylabel("")
        ax.set_xlabel(field)
        ax.set_title(settings["labels"]["environments"].get(environment, environment))
        ax.text(
            -0.5,
            1.05,
            chr(97 + i),
            transform=ax.transAxes,
            fontsize=16,
            fontweight="bold",
        )

    return fig


def plot_detonate(fp_orb_basedir: str, fp_detonate_basedir: str, settings, num_columns: int = 3, verbose=True):
    # which environments to plot and in which order
    environments = get_environments(fp_orb_basedir, settings)

    # load data
    detonate = getdata_detonate(fp_orb_basedir, fp_detonate_basedir, settings)
    data_recovery = getdata_recovery(fp_orb_basedir, settings, verbose)

    fig, axes = plt.subplots(
        int(np.ceil(len(environments) / num_columns)),
        num_columns,
        figsize=(num_columns * 5, np.ceil(len(environments) / num_columns) * 2.5),
        gridspec_kw={"wspace": 0.7, "hspace": 0.6},
    )

    for i, environment in tqdm(enumerate(environments), disable=not verbose, desc="Draw RNAquast panels"):
        order = list(data_recovery[data_recovery["environment"] == environment].sort_values(by="recovery_rank").index)
        ax = axes[i // num_columns, i % num_columns]
        sns.barplot(data=detonate.loc[environment, :], x="score_{KC}", y="assembler", ax=ax, order=order)
        ax.set_ylabel("")
        ax.set_xlabel(r"$%s$" % ax.get_xlabel())
        # ax.set_xlabel(field)
        ax.set_title(settings["labels"]["environments"].get(environment, environment))
        ax.text(
            -0.5,
            1.05,
            chr(97 + i),
            transform=ax.transAxes,
            fontsize=16,
            fontweight="bold",
        )

    return fig


def plot_block_correlation(
    fp_orb_basedir,
    fp_marbel_basedir,
    sequence_file,
    settings,
    field,
    verbose=True,
    test="Mann-Whitney-gt",
    data_block_recovery=None,
    num_columns: int = 3,
    logscale=False,
):
    # field: block_length, read_mean_count, mean_identity, Coverage
    environments = get_environments(fp_orb_basedir, settings)

    fig, axes = plt.subplots(
        int(np.ceil(len(environments) / num_columns)),
        num_columns,
        figsize=(int(np.ceil(len(environments) / num_columns)) * 4, num_columns * 6),
        gridspec_kw={"wspace": 0.6, "hspace": 0.4},
    )

    axes = axes.flatten()
    hue_order = ["recovered", "missed"]
    if data_block_recovery is None:
        data_block_recovery = getdata_block_recovery(
            fp_orb_basedir, fp_marbel_basedir, sequence_file, settings, verbose=verbose
        )
    assemblers = settings["labels"]["assemblers"]
    # 'Mann-Whitney', 't-test_ind', 'Wilcoxon'

    block_colors = {"recovered": "lightgreen", "missed": "#dddddd"}
    for i, environment in tqdm(enumerate(environments), disable=not verbose, desc=f"Drawing panels for {field} plot"):
        ax = axes[i]
        plotdata = data_block_recovery[environment]

        # assembler
        ordered_assembler = sorted(list(plotdata["assembler"].unique()), key=lambda x: x.lower())
        if field == "block_length":
            plotdata = plotdata[plotdata[field] > settings["read_length"]]

        sns.boxplot(
            data=plotdata,
            orient="h",
            y="assembler",
            x=field,
            hue="category",
            ax=ax,
            palette=block_colors,
            showfliers=False,
            hue_order=hue_order,
            order=ordered_assembler,
        )
        if logscale:
            ax.set_xscale("log")
        if environment != settings["environment_for_legend_display"]:
            ax.legend().remove()
        else:
            ax.legend(bbox_to_anchor=(1.6, -0.2), title="blocks", ncols=2)
        # fix_ass_labels(ax, dimension='Y')
        ax.set_ylabel("")
        ax.set_xlabel(
            {
                "read_mean_count": "coverage: mean read number",
                "block_length": "block length (bps)",
                "mean_identity": "sequence similarity",
            }.get(field, field)
        )
        ax.set_title(settings["labels"]["environments"].get(environment, environment))
        ax.text(
            -0.43,
            1.05,
            chr(97 + i),
            transform=ax.transAxes,
            fontsize=16,
            fontweight="bold",
        )

        ns = plotdata.groupby(["assembler", "category"]).size()
        ax_ns = ax.twinx()
        ax_ns.set_ylim(ax.get_ylim())
        ns_labels = []
        ns_colors = []
        for ass_label in ordered_assembler:
            for cat in hue_order:
                ns_labels.append(f"n={ns.loc[ass_label, cat]:,}")
                ns_colors.append({"lightgreen": "darkgreen", "#dddddd": "darkgray"}.get(block_colors[cat]))
        ax_ns.set_yticks([x + offset for x in ax.get_yticks() for offset in [-0.2, 0.2]], labels=ns_labels, fontsize=8)
        for color, label in zip(ns_colors, ax_ns.get_yticklabels()):
            label.set_color(color)

        if test is not None:
            annotator = Annotator(
                ax,
                [((assembler, "recovered"), (assembler, "missed")) for assembler in ordered_assembler],
                data=plotdata,
                x=field,
                y="assembler",
                hue="category",
                orient="h",
                hue_order=hue_order,
                order=ordered_assembler,
                verbose=0,
            )
            annotator.configure(
                test=test,
                text_format="star",
                loc="inside",
                comparisons_correction="fdr_bh",
                correction_format="default",
            )
            testres = annotator.apply_and_annotate()
        else:
            testres = pd.DataFrame(columns=[1])

        # gather mean values with expections of those assembler that yield no sig. difference
        except_assembler = [it.structs[0]["group"][0] for it in testres[1] if it.data._corrected_significance is False]
        if len(except_assembler) == len(ordered_assembler):
            except_assembler = []
        mean_field_diff = (
            data_block_recovery[environment]
            .sort_values(["assembler"])
            .set_index(["assembler"])
            .loc[[a for a in ordered_assembler if a not in except_assembler], :]
            .groupby("category")[field]
            .median()
        )

        ax_top = ax.twiny()
        ax_top.set_xticks([])
        ax_top.set_xlabel(
            "⌀ recovered ~ %.2f, ⌀ missed ~ %.2f" % (mean_field_diff.loc["recovered"], mean_field_diff["missed"])
        )

        if (field == "block_length") and (logscale is False):
            ax.set_xticks(
                [0, 100, 200, 300, 400, 500, 1000, 2000, 3000],
                labels=["0", "", "", "", "", "500", "1000", "2000", "3000"],
            )

    return fig


def plot_nruns(fp_orb_basedir, settings, num_columns: int = 3, verbose=True, contigs=None):
    sns.set_style("ticks")
    environments = get_environments(fp_orb_basedir, settings)

    fig, axes = plt.subplots(
        int(np.ceil(len(environments) / num_columns)),
        num_columns * 2,
        figsize=(num_columns * 6, int(np.ceil(len(environments) / num_columns)) * 1.5),
        gridspec_kw={"wspace": 0.51, "hspace": 1.3},
    )

    axes = axes.flatten()
    if contigs is None:
        contigs = filter_for_assembler_with_ns(fp_orb_basedir, settings)

    for i, environment in tqdm(enumerate(environments), disable=not verbose, desc=f"Drawing panels for Ns plot"):
        n_run_data = contigs[environment]
        plot_data = n_run_data.groupby(["assembler", "class", "max_n_longer_read_length"]).size().reset_index()
        ax_bad = axes[i * 2]
        sns.barplot(
            data=plot_data[plot_data["class"] == "missed"],
            y="assembler",
            x=0,
            hue="max_n_longer_read_length",
            ax=ax_bad,
        )
        ax_bad.set_xlim((0, plot_data[0].max() * 1.1))
        ax_bad.invert_xaxis()
        ax_bad.set_ylabel("")
        # ax_bad.set_title(settings['labels']['environments'].get(environment, environment), loc='right', horizontalalignment='center')

        # add panel label
        ax_bad.text(
            -0.5,
            1.05,
            chr(97 + i),
            transform=ax_bad.transAxes,
            fontsize=16,
            fontweight="bold",
        )

        ax_bad.set_xlabel("number contigs")
        ax_bad.xaxis.set_label_coords(1, -0.4)
        ax_top_bad = ax_bad.twiny()
        ax_top_bad.set_xticks([])
        ax_top_bad.set_xlabel(settings["labels"]["recovery_plot"]["bad"])

        ax_good = axes[i * 2 + 1]
        sns.barplot(
            data=plot_data[plot_data["class"] != "missed"],
            y="assembler",
            x=0,
            hue="max_n_longer_read_length",
            ax=ax_good,
        )
        ax_good.tick_params(left=False)
        ax_good.set_yticks([])
        ax_good.set_ylabel("")
        ax_good.set_xlabel(settings["labels"]["recovery_plot"]["good"])
        ax_good.xaxis.set_label_position("top")
        ax_good.set_xlim((0, plot_data[0].max() * 1.1))

        ax_good.set_position(
            [
                ax_bad.get_position().x1,
                ax_good.get_position().y0,
                ax_good.get_position().width,
                ax_good.get_position().height,
            ]
        )

        ax_good.text(
            0.95,
            0.90,
            settings["labels"]["environments"].get(environment, environment),
            transform=ax_good.transAxes,
            fontsize=12,
            fontweight="bold",
            va="top",
            ha="right",
        )

        # create one joined legend
        if environment == settings["environment_for_legend_display"]:
            handles, lables = ax_good.get_legend_handles_labels()
            ax_good.legend(
                handles,
                [l.replace("read length", "read length (%ibp)" % settings["read_length"]) for l in lables],
                bbox_to_anchor=(2.8, -0.9),
                ncol=2,
            )
        else:
            ax_good.legend().remove()
        ax_bad.legend().remove()
    return fig


def plot_hit_genomes_mapping(mapping, ax=None):
    assert mapping["Query sequence name"].unique().shape[0] == 1
    CONTIG_GAP = 100000
    hitYpos = -2

    plotinfo = mapping.copy()

    # one line per organism
    for y, (taxon, g) in enumerate(plotinfo.groupby("assemblyAccession")):
        plotinfo.loc[g.index, "ypos"] = y
        # one segment per genome assembly contig
        x_pos = 0
        for i, (contig, g2) in enumerate(g.sort_values(by="assembly_length").groupby("assembly_accession")):
            for field in [
                "gene_genomic_start",
                "gene_genomic_stop",
                "block_genomic_start",
                "block_genomic_stop",
                "hit_start",
                "hit_stop",
            ]:
                plotinfo.loc[g2.index, "xpos_%s" % field] = x_pos + plotinfo.loc[g2.index, field]
            plotinfo.loc[g2.index, "xpos_assembly_contig_start"] = x_pos
            plotinfo.loc[g2.index, "xpos_assembly_contig_stop"] = x_pos + g2["assembly_length"].iloc[0]
            x_pos += g2["assembly_length"].iloc[0] + CONTIG_GAP

    if ax is None:
        fig, axes = plt.subplots(1, 1, figsize=(10, 3))
        ax = axes
    else:
        fig = None

    ylabels = []
    total_width = plotinfo["xpos_assembly_contig_stop"].max()
    colors = sns.color_palette(n_colors=mapping.shape[0])
    color_index = 0
    for y, (taxon, g) in enumerate(plotinfo.groupby("assemblyAccession")):
        # plot horizontal lines for the contigs of the genome assembly
        ylabels.append((y, g["organism"].iloc[0] + "\n" + taxon))
        for _, asscontig in g.groupby("assembly_accession").head(1).iterrows():
            ax.plot(
                [asscontig["xpos_assembly_contig_start"], asscontig["xpos_assembly_contig_stop"]],
                [asscontig["ypos"], asscontig["ypos"]],
                color="black",
                linewidth=1,
                zorder=-1,
            )

        # plot boxes (will be reduced to vertical lines in most cases) for genes in genomes
        for gene, g_gene in g.sort_values(by="xpos_gene_genomic_start").groupby("gene_name").head(1).iterrows():
            width = max(g_gene["xpos_gene_genomic_stop"] - g_gene["xpos_gene_genomic_start"] + 1, total_width / 500)
            ax.add_patch(mpatches.Rectangle((g_gene["xpos_gene_genomic_start"], g_gene["ypos"] - 0.2), width, 0.2 * 2))

        for _, hit in g.sort_values(by="Query start coordinate").iterrows():
            hit_start, hit_stop = hit["Query start coordinate"], hit["Query end coordinate"]
            if True:
                hit_start = hit_start / hit["Query sequence length"] * total_width
                hit_stop = hit_stop / hit["Query sequence length"] * total_width
            ax.add_patch(
                mpatches.Polygon(
                    [
                        [hit_start, hitYpos],
                        [hit_stop, hitYpos],
                        [hit["xpos_hit_start"], hit["ypos"]],
                        [hit["xpos_hit_stop"], hit["ypos"]],
                    ],
                    closed=True,
                    facecolor=colors[color_index],
                    edgecolor="k",
                    linewidth=2,
                    alpha=0.6,
                )
            )
            color_index += 1

    # plot horizontal line for hit contig
    ax.plot([0, total_width], [hitYpos, hitYpos], color="blue", linewidth=1, zorder=-1)
    ylabels.append((hitYpos, plotinfo.iloc[0, :]["Query sequence name"]))

    ax.set_yticks([pos for pos, _ in ylabels], [label for _, label in ylabels])
    ax_top = ax.twiny()

    # give reasonable xticks for hit (but in genomic coordinates)
    hit_xticks = ticker.MaxNLocator(nbins=7).tick_values(0, plotinfo.iloc[0, :]["Query sequence length"])
    ax.set_xticks([x / hit["Query sequence length"] * total_width for x in hit_xticks], list(map(int, hit_xticks)))
    ax.set_xlabel("contig position (bp)")

    genome_xticks = ticker.MaxNLocator(nbins=7).tick_values(0, total_width)
    ax_top.set_xlim(ax.get_xlim())
    ax_top.set_xticks(genome_xticks, ["%.1f" % (x / 1000000) for x in genome_xticks])
    ax_top.set_xlabel("genomic position (Mbp)")

    if fig is None:
        return ax
    else:
        return fig
