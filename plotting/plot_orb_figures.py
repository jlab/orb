from plot_include_orb import *
import sys
import os


def save_plot(fig, plot_name, out_dir, file_ending):
    if file_ending == "png":
        fig.savefig(f"{out_dir}/{plot_name}.png", dpi=300, bbox_inches="tight")
    else:
        fig.savefig(f"{out_dir}/{plot_name}.svg")


fp_orb = "/vol/jlab/tlin/all_project/nf_results/orb/recovered_fix"
fp_marbel_basedir = "/vol/jlab/tlin/marbel_benchmarking_integration/benchmarking_sets_all_sparse_fixed_libsize"
sequence_df_file = "/homes/tlin/Projects/marbel/src/marbel/data/deduplicated_pangenome_EDGAR_Microbiome_JLAB2.fas.bgz.bio_index"

if len(sys.argv) < 2:
    print("Usage: plot_orb_figures.py <fp_orb_basedir> <fp_marbel_basedir> <marbel_sequence_file> <settings> <file_ending_svg_or_png> <outdir_name> <caviar_log_files>")
    sys.exit(1)


fp_orb = sys.argv[1]

fp_marbel_basedir = sys.argv[2]

marbel_sequence_file = sys.argv[3]

style_yaml = sys.argv[4]

file_ending = sys.argv[5]

out_dir = sys.argv[6]

caviar_log = sys.argv[7]

if not os.path.exists(out_dir):
    os.mkdir(out_dir)

with open(style_yaml, "r") as f:
    settings = yaml.safe_load(f)

fig = plot_recovery(fp_orb, settings)

save_plot(fig, "recovery_plot", out_dir, file_ending)

fig = plot_recovery(fp_orb, settings, detail_view=True)

save_plot(fig, "recovery_plot_detail", out_dir, file_ending)

fig = plot_gene_recovery(fp_orb, settings)

save_plot(fig, "gene_recovery", out_dir, file_ending)

if caviar_log != "":

    fig = plotTimeMemory(caviar_log, settings)

    save_plot(fig, "caviar_log", out_dir, file_ending)

fig = plot_DEgenes(fp_orb, settings)

save_plot(fig, "de_genes", out_dir, file_ending)

fig = plot_DEgenes(fp_orb, settings, forOrthogroups=True)

save_plot(fig, "de_orthogroups", out_dir, file_ending)

fig = plot_DEvennOrtho(fp_orb, fp_marbel_basedir, settings)

save_plot(fig, "de_venn_orthogroups", out_dir, file_ending)

fig = plot_heatmap(fp_orb, settings)

save_plot(fig, "dissimilarity_heatmap", out_dir, file_ending)

fig = plot_nruns(fp_orb, settings)

save_plot(fig, "n_runs", out_dir, file_ending)

fig = plot_block_correlation(fp_orb, fp_marbel_basedir, marbel_sequence_file, settings, "block_length")

save_plot(fig, "block_recovery_correlation_block_length", out_dir, file_ending)

fig = plot_block_correlation(fp_orb, fp_marbel_basedir, marbel_sequence_file, settings, "Coverage")

save_plot(fig, "block_recovery_correlation_coverage", out_dir, file_ending)

fig = plot_block_correlation(fp_orb, fp_marbel_basedir, marbel_sequence_file, settings, "mean_identity")

save_plot(fig, "block_recovery_correlation_mean_identity", out_dir, file_ending)
