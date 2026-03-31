[![GitHub Actions Linting Status](https://github.com/jlab/refbasedassemblereval/workflows/nf-core%20linting/badge.svg)](https://github.com/jlab/refbasedassemblereval/actions?query=workflow%3A%22nf-core+linting%22)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/jlab/refbasedassemblereval)

# ORB (Omics reference based benchmarking)

## Introduction

![ORB logo](./resources/logos/orb_logo.png)

**jlab/orb** is a bioinformatics pipeline that calculates performance evaluation scores for assembled sets of contigs. Using Marbel, a researcher is enabled to create an in silico dataset resembling the characteristics of the target environment and, using ORB, test which assembler to use for the analysis of their sample. The pipeline leverages minimap2, Bowtie2, Salmon, DESeq2, edgeR, Calour and custom scripts for the score calculation and includes orthologous groups and DE benchmarking.

Recommended usage:

- Create _in silico_ datasets using [Marbel](https://anaconda.org/bioconda/marbel), name the output dirs: NAME_microbiome
- Assemble the datasets
- Run the pipeline

Fill the required parameters per run and adjust the example config files in `example/tool.config`, `example/dataset.config` and `example/resources.config`.

For each dataset create a dataset.config with the data set name. For outdir parameter chose the same parent dir, if you want to visualise the datasets together, e.g., PATH/group/dataset1, PATH/group/dataset2.

Run each dataset with:

```bash
nextflow run . \
   -profile apptainer \
   -c example/dataset.config \
   -c example/resources.config \
   -c example/tool.config
```

Afterwards there is multiple methods which can visualise the runs. A script which creates all orb files is provided:
`plotting/plot_orb_figures.py

For the dependencies you can use:

```
conda env create -f plotting/orb_plots.yaml

conda activate orb_plots
```

The script can be startet with:

python plotting/plot_orb_figures.py <fp_orb_basedir> <fp_marbel_basedir> <marbel_sequence_file> <settings> <file_ending_svg_or_png> <outdir_name> <caviar_log_files>

`fp_orb_basedir`: Path of the results for orb.

`fp_marbel_basedir`: Path to the _in silico_ datasets. Datasets folders require \_microbiome suffix.

`marbel_sequence_file`: Path to the bio index file of the [Marbel repository](https://github.com/jlab/marbel): `src/marbel/data/deduplicated_pangenome_EDGAR_Microbiome_JLAB2.fas.bgz.bio_index`

`settings`: Style yaml, for example style settings please see plotting/style.yaml

`file_ending_svg_or_png`: file ending, supported: svg or png

`outdir_name`: name of the directory where the plots should be saved

`caviar_log_files`: directory of the log files, if assembled with [Caviar](https://github.com/jlab/caviar), can be left blank with: ""

This script will take some time for the first run. If you recreate environments, remove cache folder: data_cache.

Individual plots may also be imported in a notebook from plot_include_orb and be created there.

## Credits

jlab/refbasedassemblereval was originally written by Timo Wentong Lin.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE). Some code from the nf-core community was modified and is posted alongside own code in `modules/local`.

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
