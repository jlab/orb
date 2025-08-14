[![GitHub Actions CI Status](https://github.com/jlab/refbasedassemblereval/workflows/nf-core%20CI/badge.svg)](https://github.com/jlab/refbasedassemblereval/actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting Status](https://github.com/jlab/refbasedassemblereval/workflows/nf-core%20linting/badge.svg)](https://github.com/jlab/refbasedassemblereval/actions?query=workflow%3A%22nf-core+linting%22)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/jlab/refbasedassemblereval)

# ORB (Omics reference based benchmarking)

## Introduction

**jlab/orb** is a bioinformatics pipeline that calculates performance evaluation scores for assembled sets of contigs. Using Marbel, a researcher is enabled to create an in silico dataset resembling the characteristics of the target environment and, using ORB, test which assembler to use for the analysis of their sample. The pipeline leverages minimap2, Bowtie2, Salmon, DESeq2, edgeR, Calour and custom scripts for the score calculation and includes orthologous groups and DE benchmarking.

Recommended usage:

* Create an *in silico* dataset using [Marbel](https://anaconda.org/bioconda/marbel)
* Assemble the datasets with a tool to benchmark
* Run the pipeline

Fill the required parameters per run and adjust the example config files in `example/dataset.config` and `example/resources.config`

```bash
nextflow run . \
   -profile podman \
   -c example/dataset.config \
   -c example/resources.config
```

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
