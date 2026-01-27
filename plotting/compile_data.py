from glob import glob
from os.path import join, basename, relpath, sep
import yaml
import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy.spatial.distance import pdist, squareform


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

def getdata_recovery(fp_orb_basedir:str, settings, verbose=True) -> pd.DataFrame:
    """Collects data from Orb for contig recovery analysis.
    
    Parameters
    ----------
    fp_orb_basedir : str
        The filepath to Orb's base data dir
    settings : yaml-dict
        General settings for plotting Orb graphs. Here, we need
        a) the environment to skip
    verbose : boolean
        Report progress on sys.stderr
    Returns
    -------
    A pandas.DataFrame that holds performance parameters for assemblers in the different environments.
    """
    environments = get_environments(fp_orb_basedir, settings)

    data = []
    for environment in tqdm(environments, disable=not verbose, desc='Compiling data for contig recovery plot'):
        # load orb data
        orb = pd.read_csv(join(fp_orb_basedir, environment, 'mergedresults', '%s_all_scores.tsv' % environment), sep="\t", index_col=0)
        # rename metrics to names used in publication
        orb = orb.rename(index={k: c['label'] for k, c in settings['contig_classes'].items()})
        # add information about missed blocks
        trueBlocks = pd.read_csv(join(fp_orb_basedir, environment, 'mergedataframessummaries', '%s_all_contigs.tsv' % environment), sep="\t", index_col=0).loc['count', ['Blocks']].sum()
        orb.loc['missed blocks', :] = trueBlocks - orb.loc[[c['label'] for _, c in settings['contig_classes'].items() if c['class'] == 'good'], :].sum()
        # select metrics for "recovery analysis"
        orb = orb.loc[[c['label'] for _, c in settings['contig_classes'].items()], :]
        # sort assembler by amount of good contigs
        orb = orb[orb.loc[[c['label'] for _, c in settings['contig_classes'].items() if c['class'] == 'good']].sum().sort_values(ascending=True).index]
        # use pretty label for assembler
        orb = orb.rename(columns=settings['labels']['assemblers'])
        # transform axis: rows=assembler, cols=metrics+metadata
        orb = orb.T
        orb['environment'] = environment
        orb['recovery_rank'] = list(reversed(range(1, orb.shape[0] + 1)))
        data.append(orb)
    return pd.concat(data)

def get_recovered_contigs(fp_orb_basedir:str, settings, verbose=True):
    mapping_col_names = [
                "Query sequence name",
                "Query sequence length",
                "Query start coordinate",
                "Query end coordinate",
                "same strand",
                "block_id", #"Target sequence name",
                "Target sequence length",
                "Target start coordinate on the original strand",
                "Target end coordinate on the original strand",
                "Number of matching bases in the mapping",
                "Number bases, including gaps, in the mapping",
                "Mapping quality", "1", "2", "3", "4", "5", "6"]

    # which environments to plot and in which order
    environments = get_environments(fp_orb_basedir, settings)

    recovered_contigs = dict()
    for environment in tqdm(environments, disable=not verbose, desc='Compiling data for gene recovery plot'):
        recovered_contigs[environment] = dict()
        for fp_category in glob(join(fp_orb_basedir, environment, 'categorizecontigs', '*_contigs_categorised.tsv')):
            assembler = basename(fp_category).split('_contigs_categorised.tsv')[0]
            
            contig_classes = pd.read_csv(fp_category, sep="\t", index_col=0)
            contig_classes.index = list(map(str, contig_classes.index))
            print("test")
            contig_mappings = pd.read_csv(
                join(fp_orb_basedir, environment, 'minimap2', '%s_mapping.tsv' % assembler), sep="\t", header=None, names=mapping_col_names, dtype={'Query sequence name': str}
                    ).sort_values(by=['Mapping quality', 'Number of matching bases in the mapping'], ascending=[False, False]  # sort hits by mapping quality AND number of involved matching bases
                        ).groupby(['Query sequence name']).head(1  # for every contig: pick only best hit
                            ).groupby('block_id').head(1).merge(  # for every block: pick only best hit
                                contig_classes, left_on=['Query sequence name'], right_index=True, how='left')  # merge Timo's contig categorization
            # add a binary classification: good and missed ...
            contig_mappings['class'] = contig_mappings['category'].apply(lambda x: 'good' if x in [col for col, c in settings['contig_classes'].items() if c['class'] == 'good'] else 'missed')
            # ... and keep only those in class "good"
            contig_mappings = contig_mappings[contig_mappings['class'] ==  'good']
            
            # derive gene name from block_id
            contig_mappings['gene_name'] = contig_mappings['block_id'].apply(lambda x: x.split('_block')[0])        
            
            recovered_contigs[environment][assembler] = contig_mappings

    return recovered_contigs

def getdata_gene_recovery(fp_orb_basedir:str, settings, full_gene_table=False, verbose=True) -> pd.DataFrame:
    recovered_contigs = get_recovered_contigs(fp_orb_basedir, settings, verbose=True)
    # recovered_contigs = dict()
    # for environment in tqdm(environments, disable=not verbose, desc='Compiling data for gene recovery plot'):
    #     recovered_contigs[environment] = dict()
    #     for fp_category in glob(join(fp_orb_basedir, environment, 'categorizecontigs', '*_contigs_categorised.tsv')):
    #         assembler = basename(fp_category).split('_contigs_categorised.tsv')[0]
            
    #         contig_classes = pd.read_csv(fp_category, sep="\t", index_col=0)
    #         contig_classes.index = list(map(str, contig_classes.index))

    #         contig_mappings = pd.read_csv(
    #             join(fp_orb_basedir, environment, 'minimap2', '%s_mapping.tsv' % assembler), sep="\t", header=None, names=mapping_col_names, dtype={'Query sequence name': str}
    #                 ).sort_values(by=['Mapping quality', 'Number of matching bases in the mapping'], ascending=[False, False]  # sort hits by mapping quality AND number of involved matching bases
    #                     ).groupby(['Query sequence name']).head(1  # for every contig: pick only best hit
    #                         ).groupby('block_id').head(1).merge(  # for every block: pick only best hit
    #                             contig_classes, left_on=['Query sequence name'], right_index=True, how='left')  # merge Timo's contig categorization
    #         # add a binary classification: good and missed ...
    #         contig_mappings['class'] = contig_mappings['category'].apply(lambda x: 'good' if x in [col for col, c in settings['contig_classes'].items() if c['class'] == 'good'] else 'missed')
    #         # ... and keep only those in class "good"
    #         contig_mappings = contig_mappings[contig_mappings['class'] ==  'good']
            
    #         # derive gene name from block_id
    #         contig_mappings['gene_name'] = contig_mappings['block_id'].apply(lambda x: x.split('_block')[0])        
            
    #         recovered_contigs[environment][assembler] = contig_mappings

    # which environments to plot and in which order
    environments = get_environments(fp_orb_basedir, settings)

    recovered_genes = []
    for environment in tqdm(environments, disable=not verbose, desc='Compute sets of unique/shared/core genes for gene recovery plot'):
        #pd.set_option('future.no_silent_downcasting', True)
        features = pd.concat(
            [pd.Series(
                index=list(v['gene_name'].unique()),
                data=True,
                name=k,
             ).rename_axis('gene_name')
             for k, v in recovered_contigs[environment].items()
            ], axis=1).replace(np.nan, False).astype(bool)
        if full_gene_table:
            recovered_genes.append((environment,  features.rename(index=settings['labels']['assemblers'])))
            continue

        del features['idba_mt']  # as this is too close to idba_tran
        # number genes found by only one assembler
        upset = features[features.sum(axis=1) == 1].unstack().unstack().sum(axis=1)
        # number genes found by all assemblers
        upset.loc['core'] = features[features.sum(axis=1) == len(features.columns)].shape[0]
        # number genes found by two or more, but not all assemblers
        upset.loc['shared'] = features[(features.sum(axis=1) < len(features.columns)) & (features.sum(axis=1) > 1)].shape[0]
        # total number genes
        upset.loc['total'] = features.shape[0]
        upset.name = environment
        recovered_genes.append(upset)
    if full_gene_table:
        return {environment: feature_table for (environment, feature_table) in recovered_genes}
    recovered_genes = pd.concat(recovered_genes, axis=1)
    recovered_genes.index.name = 'assembler'

    # pretty assembler names
    recovered_genes = recovered_genes.rename(index=settings['labels']['assemblers'])

    return recovered_genes

def getdata_runtime_memory(fp_caviar_basedir:str, settings, verbose=True):
    fields_to_collect = ['Maximum resident set size (kbytes): ', 'User time (seconds): ', 'System time (seconds): ']
    data = []
    for fp_time in tqdm(glob('%s/**/*_time_log.txt' % fp_caviar_basedir, recursive=True), disable=not verbose, desc='Compiling data for runtime/memory footprint plot'):
        # derive environment and assembler from filenames
        environment = relpath(fp_time, fp_caviar_basedir).split(sep)[0]
        assembler = (basename(fp_time)[len(environment)+1:-1*len('_time_log.txt')])

        # parse /usr/bin/time verbose output
        with open(fp_time, 'r') as f:
            for line in f.readlines():
                for field in fields_to_collect:
                    if field in line:
                        value = float(line.split(field)[-1].strip())
                        data.append((environment, assembler, field, value))
    data = pd.DataFrame(data, columns=['environment', 'assembler', 'type', 'value'])

    # combine user and system time
    cmb_time = (data[data['type'] == 'User time (seconds): '].set_index(['environment', 'assembler'])['value'] + data[data['type'] == 'System time (seconds): '].set_index(['environment', 'assembler'])['value']).reset_index()
    cmb_time['type'] = 'CPU time (seconds)'

    timemem = pd.concat([data, cmb_time])
    timemem = timemem[timemem['assembler'] != 'convert_to_fasta']  # this is a pre-processing for some of the assembler as they not always except fastQ.gz
    
    # combine runtime for oases, i.e. velveth + velvetg + oases
    assembler = 'oases'
    oases_combined = []
    for field in fields_to_collect + ['CPU time (seconds)']:
        fct_combine = np.sum
        if field == 'Maximum resident set size (kbytes): ':
            fct_combine = np.max
        combined = timemem.set_index(['environment', 'assembler', 'type']).loc[:, ['velveth', 'velvetg', assembler], field, :].reset_index().groupby('environment')['value'].apply(fct_combine).to_frame().reset_index()
        combined['type'] = field
        combined['assembler'] = assembler
        oases_combined.append(combined)
    
    timemem = pd.concat([
        timemem[~timemem['assembler'].isin(['velveth', 'velvetg', assembler])], # result without part of oases
        pd.concat(oases_combined)]  # combined oases results
                       )
    # pretty assembler names
    timemem['assembler'] = timemem['assembler'].apply(lambda x: settings['labels']['assemblers'].get(x, x))

    # pretty environment names
    timemem['environment'] = timemem['environment'].apply(lambda x: settings['labels']['environments'].get(x, x))

    return timemem

def getdata_DEgenes(fp_orb_basedir:str, settings, verbose=True):
    recovered_contigs = get_recovered_contigs(fp_orb_basedir, settings, verbose=True)

    # which environments to plot and in which order
    environments = get_environments(fp_orb_basedir, settings)

    confusion = []
    for environment in tqdm(environments, disable=not verbose, desc='Compiling data for DE gene plot'):
        fp_truth = join(fp_orb_basedir, environment, 'refdeseq2', '%s_DESeq2_full_table.tsv' % environment)
        truth = pd.read_csv(fp_truth, sep="\t", index_col=0)
        truth.index = list(map(str, truth.index))

        # flag all genes which are DE
        truth = truth['padj'].apply(lambda x: x < 0.05).rename('truth')

        for fp_category in sorted(glob(join(fp_orb_basedir, environment, 'deseq2', '*_DESeq2_full_table.tsv'))):
            assembler = basename(fp_category).split('_DESeq2_full_table.tsv')[0]
            prediction = pd.read_csv(fp_category, sep="\t", index_col=0)
            prediction.index = list(map(str, prediction.index))
            
            # only keep contigs that have been identified as being statistically significant
            # flag all remaining contigs as DE prediction
            prediction = prediction['padj'].apply(lambda x: x < 0.05).rename('prediction')

            #pd.set_option('future.no_silent_downcasting', True)
            # combine assembled contigs with DE truth and prediction and count occurences
            conv = recovered_contigs[environment][assembler].merge(
                truth, left_on='gene_name', right_index=True, how='outer').merge(
                    prediction, left_on='Query sequence name', right_index=True, how='outer').groupby('gene_name').head(1).fillna(False).groupby(
                        ['truth', 'prediction']).size()

            # re-structure convolution data
            conv = conv.rename('num_genes').to_frame()
            # give more speaking names
            conv['class'] = list(map(lambda row: {(False, False): 'True Negative',
                                                  (False, True): 'False Positive',
                                                  (True, False): 'False Negative',
                                                  (True, True): 'True Positive'}.get(row[0], row[0]), conv.iterrows()))
            # add in environment + assembler info
            conv['environment'] = settings['labels']['assemblers'].get(environment, environment)
            conv['assembler'] = settings['labels']['assemblers'].get(assembler, assembler)
            conv = conv.reset_index()
            confusion.append(conv)
    
    # re-structure into one dataframe
    confusion = pd.concat(confusion, axis=0).set_index(['environment', 'assembler', 'class']).sort_index()
    confusion['rank'] = np.nan

    # add rank information to return dataframe
    for environment in confusion.index.levels[0]:
        order = list(reversed(
            pd.pivot_table(data=confusion.loc[environment, :], index='assembler', columns='class', values='num_genes', aggfunc="sum").sort_values(
                by=       ['True Positive', 'False Positive', 'False Negative'],
                ascending=[False,           True,             True]).index))
        confusion.loc[confusion.loc[environment, order, :].index, 'rank'] = [
            rank
            for rank in list(reversed(range(1, len(order) + 1)))
            for i in range(confusion.reset_index()['class'].unique().shape[0])]

    return confusion, recovered_contigs


def getdata_DEorthogroups(fp_orb_basedir:str, fp_marbel_basedir:str, fp_ogtruth_basedir:str, settings, verbose=True):
    recovered_contigs = get_recovered_contigs(fp_orb_basedir, settings, verbose=True)

    # which environments to plot and in which order
    environments = get_environments(fp_orb_basedir, settings)

    confusion = []
    for environment in tqdm(environments, disable=not verbose, desc='Compiling data for DE orthogroup plot'):
        # obtain orthogroup information about genes
        genes = pd.read_csv(join(fp_marbel_basedir, '%s_microbiome' % environment, "summary", "gene_summary.csv"), sep=",").set_index('gene_name')
        fp_truth = join(fp_ogtruth_basedir, environment, '%s_DESeq2_full_table.tsv' % environment)
        truth = pd.read_csv(fp_truth, sep="\t", index_col=0)
        truth.index = list(map(str, truth.index))

        # flag all genes which are DE
        truth = truth['padj'].apply(lambda x: x < 0.05).rename('truth')

        for fp_assembler in sorted(glob(join(fp_ogtruth_basedir, environment, '%s_*_DESeq2_full_table.tsv' % environment))):
            assembler = basename(fp_assembler).split('_DESeq2_full_table.tsv')[0].split('%s_' % environment)[-1]
            prediction = pd.read_csv(fp_assembler, sep="\t", index_col=0)
            prediction.index = list(map(str, prediction.index))

            prediction = prediction['padj'].apply(lambda x: x < 0.05).rename('prediction')

            # pd.set_option('future.no_silent_downcasting', True)
            conv = recovered_contigs[environment][assembler].merge(
                genes[['orthogroup']], left_on='gene_name', right_index=True, how='left').merge(
                    truth, left_on='orthogroup', right_index=True, how='outer').merge(
                        prediction, left_on='orthogroup', right_index=True, how='outer').groupby(
                            'orthogroup').head(1).fillna(False).groupby(
                                ['truth', 'prediction']).size()

            # re-structure convolution data
            conv = conv.rename('num_genes').to_frame()
            # give more speaking names
            conv['class'] = list(map(lambda row: {(False, False): 'True Negative',
                                                  (False, True): 'False Positive',
                                                  (True, False): 'False Negative',
                                                  (True, True): 'True Positive'}.get(row[0], row[0]), conv.iterrows()))
            # add in environment + assembler info
            conv['environment'] = settings['labels']['assemblers'].get(environment, environment)
            conv['assembler'] = settings['labels']['assemblers'].get(assembler, assembler)
            conv = conv.reset_index()
            confusion.append(conv)

    # re-structure into one dataframe
    confusion = pd.concat(confusion, axis=0).set_index(['environment', 'assembler', 'class']).sort_index()
    confusion['rank'] = np.nan

    # add rank information to return dataframe
    for environment in confusion.index.levels[0]:
        order = list(reversed(
            pd.pivot_table(data=confusion.loc[environment, :], index='assembler', columns='class', values='num_genes', aggfunc="sum").sort_values(
                by=       ['True Positive', 'False Positive', 'False Negative'],
                ascending=[False,           True,             True]).index))
        confusion.loc[confusion.loc[environment, order, :].index, 'rank'] = [
            rank
            for rank in list(reversed(range(1, len(order) + 1)))
            for i in range(confusion.reset_index()['class'].unique().shape[0])]

    return confusion, recovered_contigs


def getdata_DEvennOrtho(fp_orb_basedir:str, fp_ogtruth_basedir:str, fp_marbel_basedir:str, settings, verbose=True) -> pd.DataFrame:
    # which environments to plot and in which order
    environments = get_environments(fp_orb_basedir, settings)

    data = dict()
    for environment in tqdm(environments, disable=not verbose, desc='Compiling data for DE Venn diagram'):
        fp_truth = join(fp_orb_basedir, environment, 'refdeseq2', '%s_DESeq2_full_table.tsv' % environment)
        truth = pd.read_csv(fp_truth, sep="\t", index_col=0)
        truth.index = list(map(str, truth.index))
        
        # obtain orthogroup information about genes
        genes = pd.read_csv(join(fp_marbel_basedir, '%s_microbiome' % environment, "summary", "gene_summary.csv"), sep=",").set_index('gene_name')
        # add a column that reports if the gene is part of a single or multi gene orthogroup
        genes = genes.merge(genes.groupby('orthogroup').size().apply(lambda x: 'single_gene_OG' if x == 1 else 'multi_gene_OG').rename('OGsize'),
            left_on='orthogroup', right_index=True, how='left')
        truth = truth.merge(genes[['orthogroup', 'OGsize']],
                            left_index=True, right_index=True, how='left')

        # obtain DE truth if genes are collapsed into orthogroups
        fp_truth_OG = join(fp_ogtruth_basedir, environment, '%s_DESeq2_full_table.tsv' % environment)
        truthOG = pd.read_csv(fp_truth_OG, sep="\t", index_col=0)
        truthOG.index = list(map(str, truthOG.index))
        truth = truth.merge(truthOG, left_on='orthogroup', right_index=True, how='left', suffixes=('_gene', '_orthogroup'))

        # flag elements as DE
        for typ in ['gene', 'orthogroup']:
            truth['DE%s' % typ] = truth['padj_%s' % typ] < 0.05

        data[environment] = truth

    return data


def getdata_rnaquast(fp_orb_basedir:str, fp_quast_basedir:str, settings, verbose=True) -> pd.DataFrame:
    # which environments to plot and in which order
    environments = get_environments(fp_orb_basedir, settings)

    data = []
    for environment in tqdm(environments, disable=not verbose, desc='Compile RNAquast data'):
        quast = pd.read_csv(join(fp_quast_basedir, environment, "short_report.tsv"), sep="\t", index_col=0)
        quast.index.name = 'metric'
        quast.columns = list(map(lambda x: settings['labels']['assemblers'].get(x.split('_contigs')[0], x), quast.columns))
        quast = quast.stack().reset_index().rename(columns={'level_1': 'assembler', 0: 'score'})
        quast['environment'] = environment
        data.append(quast)
    return pd.concat(data).set_index(['environment', 'assembler', 'metric'])
