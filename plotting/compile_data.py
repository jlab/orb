from glob import glob
from os.path import join, basename
import yaml
import numpy as np
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

def getdata_recovery(fp_orb_basedir:str, settings) -> pd.DataFrame:
    """Collects data from Orb for contig recovery analysis.
    
    Parameters
    ----------
    fp_orb_basedir : str
        The filepath to Orb's base data dir
    settings : yaml-dict
        General settings for plotting Orb graphs. Here, we need
        a) the environment to skip

    Returns
    -------
    A pandas.DataFrame that holds performance parameters for assemblers in the different environments.
    """
    environments = get_environments(fp_orb_basedir, settings)

    data = []
    for environment in environments:
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
    
def getdata_gene_recovery(fp_orb_basedir:str, settings) -> pd.DataFrame:
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
    for environment in environments:
        recovered_contigs[environment] = dict()
        for fp_category in glob(join(fp_orb_basedir, environment, 'categorizecontigs', '*_contigs_categorised.tsv')):
            assembler = basename(fp_category).split('_contigs_categorised.tsv')[0]
            
            contig_classes = pd.read_csv(fp_category, sep="\t", index_col=0)
            contig_classes.index = list(map(str, contig_classes.index))

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
    
    recovered_genes = []
    for environment in environments:
        features = pd.concat([pd.Series(index=list(v['gene_name'].unique()), data=True, name=k).rename_axis('gene_name')
                            for k, v in recovered_contigs[environment].items()], axis=1).fillna(False)
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
    recovered_genes = pd.concat(recovered_genes, axis=1)
    recovered_genes.index.name = 'assembler'

    # pretty assembler names
    recovered_genes = recovered_genes.rename(index=settings['labels']['assemblers'])

    return recovered_genes
