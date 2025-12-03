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
    
