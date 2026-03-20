import numpy as np
import pandas as pd
from collections import defaultdict
import glob, os
from tqdm import tqdm
from pathlib import Path
from collections import defaultdict
from .config import BINDING_PATH, TRANSCRIPTOMIC_PATH, CTD_PATH

#Files for the analysis: 

def load_data(df: pd.DataFrame): 

    """
    
    A helper function to load all the necessary data files for enrichment analysis and doing the statistical tests. 
    Input: A dataframe with drug information and the correct labels for the drug which passed the permeability filter. 
    The input is filtered data for the BBB permeability.  


    """

    df = df.copy()
    df.rename(columns={"addiction_status" : "label"}, inplace=True)

    N_add_total = df[df['label'] == 1]['compound_id'].nunique()
    N_non_total = df[df['label'] == 0]['compound_id'].nunique()

    hpa = pd.read_csv(TRANSCRIPTOMIC_PATH, sep="\t")
    ctd = pd.read_csv(CTD_PATH, sep=r"\s*\|\s*", engine="python").rename(
    columns={"GeneSymbol":"target_symbol", "PathwayID":"pathway_id", "PathwayName":"pathway_name"})
    bind = pd.read_csv(BINDING_PATH, sep = r"\s+", header = None, names=["compound_id", "binding_prob", "target_symbol"] )

    #block for standardizing column names and selecting genes based on pTPM >= 1 in brain regions: 

    expr = hpa.rename(columns={"Gene name":"target_symbol", "Subregion":"region_name"})
    #expression filter for the brain data: 
    expr = expr[expr["pTPM"] >= 1]

    #filtering the binding data for only high confidence targets:

    
    bind = bind.merge(df, on ="compound_id", how="left")
    bind_keep = bind[bind["binding_prob"] >= 0.7]

    #keeping the binding data for the targets present in the expression data:
    target = expr["target_symbol"].unique()
    bind_expr = bind_keep[bind_keep["target_symbol"].isin(target)]

    return bind_expr, expr, ctd, N_add_total, N_non_total

#helper functions for the polyphramacology analysis: 

def polypharamacology(bind_expr: pd.DataFrame): 
    """
    Docstring for polypharamacology, 
    To check for the per drug unique targets
    
    :param bind_expr: Data from the filtered binding, expressed and BBB filtered drugs
    :type bind_expr: pd.DataFrame
    """

    per_drug = (bind_expr.groupby(['compound_id', 'label'])['target_symbol'].nunique().reset_index(name='n_targets'))

    return per_drug

#checking the jaccaard similarity for in class specific tests for addictive and non-addictive drugs: 



def jaccard_similarity(bind_expr: pd.DataFrame, max_pairs: int = 200000, random_state: int = 0) -> pd.DataFrame:
    """
    Pairwise Jaccard similarity between drugs within each label class.
    Caps number of pairs per class to max_pairs for speed/stability.
    """
    rng = np.random.default_rng(random_state)

    drug_targets = (
        bind_expr.groupby(['compound_id', 'label'])['target_symbol']
        .agg(lambda x: set(x))
        .reset_index()
        .rename(columns={'target_symbol': 'targets'})
    )
    
    out = []
    for label, sub in drug_targets.groupby("label"):
        compounds = sub["compound_id"].to_numpy()
        sets = sub["targets"].to_list()
        n = len(compounds)

        if n < 2:
            continue

        total_pairs = n * (n - 1) // 2

        if total_pairs <= max_pairs:
            pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
        else:
            # sample pairs of indices
            i = rng.integers(0, n, size=max_pairs)
            j = rng.integers(0, n, size=max_pairs)
            mask = i != j
            i, j = i[mask], j[mask]
            a = np.minimum(i, j)
            b = np.maximum(i, j)
            pairs = np.unique(np.stack([a, b], axis=1), axis=0)

        for i, j in pairs:
            Ti, Tj = sets[int(i)], sets[int(j)]
            inter = len(Ti & Tj)
            union = len(Ti | Tj)
            out.append((compounds[int(i)], compounds[int(j)], label, inter / union if union else 0.0))

    return pd.DataFrame(out, columns=["compound_id_1", "compound_id_2", "label", "jaccard_index"])


def enrichment(bind_expr, N_add_total, N_non_total, threshold = 0.7): 
    """
    
    A helper function to perform enrichment analysis on the binding data.

    Input: A dataframe with binding data and a threshold for the binding probability.
    """

    #binding probability threshold: 
    df = bind_expr[bind_expr["binding_prob"] >= threshold]
    
    N_add = N_add_total
    N_non = N_non_total

    counts = (
        df.groupby(['target_symbol', 'label'])['compound_id']
        .nunique()
        .unstack(fill_value = 0)
        .rename(columns = {1: "n_add", 0: "n_non"})
    ).reset_index()

    counts['N_add'] = N_add
    counts["N_non"] = N_non

    counts["EF"] = ((counts['n_add'] + 0.5 ) / (counts['N_add'] + 1.0 )) / ((counts['n_non'] + 0.5 ) / (counts['N_non'] + 1.0 ))

    return counts.sort_values("EF", ascending = False).reset_index(drop = True)

def statistical_test(enrich_df, bind_expr): 
    """
    
    A helper function to perform statistical tests on the enrichment results. 

    Input: Enrichment dataframe and binding dataframe. 
    """

    from scipy.stats import fisher_exact
    from statsmodels.stats.multitest import multipletests

    pvals = []
    oddsratios = []
    for _, row in enrich_df.iterrows():
        a = row['n_add']
        b = row['N_add'] - row['n_add']
        c = row['n_non']
        d = row['N_non'] - row['n_non']
        oddsratio, p = fisher_exact([[a, b], [c, d]], alternative='two-sided')
        oddsratios.append(oddsratio)
        pvals.append(p)

    enrich_df['p_value'] = pvals
    enrich_df['odds_ratio'] = oddsratios
    enrich_df['FDR'] = multipletests(enrich_df['p_value'], method='fdr_bh')[1]

    return enrich_df.sort_values('FDR').reset_index(drop=True)


#function to check for the tiers of the targets based on the FDR and the number of addictive and non-addictive drugs hitting: 

def targets(df: pd.DataFrame, sel_targets: int):
    """
    A helper function to extract the targets and their regions after enrichment analysis and statustical testing. 
    This function filters the targets in based on the EF and FDR values. 
    
    And returns a dataframe with the target symbols which are either addictive enriched, or non-addictive enriched. 
    
    :param df: The enrichement dataframe after statistical testing. 

    :type df: pd.Dataframe
    """ 

    df = df.copy()
    #tier2 targets for pathways and network analysis for hypothesis generation: 
    add_targets_t2 = df[(df['EF'] > 2) & (df['n_add'] >= sel_targets)]
    non_targets_t2 = df[(df['EF'] < 0.5) & (df['n_non'] >= sel_targets)]

    #tier1 targets with statistical tests and EF filters: 

    add_targets_t1 = df[(df['EF'] > 2) & (df['FDR'] <= 0.05)]
    non_targets_t1 = df[(df['EF'] < 0.5) & (df['FDR'] <= 0.05)]

    return add_targets_t1, add_targets_t2, non_targets_t1, non_targets_t2



def final_targets_run(df: pd.DataFrame, sel_targets: int): 
    """
    This function returns the tiered list with the enriched targets and their regions: 
    
    :param df: The permeability representative dataframe after BBB filtering.
    :type df: pd.DataFrame
    """
    bind_expr, expr, ctd, N_add_total, N_non_total = load_data(df)
    hits = bind_expr.copy()

    #metric1 : polypharmacology: 
    per_drug = polypharamacology(bind_expr)

    jacc_df = jaccard_similarity(bind_expr)

    print("BBB drug universe:", N_add_total, N_non_total)
    print("Drugs with ≥1 retained interaction:",
      hits.loc[hits["label"]==1, "compound_id"].nunique(),
      hits.loc[hits["label"]==0, "compound_id"].nunique())
    
    print("unique addictive targets:", hits.loc[hits['label'] == 1, 'target_symbol'].nunique())
    print("unique non-addictive targets:", hits.loc[hits['label'] == 0, 'target_symbol'].nunique())
    add_per_drug = bind_expr[bind_expr["label"]==1].groupby("compound_id")["target_symbol"].nunique()
    non_per_drug = bind_expr[bind_expr["label"]==0].groupby("compound_id")["target_symbol"].nunique()

    print("Mean targets per addictive drug:", add_per_drug.mean())
    print("Mean targets per non-addictive drug:", non_per_drug.mean())
    # average target reuse frequency

    target_freq_add = bind_expr[bind_expr["label"]==1].groupby("target_symbol")["compound_id"].nunique()
    target_freq_non = bind_expr[bind_expr["label"]==0].groupby("target_symbol")["compound_id"].nunique()

    print("Mean drugs per target (addictive):", target_freq_add.mean())
    print("Mean drugs per target (non-addictive):", target_freq_non.mean())

    enrich_df = enrichment(bind_expr, N_add_total, N_non_total, threshold = 0.7)
    final_enrich = statistical_test(enrich_df, bind_expr)
    final_enrich_path = final_enrich.merge(ctd[['target_symbol', 'pathway_id', 'pathway_name']].drop_duplicates(), on = 'target_symbol', how = 'left')
    final_enrich_path = final_enrich_path.merge(
        expr[['target_symbol', 'pTPM', 'region_name']].drop_duplicates(),
        on='target_symbol',
        how='left'
    )

    add_targets_t1, add_targets_t2, non_targets_t1, non_targets_t2 = targets(final_enrich_path, sel_targets)

    return add_targets_t1, add_targets_t2, non_targets_t1, non_targets_t2, per_drug, jacc_df










    







