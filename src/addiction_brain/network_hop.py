import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

def classify_ratio_label(add_ratio, non_ratio, n_classified,
                         add_thresh=0.7, non_thresh=0.7,
                         no_data_label="no_classified_targets",
                         add_label="addictive_dominant",
                         non_label="non_addictive_dominant",
                         mixed_label="mixed"):
    """
    Convert ratios into a continous form. 
    """
    if n_classified == 0:
        return no_data_label
    elif add_ratio >= add_thresh:
        return add_label
    elif non_ratio >= non_thresh:
        return non_label
    else:
        return mixed_label


def classify_drug_trajectory(df, add_enriched, non_enriched,
                             add_thresh=0.7, non_thresh=0.7):
    """
    For each drug:
    - classify 1st-order target composition using addictive/non-addictive ratios
    - classify 2nd-order PPI neighborhood composition using ratios
    - keep raw counts and bias scores
    """

    df = df.copy()

    ppi = df[df['interaction'] == 'PPI'].copy()
    drugs = df[df['interaction'] == 'Drug-target'].copy()

    add_enriched = set(add_enriched['target_symbol'].dropna().unique())
    non_enriched = set(non_enriched['target_symbol'].dropna().unique())

    results = []

    for drug in drugs['source'].unique():
        drug_rows = drugs[drugs['source'] == drug]
        status = drug_rows['addiction_status'].iloc[0]

        # direct targets of this drug
        direct = set(drug_rows['target'].dropna().unique())

       
        # 1st-order classification
       
        hits_add = direct & add_enriched #direct targets that are in the addictive enriched sets. 
        hits_non = direct & non_enriched #direct targets that are in the non-addictive enriched sets. 

        n_add_targets = len(hits_add)
        n_non_targets = len(hits_non)
        n_direct_total = len(direct)

        # denominator = classified targets only
        n_direct_classified = n_add_targets + n_non_targets

        if n_direct_classified == 0:
            first_order_add_ratio = np.nan
            first_order_non_ratio = np.nan
            first_order_bias = np.nan
        else:
            first_order_add_ratio = n_add_targets / n_direct_classified
            first_order_non_ratio = n_non_targets / n_direct_classified
            first_order_bias = first_order_add_ratio - first_order_non_ratio

        first_order_label = classify_ratio_label(
            add_ratio=first_order_add_ratio if not np.isnan(first_order_add_ratio) else 0,
            non_ratio=first_order_non_ratio if not np.isnan(first_order_non_ratio) else 0,
            n_classified=n_direct_classified,
            add_thresh=add_thresh,
            non_thresh=non_thresh,
            no_data_label="no_classified_targets",
            add_label="addictive_dominant",
            non_label="non_addictive_dominant",
            mixed_label="mixed"
        )

        
        # 2nd-order classification
        # STRING PPIs are undirected
       
        so_from_source = set(ppi.loc[ppi['source'].isin(direct), 'target'].dropna().unique())
        so_from_target = set(ppi.loc[ppi['target'].isin(direct), 'source'].dropna().unique())

        second_order = (so_from_source | so_from_target) - direct

        so_add = second_order & add_enriched #second order neighbors that are in the addictive enriched sets. 
        so_non = second_order & non_enriched #second order neighbors that are in the non-addictive enriched sets. 

        n_so_add = len(so_add)
        n_so_non = len(so_non)
        n_second_total = len(second_order)

        # denominator = classified second-order neighbors only
        n_second_classified = n_so_add + n_so_non

        if n_second_classified == 0:
            second_order_add_ratio = np.nan
            second_order_non_ratio = np.nan
            second_order_bias = np.nan
        else:
            second_order_add_ratio = n_so_add / n_second_classified
            second_order_non_ratio = n_so_non / n_second_classified
            second_order_bias = second_order_add_ratio - second_order_non_ratio

        second_order_label = classify_ratio_label(
            add_ratio=second_order_add_ratio if not np.isnan(second_order_add_ratio) else 0,
            non_ratio=second_order_non_ratio if not np.isnan(second_order_non_ratio) else 0,
            n_classified=n_second_classified,
            add_thresh=add_thresh,
            non_thresh=non_thresh,
            no_data_label="no_classified_neighbors",
            add_label="addictive_neighborhood",
            non_label="non_addictive_neighborhood",
            mixed_label="mixed_neighborhood"
        )

        trajectory = f"{first_order_label} → {second_order_label}"

        results.append({
            'drug': drug,
            'drug_name': drug_rows['compound_name'].iloc[0],
            'true_status': status,

            #targets found in the enriched sets for the drugs:
            'first_order_hits_add': hits_add, 
            'first_order_hits_non': hits_non,
            'second_order_hits_add': so_add,
            'second_order_hits_non': so_non,

            # labels:
            'first_order': first_order_label,
            'second_order': second_order_label,
            'trajectory': trajectory,

            # first-order counts:
            'n_direct_total': n_direct_total,
            'n_add_targets': n_add_targets,
            'n_non_targets': n_non_targets,
            'n_direct_classified': n_direct_classified,

            # first-order ratios:
            'first_order_add_ratio': first_order_add_ratio,
            'first_order_non_ratio': first_order_non_ratio,
            'first_order_bias': first_order_bias,

            # second-order counts:
            'n_second_total': n_second_total,
            'n_so_add': n_so_add,
            'n_so_non': n_so_non,
            'n_second_classified': n_second_classified,

            # second-order ratios:
            'second_order_add_ratio': second_order_add_ratio,
            'second_order_non_ratio': second_order_non_ratio,
            'second_order_bias': second_order_bias,
        })

    results_df = pd.DataFrame(results)

    summary = (
        results_df
        .groupby(['trajectory', 'true_status'])
        .size()
        .unstack(fill_value=0)
        .sort_values(by=list(results_df['true_status'].dropna().unique()), ascending=False)
    )

    print(summary)

    return results_df, summary


def trajectory_fisher(results_df):
    
    total_add = (results_df['true_status'] == 1).sum()  # full universe
    total_non = (results_df['true_status'] == 0).sum()
    
    summary = results_df.groupby(['trajectory', 'true_status']).size().unstack(fill_value=0)
    
    rows = []
    for trajectory, row in summary.iterrows():
        n_add = row.get(1, 0)
        n_non = row.get(0, 0)
        
        table = [
            [n_add,     total_add - n_add],
            [n_non,     total_non - n_non]
        ]
        or_, p = fisher_exact(table, alternative='two-sided')
        rows.append({
            'trajectory':           trajectory,
            'n_addictive':          n_add,
            'n_non_addictive':      n_non,
            'odds_ratio':           or_,
            'p_value':              p,
            'precision_addictive':  n_add / (n_add + n_non) if (n_add + n_non) > 0 else np.nan
        })
    
    result = pd.DataFrame(rows)
    result['FDR'] = multipletests(result['p_value'], method='fdr_bh')[1]
    return result.sort_values('FDR')

def run_trajectory(network_df, add, non):
    trajectory_df, summary = classify_drug_trajectory(network_df, add, non)
    fisher_results = trajectory_fisher(trajectory_df)
    return trajectory_df, summary, fisher_results