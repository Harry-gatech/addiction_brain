import pandas as pd 
import numpy as np 
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

def load(rep: pd.DataFrame, add: pd.DataFrame, non: pd.DataFrame): 
    """
    A helper function to check for the data and make the dataframes to calculate 
    the ratios based on the drugs and what they are hitting
    """
    drugs = rep[['compound_id', 'compound_name', 'addiction_status']]

    bind = pd.read_csv(r"Data\Binding.txt", sep=r"\s+", header=None, 
                       names=["compound_id", "binding_prob", "target_symbol"])

    add_targets = bind[(bind['target_symbol'].isin(add['target_symbol'])) & (bind['binding_prob'] >= 0.7)]
    non_targets = bind[(bind['target_symbol'].isin(non['target_symbol'])) & (bind['binding_prob'] >= 0.7)]

    return drugs, add_targets, non_targets


def targets_per_drug(drugs, add_targets, non_targets): 

    add_drugs = drugs[drugs['addiction_status'] == 1]
    non_drugs = drugs[drugs['addiction_status'] == 0]

    add_drugs_adtargets = add_drugs.merge(add_targets, on='compound_id', how='left')
    add_drugs_nontargets = add_drugs.merge(non_targets, on='compound_id', how='left')
    non_drugs_adtargets = non_drugs.merge(add_targets, on='compound_id', how='left')
    non_drugs_nontargets = non_drugs.merge(non_targets, on='compound_id', how='left')

    return add_drugs_adtargets, add_drugs_nontargets, non_drugs_adtargets, non_drugs_nontargets


def ratios(add_drugs_adtargets, add_drugs_nontargets, non_drugs_adtargets, non_drugs_nontargets): 
    """
    Calculates the ratio (addictive target count / non-addictive target count + 1) per drug.
    Merges counts by compound_id before dividing to ensure alignment.
    """

    # Count targets per compound for each group
    add_count = (add_drugs_adtargets.groupby('compound_id')['target_symbol']
                 .count().reset_index(name='add_target_count'))
    
    non_add_count = (add_drugs_nontargets.groupby('compound_id')['target_symbol']
                     .count().reset_index(name='non_add_target_count'))

    non_drug_add_count = (non_drugs_adtargets.groupby('compound_id')['target_symbol']
                          .count().reset_index(name='add_target_count'))
    
    non_drug_non_count = (non_drugs_nontargets.groupby('compound_id')['target_symbol']
                          .count().reset_index(name='non_add_target_count'))

    # Merge counts by compound_id before calculating ratio
    add_merged = add_count.merge(non_add_count, on='compound_id', how='left')
    add_merged['non_add_target_count'] = add_merged['non_add_target_count'].fillna(0)
    add_merged['ratio'] = add_merged['add_target_count'] / (add_merged['non_add_target_count'] + 1)

    non_merged = non_drug_add_count.merge(non_drug_non_count, on='compound_id', how='left')
    non_merged['non_add_target_count'] = non_merged['non_add_target_count'].fillna(0)
    non_merged['ratio'] = non_merged['add_target_count'] / (non_merged['non_add_target_count'] + 1)

    return add_merged, non_merged

def statistical_tests(add_merged, non_merged):
    """
    Tests whether the ratio of addictive to non-addictive target hits
    differs significantly between addictive and non-addictive drugs.
    Uses Mann-Whitney U since ratios are not normally distributed.
    """
    add_ratios = add_merged['ratio'].dropna()
    non_ratios = non_merged['ratio'].dropna()

    stat, p = mannwhitneyu(add_ratios, non_ratios, alternative='two-sided')

    print(f"Addictive drugs   — mean ratio: {add_ratios.mean():.3f}, median: {add_ratios.median():.3f}, n={len(add_ratios)}")
    print(f"Non-addictive drugs — mean ratio: {non_ratios.mean():.3f}, median: {non_ratios.median():.3f}, n={len(non_ratios)}")
    print(f"Mann-Whitney U statistic: {stat:.1f}, p-value: {p:.4e}")

    return {"stat": stat, "p_value": p, "add_mean": add_ratios.mean(), "non_mean": non_ratios.mean()}


def run_ratios(rep, add, non):
    drugs, add_targets, non_targets = load(rep, add, non)
    add_drugs_adtargets, add_drugs_nontargets, non_drugs_adtargets, non_drugs_nontargets = targets_per_drug(drugs, add_targets, non_targets)
    add_ratios, non_ratios = ratios(add_drugs_adtargets, add_drugs_nontargets, non_drugs_adtargets, non_drugs_nontargets)
    stats = statistical_tests(add_ratios, non_ratios)
    return add_ratios, non_ratios, stats