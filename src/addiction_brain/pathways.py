import pandas as pd 
import numpy as np
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from .config import BINDING_PATH
def load_data(df:pd.DataFrame): 

    """
    Load the data for the engagement scores and everything in between. 
    
    :param df: Description
    :type df: pd.DataFrame
    """

    #load the data for the binding data to merge with rep data to make a drug -> target -> label -> binding dataframe.
    df = df.copy() 

    df = df[['compound_id', 'compound_name', 'addiction_status']]

    bind = pd.read_csv(BINDING_PATH, sep = r"\s+", header = None, names=["compound_id", "binding_prob", "target_symbol"] )

    df = df.merge(bind, on= 'compound_id', how= 'inner')
    return df 

def engagement(df_drug:pd. DataFrame, df_region: pd.DataFrame): 
    """
    Calculate the engagement score for each drug and then summarize at the pathway level with drug dataframe
    and the regional level dataframe for the specific targets and the regions in which they are. 
    
    :param df_drug: Description
    :type df_drug: pd.DataFrame
    :param df_region: Description
    :type df_region: pd.DataFrame

    Returns: 
    Dp: drug-pathway dataframe which contains the specifics for the drug and pathway level information. 
    pathway_summary: the summary for each pathway and what engagement scores they have, this is going to be used for discussion section.
    """

    df_drug = df_drug.copy()
    df_region = df_region.copy()
    df_region_use = df_region[["target_symbol", "pathway_id", "pathway_name", "target_class"]].drop_duplicates()


    df_drug = df_drug[df_drug['binding_prob'] >= 0.70].copy()
    bt = df_drug.merge(df_region_use, on="target_symbol", how="inner")
    dp = (bt.groupby(["compound_id","addiction_status","pathway_id", "pathway_name","target_class"])['target_symbol'].nunique().unstack(fill_value = 0).reset_index())
    if "add" not in dp.columns:
        dp["add"] = 0
    if "non" not in dp.columns:
        dp["non"] = 0

    dp = dp.rename(columns={"add": "t_add", "non": "t_non"})

    denom = (dp["t_add"] + dp["t_non"]).astype(float)
    dp["E_drug"] = np.where(denom > 0, (dp["t_add"] - dp["t_non"]) / denom, np.nan)

    p_comp = (
        df_region_use.groupby(
            ['pathway_id', 'pathway_name', 'target_class']
        )['target_symbol'].nunique().unstack(fill_value = 0).reset_index()
    )
    if "add" not in p_comp.columns:
        p_comp["add"] = 0
    if "non" not in p_comp.columns:
        p_comp["non"] = 0

    p_comp = p_comp.rename(columns = {"add": "n_add_p", "non": "n_non_p"})


    dp = dp.merge(p_comp, on = ['pathway_id', 'pathway_name'], how = 'inner')
    dp["E_comp"] = np.where(
    (dp["n_add_p"] > 0) & (dp["n_non_p"] > 0),
    (dp["t_add"] / dp["n_add_p"]) - (dp["t_non"] / dp["n_non_p"]),
    np.nan)

    

    pathway_summary = (
    dp.groupby(["pathway_id", "pathway_name"])
    .agg(
        mean_E_drug=("E_drug", "mean"),
        mean_E_comp=("E_comp", "mean"),
        frac_tilt_edrug=("E_drug", lambda s: (s > 0).mean()),
        frac_tilt_ecomp = ("E_comp", lambda s: (s > 0).mean()),
        n_drugs_addictive=("addiction_status", lambda x: (x == 1).sum()),
        n_drugs_nonaddictive=("addiction_status", lambda x: (x == 0).sum()),
        n_add_p=("n_add_p", "first"),
        n_non_p=("n_non_p", "first"),
    )
    .reset_index()
    )

    return dp, pathway_summary

def statistical_tests(df: pd.DataFrame): 
    """
    Run the statistical tests for the engagement scores for the drugs and the compositional level data for each pathway to determine significance. 
    :param df: Description
    :type df: pd.DataFrame
    """

    df = df.copy()
    rows = []

    for (pid, pname), sub in df.groupby(['pathway_id', 'pathway_name'], dropna = False): 

        x = sub.loc[sub['addiction_status'] == 1 , 'E_drug'].dropna()
        y = sub.loc[sub['addiction_status'] == 0 , 'E_drug'].dropna()

        if len(x) >= 2 and len(y) >= 2: 
            p_drug = mannwhitneyu(x,y, alternative = 'two-sided').pvalue
        else: 
            p_drug = np.nan

        x1 = sub.loc[sub['addiction_status'] == 1 , 'E_comp'].dropna()
        y1 = sub.loc[sub['addiction_status'] == 0 , 'E_comp'].dropna()

        if len(x1) >= 2 and len(y1) >= 2: 
            p_comp = mannwhitneyu(x1,y1, alternative = 'two-sided').pvalue
        else: 
            p_comp = np.nan

        rows.append({
            "pathway_id": pid,
            "pathway_name": pname,
            "n_add_drugs": int(len(x)),
            "n_non_drugs": int(len(y)),
            "p_value_drug": p_drug,
            "p_value_comp": p_comp,
        })
    stats_df = pd.DataFrame(rows)
    stats_df["FDR_drug"] = np.nan
    m = stats_df["p_value_drug"].notna()
    if m.sum() > 0:
        stats_df.loc[m, "FDR_drug"] = multipletests(stats_df.loc[m, "p_value_drug"].values, method="fdr_bh")[1]

    stats_df["FDR_comp"] = np.nan
    m2 = stats_df["p_value_comp"].notna()
    if m2.sum() > 0:
        stats_df.loc[m2, "FDR_comp"] = multipletests(stats_df.loc[m2, "p_value_comp"].values, method="fdr_bh")[1]

    return stats_df

def run_pathways(df_rep: pd.DataFrame, df_region_data: pd.DataFrame): 
    """
    run the analysis for the engagement scores for the pathways and statistical tests on the engagement scores. 
    This is in concordance with the tiered lists for the overlap and region specific lists. 
    
    :param df: Description
    :type df: pd.Dataframe
    """
    df_rep = df_rep.copy()

    df_region_data = df_region_data.copy()

    reps = load_data(df_rep)

    dp, path_summary = engagement(reps, df_region_data)

    dp_stats = statistical_tests(dp)

    path_summary = path_summary.merge(dp_stats, on = ["pathway_id", "pathway_name"], how = "left")

    return dp, path_summary








    



        



        





















 