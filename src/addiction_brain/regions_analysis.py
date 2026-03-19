import pandas as pd




def regions_overlap(add_targets_t1 : pd.DataFrame, add_targets_t2: pd.DataFrame, non_targets_t1: pd.DataFrame, non_targets_t2: pd.DataFrame):
    """
    A function to perform regions analysis on the targets after enrichment analysis and 
    get the top 5 regions for addictive enriched and non-addictive enriched targets. 
    This does the analysis for the tier1 and tier2 level targets to check for the overlapping, addictive 
    specific and non-addictive specific regions for each tier. 
    :param df: The statistical enrichment dataframe with statistical testing results. 
    :type df: pd.DataFrame
    """


    #getting the highest expressing regions for the targets for tier1: 
    top_regions_addt1 = set(
    add_targets_t1.sort_values(["target_symbol", "pTPM"], ascending=[True, False])
    .groupby("target_symbol")
    .head(5)['region_name']
    )

    top_regions_nont1 = set(
        non_targets_t1.sort_values(['target_symbol', 'pTPM'], ascending = [True, False])
        .groupby("target_symbol")
        .head(5)['region_name']
    )

    #getting the highest expressing regions for the targets of tier2: 

    top_regions_addt2  = set(
    add_targets_t2.sort_values(["target_symbol", "pTPM"], ascending=[True, False])
    .groupby("target_symbol")
    .head(5)['region_name']
    )

    top_regions_nont2 = set(
        non_targets_t2.sort_values(['target_symbol', 'pTPM'], ascending = [True, False])
        .groupby("target_symbol")
        .head(5)['region_name']
    )

    #getting the addiction specific and nonaddictive specific and overlapping regions: 

    add_specific = top_regions_addt1 - top_regions_nont1
    non_specific = top_regions_nont1 - top_regions_addt1
    overlap = top_regions_addt1 & top_regions_nont1
    
    add_specifict2 = top_regions_addt2 - top_regions_nont2
    non_specifict2 = top_regions_nont2 - top_regions_addt2
    overlapt2 = top_regions_addt2 & top_regions_nont2

    #making the target list for each of the specific regions, overlapping and the add and non-add specific regions: 
    
    #overlapping regions targets: 
    add_overt1 = add_targets_t1[(add_targets_t1['region_name'].isin(overlap)) & (add_targets_t1['pTPM'] >=1)]
    non_overt1 = non_targets_t1[(non_targets_t1['region_name'].isin(overlap)) & (non_targets_t1['pTPM'] >=1)]

    add_overt2 = add_targets_t2[(add_targets_t2['region_name'].isin(overlapt2)) & (add_targets_t2['pTPM'] >=1)]
    non_overt2 = non_targets_t2[(non_targets_t2['region_name'].isin(overlapt2)) & (non_targets_t2['pTPM'] >=1)]
    
    #specific regions targets: 
    add_regions = add_targets_t1[(add_targets_t1['region_name'].isin(add_specific)) & (add_targets_t1['pTPM'] >=1)]
    non_regions = non_targets_t1[(non_targets_t1['region_name'].isin(non_specific)) & (non_targets_t1['pTPM'] >=1)]

    add_regionst2 = add_targets_t2[(add_targets_t2['region_name'].isin(add_specifict2)) & (add_targets_t2['pTPM'] >=1)]
    non_regionst2 = non_targets_t2[(non_targets_t2['region_name'].isin(non_specifict2)) & (non_targets_t2['pTPM'] >=1)]
    
    
    
    return add_overt1, add_overt2, non_overt1, non_overt2, add_regions, add_regionst2, non_regions, non_regionst2


def run_region(add_targets_t1 : pd.DataFrame, add_targets_t2: pd.DataFrame, non_targets_t1: pd.DataFrame, non_targets_t2: pd.DataFrame): 
    """
    For running the full regional analysis for the enriched targets. 

    This gets the overallping, class specific tiered regions for addcitive and non-addictive classes. 

    Returns: 

    Returns 8 dataframes, with targets in the overlapping, specific regions (tier1 and tier2)
    
    :param df: The enriched dataframe for the targets. 
    :type df: pd.Dataframe
    """

    add_targets_t1, add_targets_t2, non_targets_t1, non_targets_t2 = add_targets_t1.copy(), add_targets_t2.copy(), non_targets_t1.copy(), non_targets_t2.copy()

    add_overt1, add_overt2, non_overt1, non_overt2, add_regions, add_regionst2, non_regions, non_regionst2 = regions_overlap(add_targets_t1, add_targets_t2, non_targets_t1, non_targets_t2)

    return add_overt1, add_overt2, non_overt1, non_overt2, add_regions, add_regionst2, non_regions, non_regionst2

    








