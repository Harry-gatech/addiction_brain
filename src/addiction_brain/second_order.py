import pandas as pd 
import numpy as np 


def second_order(df): 
    add_drugs = df[df['true_status'] == 1]
    non_drugs = df[df['true_status'] == 0]

    # Use the pre-filtered DataFrames
    targets_add = set().union(*add_drugs["second_order_hits_non"].dropna())
    targets_non_first = set().union(*non_drugs["first_order_hits_non"].dropna())
    targets_second_non = set().union(*non_drugs["second_order_hits_non"].dropna())
    
    return targets_add, targets_non_first, targets_second_non

def overlap(targets_add, targets_non_first, targets_second_non):
    #get the overlaps for the first and second order addictive and non-addictive targets interactions: 
    overlap_add_non_first = targets_add.intersection(targets_non_first)
    overlap_add_second_non = targets_add.intersection(targets_second_non)

    # check for the percentage of the overlap in the first and second order interactions: 

    first_overlap = len(overlap_add_non_first) / len(targets_non_first) if len(targets_non_first) > 0 else 0
    second_overlap = len(overlap_add_second_non) / len(targets_second_non) if len(targets_second_non) > 0 else 0

    return overlap_add_non_first, overlap_add_second_non, first_overlap, second_overlap

def run_second(df):
    targets_add, targets_non_first, targets_second_non = second_order(df)
    overlap_add_non_first, overlap_add_second_non, first_overlap, second_overlap = overlap(
        targets_add, targets_non_first, targets_second_non
    )

    print(f"Set-level overlap (addictive 2nd-order non-add ∩ non-addictive 1st-order non-add): "
          f"{len(overlap_add_non_first)} proteins ({first_overlap:.2%})")
    print(f"Set-level overlap (addictive 2nd-order non-add ∩ non-addictive 2nd-order non-add): "
          f"{len(overlap_add_second_non)} proteins ({second_overlap:.2%})")

    return overlap_add_non_first,overlap_add_second_non, first_overlap, second_overlap

