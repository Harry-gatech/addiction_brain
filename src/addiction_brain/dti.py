import pandas as pd
import numpy as np
import re


def make_data(bind, add, non, rep): 
    
    bind = bind.copy()
    add = add.copy()
    non = non.copy()
    rep = rep.copy()

    bind = bind[bind['binding_prob'] >= 0.7]

    add_targets = add['target_symbol'].unique()
    non_targets = non['target_symbol'].unique()

    rep = rep[['compound_id','compound_name', 'addiction_status']]

    bind = bind[bind['target_symbol'].isin(add_targets) | bind['target_symbol'].isin(non_targets)]

    bind = bind.merge(rep, on='compound_id', how = 'inner')

    return bind

def parse_string_edge(shared_name):
    parts = re.findall(r"9606\.\w+", shared_name)
    if len(parts) == 2:
        return parts[0], parts[1]
    return None, None


def string_data(node, edge, dti): 
    
    node = node.copy()
    edge = edge.copy()

    id_map = dict(zip(node["name"], node["display name"]))
    #extracting the ensp1 and ensp2 from the shared name.

    edge[["ensp1", "ensp2"]] = edge['shared name'].apply(lambda x: pd.Series(parse_string_edge(x)))
    #mapping ensps to the gene ids. 

    edge['source'] = edge['ensp1'].map(id_map)
    edge['target'] = edge['ensp2'].map(id_map)

    string_edge = edge.dropna(subset = ['source', 'target'])

    string_fmt = string_edge[['source', 'target', 'stringdb::score']].rename(columns = {"stringdb::score": "weight"})

    string_fmt['interaction'] = "PPI"
    string_fmt['addiction_status'] = None
    string_fmt['compound_name'] = None

    dti_edges_fmt = dti.rename(columns={
    "compound_id"    : "source",
    "target_symbol"  : "target",
    "binding_prob"   : "weight"
})[["source", "target", "weight", "compound_name", "addiction_status"]]
    
    dti_edges_fmt['interaction'] = "Drug-target"

    combined_edges = pd.concat([string_fmt, dti_edges_fmt], ignore_index=True)

    return combined_edges

def run_dti(add, non, rep): 
    
    bind = pd.read_csv(r"Data\Binding.txt", sep = r"\s+", header = None, names=["compound_id", "binding_prob", "target_symbol"])
    node = pd.read_csv(r"Data\STRING_node.csv")
    edge = pd.read_csv(r"Data\STRING_edge.csv")
    dti = make_data(bind, add, non, rep)
    combined_edges = string_data(node, edge, dti)
    return combined_edges





    




