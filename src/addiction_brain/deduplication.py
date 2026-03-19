import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold
import re
from rdkit.Chem import SaltRemover
from rdkit.ML.Cluster import Butina


#Standardization helpers
_remover = SaltRemover.SaltRemover()  # default RDKit salt definitions
_lfc = rdMolStandardize.LargestFragmentChooser()


#helper functions: 

def load_mapping(map_file):
    rows = []
    with open(map_file, 'r') as f:
        for line in f:
            # normalize the data: 
            cleaned_line = line.strip().replace('"', '')
            parts = re.split(r'\s+', cleaned_line)
            if len(parts) < 4:
                continue
            smiles = parts[0]
            compound_id = parts[1]

            # Addiction status is first numeric token after compound_id: 
            addiction_status = None
            name_tokens = []
            for token in parts[2:]:
                if addiction_status is None and re.fullmatch(r'\d+', token):
                    addiction_status = int(token)
                elif addiction_status is not None:
                    name_tokens.append(token)
            if addiction_status is None:
                print(f"[WARNING] No addiction status found for line: {line}")
                continue

            compound_name = ' '.join(name_tokens).strip()
            rows.append({
                'compound_id': compound_id,
                'smiles': smiles,
                'addiction_status': addiction_status,
                'compound_name': compound_name
            })
    return pd.DataFrame(rows)

#chemistry functions:

def to_parent_mol(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    # Strip common salts/counterions; avoid empty-mol outcomes
    mol = _remover.StripMol(mol, dontRemoveEverything=True, sanitize=True)
    if mol is None or mol.GetNumAtoms() == 0:
        return None

    # Keep the largest remaining fragment (guarantees a single "parent")
    mol = _lfc.choose(mol)
    if mol is None or mol.GetNumAtoms() == 0:
        return None

    return mol

def canonicalize_parent_smiles(smiles: str):
    mol = to_parent_mol(smiles)
    return Chem.MolToSmiles(mol, canonical=True) if mol else None

def get_inchikey_parent(parent_smiles: str):
    mol = Chem.MolFromSmiles(parent_smiles) if parent_smiles else None
    return Chem.MolToInchiKey(mol) if mol else None

def get_fingerprint_parent(parent_smiles: str):
    mol = Chem.MolFromSmiles(parent_smiles) if parent_smiles else None
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048) if mol else None

def get_murcko_parent(parent_smiles: str):
    mol = Chem.MolFromSmiles(parent_smiles) if parent_smiles else None
    if not mol:
        return None
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    return Chem.MolToSmiles(scaffold, canonical=True) if scaffold else None


#inchikey removal approach for the label conflict (Stage 1): 

def remove_exact_duplicates(df: pd.DataFrame): 
    """
    Docstring for remove_exact_duplicates

    Remove exact duplicates based on inchikey for the DEA and the phase 4 drugs. THis goes for the stage 1
    for clearing and cleaning the dataset. 

    Policy: 

    1. Keep the addictive label and remove the non-addictive one if similar are found. 
    
    Returns: 
    pd.DataFrame: Cleaned_df -> this is the cleaned dataframe after removing the exact duplicates. 
    """

    df = df.copy()

    if "parent_smiles" not in df.columns:
        df['parent_smiles'] = df['smiles'].apply(canonicalize_parent_smiles)
    df = df[df['parent_smiles'].notnull()].copy()

    if "inchikey" not in df.columns:
        df['inchikey'] = df['parent_smiles'].apply(get_inchikey_parent)
    df = df[df['inchikey'].notnull()].copy()

    #inchikey conflict: 

    label_nunique = df.groupby('inchikey')['addiction_status'].nunique()
    conflict_inchikeys = set(label_nunique[label_nunique > 1].index)

    df['mask_label'] = df['inchikey'].isin(conflict_inchikeys)

    drop = df['mask_label'] & (df['addiction_status'] == 0)
    cleaned_df = df[~drop].copy()
    cleaned_df = cleaned_df.drop(columns=['mask_label'])

    return cleaned_df.reset_index(drop=True)

#stage 2: Tanimoto based clustering:


def butina_clusters_fromfps(fingerprints, dist_cutoff=0.2):
    dist = []
    n = len(fingerprints)
    for i in range(1, n):
        sims = DataStructs.BulkTanimotoSimilarity(fingerprints[i], fingerprints[:i])
        dist.extend([1.0 - s for s in sims])
    clusters = Butina.ClusterData(dist, n, dist_cutoff, isDistData=True)
    return sorted(clusters, key=len, reverse=True)



def annotate_butina_clusters(df: pd.DataFrame, clusters): 

    """
    Docstring for annotate_butina_clusters
    
    Annotate butina clusters based on the addiction liability they are showing.
    Add cluster id and cluster type based on the addiction liability.
    """

    df = df.copy()
    cluster_ids = [-1] * len(df)
    cluster_types = [""] * len(df)

    for cluster_id, member in enumerate(clusters):

        labels = set(df.iloc[list(member)]['addiction_status'].unique())

        if labels == {1}:
            cluster_type = "addictive_only"
        elif labels == {0}:
            cluster_type = "non-addictive_only"
        else:
            cluster_type = "mixed"

        for m in member:
            cluster_ids[m] = cluster_id
            cluster_types[m] = cluster_type
        
    df['butina_cluster_id'] = cluster_ids
    df['butina_cluster_type'] = cluster_types
    return df

def select_one_per_label_per_cluster(df: pd.DataFrame) -> pd.DataFrame:
    """
    For each cluster_id:
      - if mixed (labels {0,1}), pick 1 rep from label 0 and 1 rep from label 1
      - else pick 1 rep from the only label present
    Deterministic: first row within each label-group.
    """
    reps = []
    for cid, g in df.groupby("butina_cluster_id", sort=True):
        labels = sorted([int(x) for x in g["addiction_status"].dropna().unique()])

        if len(labels) == 0:
            continue

        # Works for both mixed and non-mixed
        for label in labels:
            reps.append(g[g["addiction_status"] == label].iloc[0].copy())

    return pd.DataFrame(reps).reset_index(drop=True)

#end to end classification: 


def clean_inchi_butina(map_file, tanimoto_threshold = 0.8, compute_murcko: bool = False): 

    """
    End-to-end chemical deduplication and diversity filtering pipeline using
    InChIKey normalization followed by Butina (Tanimoto-based) clustering.

    This function performs a two-stage cleanup of a drug–label mapping file:

    Stage 1 — Exact structure deduplication (InChIKey-based):
        • SMILES are standardized by salt stripping and largest-fragment selection.
        • Parent molecules are canonicalized and converted to InChIKeys.
        • If the same InChIKey appears with conflicting addiction labels,
          the addictive label (1) is retained and the non-addictive label (0) is removed.
        • Remaining exact duplicates are collapsed to a single representative.

    Stage 2 — Structural similarity clustering (Butina clustering):
        • Morgan fingerprints (radius=2, 2048 bits) are computed on parent SMILES.
        • Compounds are clustered using Butina clustering with a Tanimoto
          similarity threshold (default = 0.8).
        • Each cluster is annotated as:
            - 'addictive_only'
            - 'non-addictive_only'
            - 'mixed'
        • One representative per label per cluster is selected deterministically
          (first occurrence within each label).

    The output consists of:
        1) A reduced representative set suitable for model training or analysis.
        2) A fully annotated dataset with cluster assignments for inspection or QC.

    Parameters
    ----------
    map_file : str
        Path to the compound mapping file containing SMILES, compound IDs,
        addiction labels, and compound names.
    tanimoto_threshold : float, optional (default=0.8)
        Tanimoto similarity threshold used for Butina clustering.
        Internally converted to a distance cutoff as (1 - threshold).
    compute_murcko : bool, optional (default=False)
        If True, compute and store Murcko scaffolds for each parent structure.

    Returns
    -------
    reps : pandas.DataFrame
        Representative compounds after InChIKey cleanup and one-per-label
        selection within each Butina cluster.
    annotated : pandas.DataFrame
        Full cleaned dataset annotated with Butina cluster IDs and cluster types.

    Notes
    -----
    • The pipeline is deterministic given fixed input ordering.
    • Designed for label-conflict resolution and scaffold-aware dataset reduction
      prior to machine learning or statistical analysis.
    • Fingerprint columns are removed from outputs to reduce memory footprint.
    """

    df = load_mapping(map_file)

    cleaned = remove_exact_duplicates(df)
    cleaned = cleaned.drop_duplicates(subset="inchikey", keep="first").reset_index(drop=True)

    # fingerprints
    cleaned["fingerprint"] = cleaned["parent_smiles"].apply(get_fingerprint_parent)
    cleaned = cleaned[cleaned["fingerprint"].notnull()].reset_index(drop=True)
    
    if compute_murcko:
        cleaned["murcko"] = cleaned["parent_smiles"].apply(get_murcko_parent)

    # Butina clustering (distance cutoff)
    dist_cutoff = 1.0 - float(tanimoto_threshold)
    clusters = butina_clusters_fromfps(cleaned["fingerprint"].tolist(), dist_cutoff=dist_cutoff)


    # Annotate
    annotated = annotate_butina_clusters(cleaned, clusters)

    # One-per-label representatives per cluster
    reps = select_one_per_label_per_cluster(annotated)

    # drop heavy fp column for downstream
    annotated = annotated.drop(columns=["fingerprint"], errors="ignore")
    reps = reps.drop(columns=["fingerprint"], errors="ignore")

    # QC
    print("=== InChIKey conflict cleanup ===")
    print(f"Loaded: {len(df)}")
    #print(f"Conflicted inchikeys: {len(conflict_inchikeys)}")
    #print(f"Dropped (label=0 in conflicted inchikeys): {len(dropped)}")
    print(f"Remaining after cleanup: {len(cleaned)}")
    print()
    print("=== Butina clustering ===")
    #print(f"Tanimoto threshold: {tanimoto_threshold} (dist_cutoff={dist_cutoff})")
    print(f"Clusters: {len(clusters)}")
    print(f"Representatives (one-per-label per cluster): {len(reps)}")
    print("Cluster types in reps:")
    print(reps["butina_cluster_type"].value_counts(dropna=False))

    return reps, annotated








