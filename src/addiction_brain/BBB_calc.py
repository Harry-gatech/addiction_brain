# src/addiction_brain/BBB_calc.py
# RDKit-only BBB permeability annotation utilities (no external models)

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski
from rdkit.Chem.rdMolDescriptors import CalcTPSA



# RDKit mol + properties


def mol_from_smiles(smiles: str) -> Optional[Chem.Mol]:
    """Safe RDKit MolFromSmiles wrapper."""
    if not isinstance(smiles, str) or not smiles.strip():
        return None
    return Chem.MolFromSmiles(smiles)


def compute_bbb_properties_from_smiles(smiles: str) -> dict:
    """
    Compute a standard set of RDKit descriptors used as BBB proxies.
    Expects a single SMILES (ideally your parent/canonical SMILES).

    Returns keys:
      MW, TPSA, HBD, HBA, RotB, clogP, rings, fr_ArRings
    """
    mol = mol_from_smiles(smiles)
    if mol is None:
        return {
            "MW": None,
            "TPSA": None,
            "HBD": None,
            "HBA": None,
            "RotB": None,
            "clogP": None,
            "rings": None,
            "fr_ArRings": None,
        }

    return {
        "MW": float(Descriptors.MolWt(mol)),
        "TPSA": float(CalcTPSA(mol)),
        "HBD": int(Lipinski.NumHDonors(mol)),
        "HBA": int(Lipinski.NumHAcceptors(mol)),
        "RotB": int(Lipinski.NumRotatableBonds(mol)),
        "clogP": float(Crippen.MolLogP(mol)),
        "rings": int(Lipinski.RingCount(mol)),
        "fr_ArRings": int(Lipinski.NumAromaticRings(mol)),
    }



# BBB scoring: simple heuristic


@dataclass(frozen=True)
class BBBHeuristicParams:
    """
    Deterministic BBB-likeness proxy (0..5).
    Interpretable and reproducible; not a perfect BBB classifier.
    """
    max_mw: float = 480.0
    max_tpsa: float = 90.0
    max_hbd: int = 2
    max_rotb: int = 10
    logp_range: Tuple[float, float] = (0.0, 5.0)
    min_score_for_permeable: int = 3


def bbb_heuristic_score(
    MW: Optional[float],
    TPSA: Optional[float],
    HBD: Optional[int],
    RotB: Optional[int],
    clogP: Optional[float],
    params: BBBHeuristicParams = BBBHeuristicParams(),
) -> Tuple[Optional[int], Optional[int]]:
    """
    Returns (score, permeable_flag).
    If any required input is missing, returns (None, None).
    """
    if MW is None or TPSA is None or HBD is None or RotB is None or clogP is None:
        return None, None

    score = 0
    if MW <= params.max_mw:
        score += 1
    if TPSA <= params.max_tpsa:
        score += 1
    if HBD <= params.max_hbd:
        score += 1
    if RotB <= params.max_rotb:
        score += 1
    if params.logp_range[0] <= clogP <= params.logp_range[1]:
        score += 1

    permeable = 1 if score >= params.min_score_for_permeable else 0
    return int(score), int(permeable)


# BOILED-Egg-ish gate (simple)


@dataclass(frozen=True)
class EggGateParams:
    """
    Transparent approximation of a BOILED-Egg-like gate in TPSA/logP space.
    Not the exact SwissADME ellipse; useful as a second view.
    """
    max_tpsa: float = 90.0
    logp_min: float = 0.0
    logp_max: float = 5.0


def egg_gate_label(
    TPSA: Optional[float],
    clogP: Optional[float],
    params: EggGateParams = EggGateParams(),
) -> Tuple[Optional[int], Optional[int]]:
    """Returns (egg_pass, bbb_permeable_egg). If missing values, returns (None, None)."""
    if TPSA is None or clogP is None:
        return None, None
    inside = (TPSA <= params.max_tpsa) and (params.logp_min <= clogP <= params.logp_max)
    val = 1 if inside else 0
    return val, val


# Public: annotate a DataFrame


def annotate_bbb(
    df: pd.DataFrame,
    smiles_col: str = "parent_smiles",
    heuristic_params: BBBHeuristicParams = BBBHeuristicParams(),
    egg_params: EggGateParams = EggGateParams(),
    prefix: str = "bbb_",
) -> pd.DataFrame:
    """
    Adds BBB-related columns to df:
      {prefix}MW, TPSA, HBD, HBA, RotB, clogP, rings, fr_ArRings
      {prefix}heuristic_score
      {prefix}permeable_heuristic
      {prefix}permeable_egg

    Designed to run after your parent SMILES standardization step.
    """
    if smiles_col not in df.columns:
        raise ValueError(f"'{smiles_col}' not found in df columns")

    out = df.copy()

    props_series = out[smiles_col].apply(compute_bbb_properties_from_smiles)
    props_df = pd.DataFrame(list(props_series))

    # attach with prefix
    for col in props_df.columns:
        out[f"{prefix}{col}"] = props_df[col].values

    # heuristic score + label
    hs = out.apply(
        lambda r: bbb_heuristic_score(
            r[f"{prefix}MW"],
            r[f"{prefix}TPSA"],
            r[f"{prefix}HBD"],
            r[f"{prefix}RotB"],
            r[f"{prefix}clogP"],
            params=heuristic_params,
        ),
        axis=1,
    )
    out[f"{prefix}heuristic_score"] = [x[0] for x in hs]
    out[f"{prefix}permeable_heuristic"] = [x[1] for x in hs]

    # egg gate label
    eg = out.apply(
        lambda r: egg_gate_label(
            r[f"{prefix}TPSA"],
            r[f"{prefix}clogP"],
            params=egg_params,
        ),
        axis=1,
    )
    out[f"{prefix}permeable_egg"] = [x[1] for x in eg]

    return out


def bbb_summary(df: pd.DataFrame, label_col: str = "bbb_permeable_heuristic") -> pd.DataFrame:
    """
    Convenience table: counts by (addiction_status, BBB label).
    """
    if "addiction_status" not in df.columns:
        raise ValueError("df must contain 'addiction_status'")
    if label_col not in df.columns:
        raise ValueError(f"df must contain '{label_col}'")

    return (
        df.groupby(["addiction_status", label_col])
          .size()
          .rename("count")
          .reset_index()
          .sort_values(["addiction_status", label_col])
    )


def bbb_qc_report(df: pd.DataFrame, prefix: str = "bbb_") -> dict:
    """
    Quick QC metrics for sanity checking the descriptor ranges.
    Returns a small dict you can print/log.
    """
    cols = [f"{prefix}MW", f"{prefix}TPSA", f"{prefix}HBD", f"{prefix}HBA", f"{prefix}RotB", f"{prefix}clogP"]
    present = [c for c in cols if c in df.columns]
    report = {
        "n_rows": int(len(df)),
        "n_missing_any_core": int(df[present].isna().any(axis=1).sum()) if present else None,
    }
    for c in present:
        s = df[c]
        report[f"{c}_missing"] = int(s.isna().sum())
        if s.notna().any():
            report[f"{c}_min"] = float(s.min())
            report[f"{c}_max"] = float(s.max())
    return report
