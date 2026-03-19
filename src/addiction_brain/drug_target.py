import pandas as pd

def drug_target_mapping(bind_path: str, rep: pd.DataFrame, prob_thresh: float = 0.7) -> pd.DataFrame:
    # Load bindings file (no header)
    bind_df = pd.read_csv(
        bind_path,
        sep=r"\s+",
        header=None,
        names=["compound_id", "prob", "target_symbol"],
        engine="python"
    )

    # Filter
    bind_filt = bind_df[bind_df["prob"] >= prob_thresh].copy()

    # Merge name + status
    bind_drug = bind_filt.merge(
        rep[["compound_id", "compound_name", "addiction_status"]],
        on="compound_id",
        how="inner"
    )

    return bind_drug


def run_mapping(bind_path: str, rep: pd.DataFrame) -> pd.DataFrame:
    return drug_target_mapping(bind_path, rep)
