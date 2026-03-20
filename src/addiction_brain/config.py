from pathlib import Path

_DATA_DIR = Path(__file__).parent.parent.parent / "Data"

BINDING_PATH       = _DATA_DIR / "Binding.txt"
TRANSCRIPTOMIC_PATH = _DATA_DIR / "transcriptomic_data.tsv"
CTD_PATH           = _DATA_DIR / "CTD_genes_pathways.txt"
STRING_NODE_PATH   = _DATA_DIR / "STRING_node.csv"
STRING_EDGE_PATH   = _DATA_DIR / "STRING_edge.csv"