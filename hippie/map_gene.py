from mygene import MyGeneInfo

def normalize_token(tok: str) -> str:
    tok = tok.strip()
    if not tok:
        return ""
    # handle "K0240_HUMAN,GSC1L_HUMAN"
    # splitting happens outside too, but keep safe
    tok = tok.replace(" ", "")
    # strip common UniProt-ish suffix
    if tok.endswith("_HUMAN"):
        tok = tok[:-6]
    return tok

# read + normalize
genes = []
with open("protein_T2.txt") as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        for part in line.split(","):
            g = normalize_token(part)
            if g:
                genes.append(g)

# de-duplicate while preserving order
seen = set()
genes = [g for g in genes if not (g in seen or seen.add(g))]

mg = MyGeneInfo()

# Query: return symbols + entrez + name (handy for QC)
results = mg.querymany(
    genes,
    scopes=["symbol", "alias"],   # a bit more forgiving than symbol-only
    fields=["symbol", "entrezgene", "name"],
    species="human",
    as_dataframe=False,
    returnall=False
)

# Write mapping table
with open("gene_symbol_to_entrez.tsv", "w") as f:
    f.write("input\tsymbol\tentrez_id\tname\tstatus\n")
    for r in results:
        inp = r.get("query", "")
        sym = r.get("symbol", "")
        eid = r.get("entrezgene", "")
        name = r.get("name", "")
        status = "ok" if eid else ("notfound" if r.get("notfound") else "missing_entrez")
        f.write(f"{inp}\t{sym}\t{eid}\t{name}\t{status}\n")

# Also write just Entrez list (for HIPPIE jar input etc.)
with open("genes_entrez.txt", "w") as f:
    for r in results:
        if "entrezgene" in r:
            f.write(str(r["entrezgene"]) + "\n")

print("done")
