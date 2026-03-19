import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns 

def load_data(drugs_df: pd.DataFrame)-> pd.DataFrame: 
    
    bind = pd.read_csv(r"Data\Binding.txt", sep = r"\s+", header = None, names=["compound_id", "binding_prob", "target_symbol"] )
    bind_drug = drugs_df.merge(bind, on = 'compound_id', how ='left')

    expr = pd.read_csv(r"Data\transcriptomic_data.tsv", sep = '\t').rename(columns ={"Gene name": "target_symbol"})

    expr = expr[expr['pTPM'] >= 1]

    bind_expr = bind_drug[bind_drug['target_symbol'].isin(expr['target_symbol'].unique())]

    return bind_expr

def histogram(bind_expr: pd.DataFrame): 
    
    bind_expr = bind_expr.copy()

    print(bind_expr.columns)

def plot_fig(drugs_df: pd.DataFrame): 
    
    bind_expr = load_data(drugs_df)

    histogram(bind_expr)
    


