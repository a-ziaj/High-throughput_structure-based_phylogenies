import pandas as pd
from scipy.stats import spearmanr


df = pd.read_csv('results.csv')
df['TM_score'] = pd.to_numeric(df['TM_score'], errors='coerce')
clean_df = df.dropna(subset=['TM_score', 'FATCAT_score'])
correlation, p_value = spearmanr(clean_df['TM_score'], clean_df['FATCAT_score'])

print(f"valid pairs: {len(clean_df)}")
print(f"spearmans rho correlation coef: {correlation:.4f}")
print(f"p-value: {p_value:.4g}")

