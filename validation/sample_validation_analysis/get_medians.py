import pandas as pd
import os

df1 = pd.read_csv('PLG_Stage1_Significance.csv')

print(f"Activation Median: {df1['activation'].median():.2f}")
print(f"Occurrence Median: {df1['occurrence_count'].median():.2f}")
