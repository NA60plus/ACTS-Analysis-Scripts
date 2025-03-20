true_match_rate = 0.55
import math
def ratio(x):
    return 1/x**2-1

def get_eff(x):
    return 1/math.sqrt(x+1)
print(ratio(0.55))
print(ratio(0.62))
print(ratio(0.9))

print(get_eff(0.3))

import pandas as pd

# Read the CSV file
df = pd.read_csv("chi2NA60.csv")  # Replace with your actual file path

# Compute the sum of column 'y'
y_sum = df['y'].sum()
chi2max = 1.5
y_sel = 0
for x,y in zip(df['x'], df['y']):
    if x < 1.5:
        y_sel += y

print("Eff (chi2<3):", y_sum)
print("Eff (chi2<1.5):", y_sel)
