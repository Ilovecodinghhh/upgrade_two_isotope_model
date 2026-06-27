import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from scipy.stats import linregress

# ============================================================
# SSA decomposition
# ============================================================
def ssa_decompose(ts, window):
    """Return reconstructed components using SSA."""
    N = len(ts)
    K = N - window + 1
    X = np.column_stack([ts[i:i+K] for i in range(window)])
    pca = PCA()
    T = pca.fit_transform(X)
    X_recon = pca.inverse_transform(T)
    recon = np.zeros(N)
    counts = np.zeros(N)

    for i in range(window):
        recon[i:i+K] += X_recon[:, i]
        counts[i:i+K] += 1

    return recon / counts

# ============================================================
# Extract seasonal component (monthly data)
# ============================================================
def extract_seasonal_component(ts):
    """
    For monthly data, seasonal band is ~6–14 months.
    Approximate using SSA with 24-month window.
    """
    window = 24
    return ssa_decompose(ts, window)

# ============================================================
# Main script
# ============================================================
folder = os.path.dirname(os.path.abspath(__file__))
output_folder = os.path.join(folder, "sca_trends_plots")
os.makedirs(output_folder, exist_ok=True)

summary_rows = []

all_files = [f for f in os.listdir(folder) if f.endswith(".txt")]

for file in all_files:
    path = os.path.join(folder, file)
    df = pd.read_csv(path, comment="#", delim_whitespace=True,
                     names=["site", "year", "month", "value"])

    df["date"] = pd.to_datetime(df["year"].astype(str) + "-" +
                                df["month"].astype(str) + "-15")
    df = df.sort_values("date")
    df = df.dropna(subset=["value"])

    site = df["site"].iloc[0]

    # Extract seasonal component
    ts = df.set_index("date")["value"]
    seasonal = extract_seasonal_component(ts.values)
    df["seasonal"] = seasonal

    # Compute annual SCA
    df["year"] = df["date"].dt.year
    annual = df.groupby("year")["seasonal"].agg(["max", "min"])
    annual["SCA"] = annual["max"] - annual["min"]

    # Linear regression
    years = annual.index.values
    sca = annual["SCA"].values

    if len(years) < 2:
        print(f"Skipping {site}: not enough years for regression.")
        continue

    slope, intercept, r_value, p_value, stderr = linregress(years, sca)

    summary_rows.append([
        site,
        years.min(),
        years.max(),
        len(years),
        slope,
        p_value,
        r_value**2
    ])

    # Plot SCA vs year
    plt.figure(figsize=(8, 5))
    plt.scatter(years, sca, color="black", label="Annual SCA")
    plt.plot(years, intercept + slope * years, color="red",
             label=f"Trend = {slope:.3f} per yr\np = {p_value:.3g}")

    plt.title(f"{site} — SCA Trend")
    plt.xlabel("Year")
    plt.ylabel("Seasonal Cycle Amplitude (‰)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder, f"{site}_sca_trend.png"))
    plt.close()

# Save summary
summary_df = pd.DataFrame(summary_rows, columns=[
    "Site", "StartYear", "EndYear", "N_years",
    "Slope_per_year", "P_value", "R2"
])

summary_df.to_csv(os.path.join(folder, "sca_trend_summary.csv"), index=False)

print("Finished computing SCA trends for all sites.")
