import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

# ============================================================
# Helper: SSA decomposition
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
# Step 1 — Scan folder and summarize sites
# ============================================================
folder = os.path.dirname(os.path.abspath(__file__))

summary = []
all_files = [f for f in os.listdir(folder) if f.endswith(".txt")]

for file in all_files:
    path = os.path.join(folder, file)
    df = pd.read_csv(path, comment="#", delim_whitespace=True,
                     names=["site", "year", "month", "value"])

    df["date"] = pd.to_datetime(df["year"].astype(str) + "-" +
                                df["month"].astype(str) + "-15")

    df = df.sort_values("date")
    df_valid = df.dropna(subset=["value"])

    site = df["site"].iloc[0]
    start = df_valid["date"].min().date()
    end = df_valid["date"].max().date()
    n = len(df_valid)

    summary.append([site, start, end, n])

summary_df = pd.DataFrame(summary, columns=["Site", "Start", "End", "N_points"])
print(summary_df)

# ============================================================
# Step 2 — Plot seasonal cycle using ALL data
# ============================================================
output_folder = os.path.join(folder, "plots_all_data")
os.makedirs(output_folder, exist_ok=True)

for file in all_files:
    path = os.path.join(folder, file)
    df = pd.read_csv(path, comment="#", delim_whitespace=True,
                     names=["site", "year", "month", "value"])

    df["date"] = pd.to_datetime(df["year"].astype(str) + "-" +
                                df["month"].astype(str) + "-15")
    df = df.sort_values("date")
    df = df.dropna(subset=["value"])

    site = df["site"].iloc[0]

    # Monthly time series
    ts = df.set_index("date")["value"]

    # Extract seasonal component
    seasonal = extract_seasonal_component(ts.values)
    df["seasonal"] = seasonal

    # Compute mean seasonal cycle across ALL years
    df["month"] = df["date"].dt.month
    mean_cycle = df.groupby("month")["seasonal"].mean()
    std_cycle = df.groupby("month")["seasonal"].std()

    # Plot
    fig, ax = plt.subplots(figsize=(10, 5))

    months = np.arange(1, 13)

    ax.plot(months, mean_cycle, color="blue", label="Mean seasonal cycle")
    ax.fill_between(months,
                    mean_cycle - std_cycle,
                    mean_cycle + std_cycle,
                    color="blue", alpha=0.2,
                    label="±1 std (all years)")

    ax.set_title(f"{site} δ13C–CH₄ Seasonal Cycle (All Data)")
    ax.set_xlabel("Month")
    ax.set_ylabel("δ13C–CH₄ (‰)")
    ax.set_xticks(months)
    ax.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(output_folder, f"{site}_seasonal_cycle_all_data.png"))
    plt.close()

print("All plots saved using full dataset for each site.")
