from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def summarize_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    years = df.index.to_numpy()
    a = df.to_numpy(dtype=float)

    neg_pct = np.mean(a < 0, axis=1) * 100.0
    nan_pct = np.mean(~np.isfinite(a), axis=1) * 100.0

    # Use nan-safe reductions
    row_min = np.nanmin(a, axis=1)
    row_p01 = np.nanpercentile(a, 1, axis=1)
    row_p50 = np.nanpercentile(a, 50, axis=1)
    row_p99 = np.nanpercentile(a, 99, axis=1)
    row_max = np.nanmax(a, axis=1)

    out = pd.DataFrame(
        {
            "Year": years,
            "neg_pct": neg_pct,
            "nan_pct": nan_pct,
            "min": row_min,
            "p01": row_p01,
            "p50": row_p50,
            "p99": row_p99,
            "max": row_max,
        }
    )
    return out


def first_issue(summary: pd.DataFrame, *, neg_threshold_pct: float = 0.1, max_threshold: float = 1e3) -> str:
    # “non-physical” heuristics:
    # - any meaningful negative fraction
    # - extremely large values (default 1000 Tg/yr)
    neg_mask = summary["neg_pct"] > neg_threshold_pct
    big_mask = summary["max"] > max_threshold

    if not (neg_mask.any() or big_mask.any()):
        return "No negatives or extreme values detected by current thresholds."

    first_idx = None
    if neg_mask.any():
        first_idx = int(np.argmax(neg_mask.to_numpy()))
        reason = f"negatives appear (neg_pct={summary.loc[first_idx,'neg_pct']:.2f}%)"
    if big_mask.any():
        first_big = int(np.argmax(big_mask.to_numpy()))
        if first_idx is None or first_big < first_idx:
            first_idx = first_big
            reason = f"extreme max appears (max={summary.loc[first_idx,'max']:.2f})"

    row = summary.loc[first_idx]
    return (
        f"First flagged year: {int(row['Year'])} ({reason}). "
        f"min={row['min']:.2f}, p01={row['p01']:.2f}, p50={row['p50']:.2f}, p99={row['p99']:.2f}, max={row['max']:.2f}."
    )


def main() -> None:
    out_dir = Path("Output")
    files = {
        "BB": out_dir / "BB_3source_MC_alliterations.csv",
        "FF": out_dir / "FF_3source_MC_alliterations.csv",
        "Mic": out_dir / "Mic_3source_MC_alliterations.csv",
    }

    for name, path in files.items():
        if not path.exists():
            raise SystemExit(f"Missing: {path}")

        summary = summarize_matrix(path)
        print(f"\n=== {name} ===")
        print(first_issue(summary))

        overall_neg = float((summary["neg_pct"] > 0).mean() * 100.0)
        print(f"Years with any negatives: {overall_neg:.1f}% ({int((summary['neg_pct']>0).sum())}/{len(summary)})")


if __name__ == "__main__":
    main()

