#!/usr/bin/env python3
"""
SwanLab metric line-chart plotter.

Fetch scalar metrics via the SwanLab OOP Api, then render line charts
with matplotlib.  Each metric key gets its own subplot with trend line
and summary statistics (min / max / mean / std / last).

Usage:
    # Basic (requires `swanlab login` first)
    python plot_metrics.py username/project_name/run_id --keys loss,acc

    # Custom sample size and output path
    python plot_metrics.py username/project_name/run_id --keys loss --sample 500 -o loss_chart.png

    # Self-hosted server
    python plot_metrics.py username/project_name/run_id --keys loss --host https://your-swanserver

    # Explicit API key (skip login)
    python plot_metrics.py username/project_name/run_id --keys loss --api-key YOUR_KEY

    # Plot from a pre-fetched JSON file (skip API calls entirely)
    python plot_metrics.py --data metrics.json --keys loss,acc -o chart.png
"""

import argparse
import json
import math
import sys
from typing import Any, Dict, List, Optional, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from swanlab.api import Api

# ---------------------------------------------------------------------------
#  Data fetching
# ---------------------------------------------------------------------------


def fetch_metrics(
    path: str,
    keys: List[str],
    sample: int = 1500,
    api_key: Optional[str] = None,
    host: Optional[str] = None,
) -> Dict[str, Any]:
    """Fetch scalar metric data via the SwanLab Api."""
    api = Api(api_key=api_key, host=host)
    experiment = api.run(path)
    if not experiment.run_id:
        raise ValueError(f"Failed to fetch experiment at path '{path}'. Please verify the path and credentials.")
    return experiment.metrics(keys=keys, sample=sample, ignore_timestamp=True)


def fetch_summary(
    path: str,
    keys: List[str],
    api_key: Optional[str] = None,
    host: Optional[str] = None,
) -> Dict[str, Any]:
    """Fetch metric summary statistics via the SwanLab Api."""
    api = Api(api_key=api_key, host=host)
    experiment = api.run(path)
    if not experiment.run_id:
        raise ValueError(f"Failed to fetch experiment at path '{path}'. Please verify the path and credentials.")
    return experiment.summary(keys=keys)


# ---------------------------------------------------------------------------
#  Data extraction — pull (steps, values) pairs from the raw API dict
# ---------------------------------------------------------------------------


def extract_series(metric_data: Dict[str, Any], key: str) -> Tuple[List[int], List[float]]:
    """
    Extract (steps, values) for a single metric key.

    Handles two common response shapes:
    - ``{key: [{"step": s, "value": v}, ...]}``          (flat list)
    - ``{key: {"data": [{"step": s, "value": v}, ...]}}`` (nested under "data")
    """
    entry = metric_data.get(key)
    if entry is None:
        return [], []

    if isinstance(entry, dict) and "data" in entry:
        points = entry["data"]
    elif isinstance(entry, list):
        points = entry
    else:
        return [], []

    steps: List[int] = []
    values: List[float] = []
    for pt in points:
        if not isinstance(pt, dict):
            continue
        step = pt.get("step")
        value = pt.get("value")
        if step is not None and value is not None:
            steps.append(int(step))
            values.append(float(value))
    return steps, values


# ---------------------------------------------------------------------------
#  Statistics helpers
# ---------------------------------------------------------------------------


def compute_stats(values: List[float]) -> Dict[str, float]:
    """Compute basic descriptive statistics for a value series."""
    if not values:
        return {}
    arr = np.array(values)
    return {
        "min": float(np.min(arr)),
        "max": float(np.max(arr)),
        "mean": float(np.mean(arr)),
        "std": float(np.std(arr)),
        "last": float(arr[-1]),
    }


def format_stats_text(stats: Dict[str, float]) -> str:
    """Format statistics dict into a compact annotation string."""
    if not stats:
        return "N/A"
    parts = [
        f"min={stats['min']:.6g}",
        f"max={stats['max']:.6g}",
        f"mean={stats['mean']:.6g}",
        f"std={stats['std']:.6g}",
        f"last={stats['last']:.6g}",
    ]
    return "  ".join(parts)


# ---------------------------------------------------------------------------
#  Plotting
# ---------------------------------------------------------------------------


def plot_line_chart(
    metric_data: Dict[str, Any],
    keys: List[str],
    output: str = "metrics_chart.png",
    title: Optional[str] = None,
    dpi: int = 150,
) -> None:
    """
    Render one subplot per metric key with a line chart and stats annotation.

    :param metric_data: Raw dict returned by ``Api.metrics()``.
    :param keys: Metric keys to plot.
    :param output: Output image file path.
    :param title: Overall chart title (auto-generated if None).
    :param dpi: Image resolution.
    """
    n = len(keys)
    if n == 0:
        print("No keys to plot.", file=sys.stderr)
        return

    # Layout: at most 2 subplots per row
    ncols = min(n, 2)
    nrows = math.ceil(n / ncols)

    fig, axes = plt.subplots(nrows, ncols, figsize=(7 * ncols, 4.5 * nrows), squeeze=False)
    fig.suptitle(title or "SwanLab Metrics", fontsize=14, fontweight="bold", y=0.98)

    for idx, key in enumerate(keys):
        row, col = divmod(idx, ncols)
        ax = axes[row][col]

        steps, values = extract_series(metric_data, key)

        if not steps:
            ax.set_title(key, fontsize=11)
            ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes, fontsize=12, color="gray")
            ax.set_xlabel("step")
            ax.set_ylabel(key)
            continue

        # Line chart
        ax.plot(steps, values, linewidth=1.0, alpha=0.85)

        # Highlight the last data point
        ax.scatter([steps[-1]], [values[-1]], s=20, zorder=5)

        # Stats annotation box
        stats = compute_stats(values)
        stats_text = format_stats_text(stats)
        ax.text(
            0.02,
            0.02,
            stats_text,
            transform=ax.transAxes,
            fontsize=7,
            verticalalignment="bottom",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.8),
        )

        ax.set_title(key, fontsize=11)
        ax.set_xlabel("step")
        ax.set_ylabel(key)
        ax.grid(True, alpha=0.3)

    # Hide unused subplots
    for idx in range(n, nrows * ncols):
        row, col = divmod(idx, ncols)
        axes[row][col].set_visible(False)

    plt.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Chart saved to: {output}")


# ---------------------------------------------------------------------------
#  CLI entry point
# ---------------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fetch SwanLab experiment metrics and plot line charts.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "path",
        nargs="?",
        default=None,
        help="Experiment path: username/project_name/run_id (optional when --data is used)",
    )
    parser.add_argument(
        "--data", default=None, help="Path to a JSON file with pre-fetched metric data. Skips API calls entirely."
    )
    parser.add_argument("--keys", "-k", required=True, help="Comma-separated metric keys, e.g. 'loss,acc'")
    parser.add_argument("--sample", "-s", type=int, default=1500, help="Sample size (default: 1500)")
    parser.add_argument(
        "--output", "-o", default="metrics_chart.png", help="Output image path (default: metrics_chart.png)"
    )
    parser.add_argument("--title", "-t", default=None, help="Chart title (default: auto)")
    parser.add_argument("--dpi", type=int, default=150, help="Image DPI (default: 150)")
    parser.add_argument("--api-key", default=None, help="SwanLab API key (or use swanlab login)")
    parser.add_argument("--host", default=None, help="SwanLab API host URL (for self-hosted)")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    keys = [k.strip() for k in args.keys.split(",") if k.strip()]

    if not keys:
        print("Error: --keys must contain at least one key.", file=sys.stderr)
        sys.exit(1)

    if not args.path and not args.data:
        print("Error: provide either a PATH argument or --data <json_file>.", file=sys.stderr)
        sys.exit(1)

    if args.data:
        # Load metric data from a pre-fetched JSON file
        with open(args.data, "r", encoding="utf-8") as f:
            metric_data = json.load(f)
        print(f"Loaded metric data from: {args.data}")
    else:
        print(f"Fetching metrics for: {args.path}")
        print(f"Keys: {keys}, sample: {args.sample}")
        metric_data = fetch_metrics(
            path=args.path,
            keys=keys,
            sample=args.sample,
            api_key=args.api_key,
            host=args.host,
        )

    plot_line_chart(
        metric_data=metric_data,
        keys=keys,
        output=args.output,
        title=args.title or f"Metrics: {args.path}",
        dpi=args.dpi,
    )


if __name__ == "__main__":
    main()
