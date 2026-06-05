#!/usr/bin/env python3
"""
Cross-experiment benchmark plotter (runs_benchmark.py).

Fetch the same metric from multiple SwanLab experiments via the OOP Api,
overlay line charts (raw trend + running-best step line) for visual
comparison, and print a summary statistics table.

Usage:
    # Compare one metric across 3 experiments (requires `swanlab login`)
    python runs_benchmark.py user/proj/run1 user/proj/run2 user/proj/run3 -k loss

    # Lower-is-better metric with custom labels
    python runs_benchmark.py user/proj/run1 user/proj/run2 -k latency_ms \\
           --direction lower --labels "v1,v2"

    # Multiple keys (one subplot per key)
    python runs_benchmark.py user/proj/run1 user/proj/run2 -k loss,acc

    # Normalize x-axis to 0-100% for unequal-length experiments
    python runs_benchmark.py user/proj/run1 user/proj/run2 -k loss --normalize

    # Self-hosted server
    python runs_benchmark.py user/proj/run1 user/proj/run2 -k loss \\
           --host https://your-swanserver --api-key YOUR_KEY

    # Benchmark from a pre-fetched JSON file (skip API calls entirely)
    python runs_benchmark.py --data benchmark_data.json -k loss
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from swanlab.api import Api

# ---------------------------------------------------------------------------
#  Data model
# ---------------------------------------------------------------------------


@dataclass
class ExperimentSeries:
    """Metric time-series for a single experiment."""

    name: str
    path: str
    steps: List[int]
    values: List[float]


# ---------------------------------------------------------------------------
#  Data fetching — one Api call per experiment, no repeated requests
# ---------------------------------------------------------------------------


def _create_api(api_key: Optional[str] = None, host: Optional[str] = None) -> Api:
    return Api(api_key=api_key, host=host)


def fetch_series(
    path: str,
    key: str,
    sample: int = 1500,
    api_key: Optional[str] = None,
    host: Optional[str] = None,
) -> ExperimentSeries:
    """Fetch a single metric key from one experiment."""
    api = _create_api(api_key, host)
    experiment = api.run(path)
    name = experiment.name
    raw = experiment.metrics(keys=[key], sample=sample, ignore_timestamp=True)
    steps, values = _extract(raw, key)
    return ExperimentSeries(name=name, path=path, steps=steps, values=values)


def _extract(metric_data: Dict[str, Any], key: str) -> Tuple[List[int], List[float]]:
    """Pull (steps, values) from the raw API response dict."""
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
    """Return {min, max, mean, std, last} for a value series."""
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


def running_best(values: List[float], direction: str) -> List[float]:
    """
    Cumulative best-so-far.

    direction="higher" → running max (higher is better).
    direction="lower"  → running min (lower is better).
    """
    if not values:
        return []
    best = float("-inf") if direction == "higher" else float("inf")
    result: List[float] = []
    for v in values:
        if direction == "higher":
            best = max(best, v)
        else:
            best = min(best, v)
        result.append(best)
    return result


def normalize_steps(steps: List[int]) -> List[float]:
    """Map steps to [0, 100] percentage scale."""
    if not steps:
        return []
    total = max(steps) if max(steps) != 0 else 1
    return [s / total * 100.0 for s in steps]


# ---------------------------------------------------------------------------
#  Summary table
# ---------------------------------------------------------------------------


def print_comparison_table(
    series_list: List[ExperimentSeries],
    key: str,
) -> None:
    """Print a formatted comparison table to stdout."""
    header = f"{'Experiment':<24s} | {'Steps':>5s} | {'Min':>10s} | {'Max':>10s} | {'Mean':>10s} | {'Std':>10s} | {'Last':>10s}"
    sep = "-" * len(header)
    print(f"\n  Benchmark: {key}")
    print(f"  {sep}")
    print(f"  {header}")
    print(f"  {sep}")
    for s in series_list:
        stats = compute_stats(s.values)
        if not stats:
            print(
                f"  {s.name:<24s} | {'N/A':>5s} | {'N/A':>10s} | {'N/A':>10s} | {'N/A':>10s} | {'N/A':>10s} | {'N/A':>10s}"
            )
            continue
        print(
            f"  {s.name:<24s} | {len(s.steps):>5d} | {stats['min']:>10.6g} | {stats['max']:>10.6g} "
            f"| {stats['mean']:>10.6g} | {stats['std']:>10.6g} | {stats['last']:>10.6g}"
        )
    print(f"  {sep}\n")


# ---------------------------------------------------------------------------
#  Plotting
# ---------------------------------------------------------------------------

# Default color cycle
_COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]


def _pick_color(idx: int) -> str:
    return _COLORS[idx % len(_COLORS)]


def plot_benchmark(
    all_series: Dict[str, List[ExperimentSeries]],
    keys: List[str],
    output: str,
    title: Optional[str],
    direction: str,
    normalize: bool,
    dpi: int,
) -> None:
    """
    Render benchmark chart — one subplot per key, overlays all experiments.

    Each subplot draws:
      - A thin line for the raw metric trend.
      - A thicker step line for the running-best envelope.
    """
    n_keys = len(keys)
    if n_keys == 0:
        print("No keys to plot.", file=sys.stderr)
        return

    ncols = min(n_keys, 2)
    nrows = math.ceil(n_keys / ncols)

    fig, axes = plt.subplots(nrows, ncols, figsize=(9 * ncols, 5.5 * nrows), squeeze=False)
    fig.suptitle(
        title or "Cross-Experiment Benchmark",
        fontsize=14,
        fontweight="bold",
        y=0.98,
    )

    for k_idx, key in enumerate(keys):
        row_idx, col_idx = divmod(k_idx, ncols)
        ax = axes[row_idx][col_idx]
        series_list = all_series.get(key, [])

        if not series_list or all(len(s.values) == 0 for s in series_list):
            ax.set_title(key, fontsize=12)
            ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes, fontsize=13, color="gray")
            ax.set_xlabel("step" if not normalize else "progress (%)")
            ax.set_ylabel(key)
            continue

        for exp_idx, s in enumerate(series_list):
            if not s.steps:
                continue

            x = normalize_steps(s.steps) if normalize else s.steps
            color = _pick_color(exp_idx)
            label_base = s.name

            # Raw trend line
            ax.plot(
                x,
                s.values,
                linewidth=1.0,
                alpha=0.6,
                color=color,
                label=label_base,
            )

            # Running-best step line
            rb = running_best(s.values, direction)
            ax.step(
                x,
                rb,
                where="post",
                linewidth=2.2,
                alpha=0.85,
                color=color,
                linestyle="--",
                label=f"{label_base} (best)",
            )

            # Last-point marker
            ax.scatter([x[-1]], [s.values[-1]], s=24, color=color, zorder=5, edgecolors="white", linewidths=0.6)

        ax.set_title(key, fontsize=12)
        ax.set_xlabel("step" if not normalize else "progress (%)", fontsize=11)
        ax.set_ylabel(key, fontsize=11)
        ax.grid(True, alpha=0.3)

        # De-duplicate legend entries: keep only one per experiment (raw line),
        # add a single generic "best" proxy for the step lines.
        handles, labels = ax.get_legend_handles_labels()
        seen: set[str] = set()
        unique_handles = []
        unique_labels = []
        best_proxy = None
        for h, lbl in zip(handles, labels):
            if lbl.endswith(" (best)"):
                if best_proxy is None:
                    best_proxy = h
                continue
            if lbl not in seen:
                seen.add(lbl)
                unique_handles.append(h)
                unique_labels.append(lbl)

        # Append a single "running best" entry
        if best_proxy is not None:
            unique_handles.append(best_proxy)
            unique_labels.append("running best (--)")
        ax.legend(unique_handles, unique_labels, loc="best", frameon=True, fontsize=8)

        # Stats annotation box
        lines: List[str] = []
        for s in series_list:
            stats = compute_stats(s.values)
            if stats:
                lines.append(f"{s.name}: last={stats['last']:.4g}  mean={stats['mean']:.4g}")
        if lines:
            ax.text(
                0.02,
                0.02,
                "\n".join(lines),
                transform=ax.transAxes,
                fontsize=7,
                verticalalignment="bottom",
                fontfamily="monospace",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.8),
            )

    # Hide unused subplots
    for idx in range(n_keys, nrows * ncols):
        r, c = divmod(idx, ncols)
        axes[r][c].set_visible(False)

    plt.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Chart saved to: {output}")


# ---------------------------------------------------------------------------
#  CLI
# ---------------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Cross-experiment benchmark: compare the same metric across multiple runs.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "paths",
        nargs="*",
        help="Experiment paths (username/project/run_id). Optional when --data is used.",
    )
    parser.add_argument(
        "--keys",
        "-k",
        required=True,
        help="Comma-separated metric keys to compare, e.g. 'loss,acc'.",
    )
    parser.add_argument(
        "--labels",
        default=None,
        help="Comma-separated display names for each experiment (matched by order).",
    )
    parser.add_argument(
        "--sample",
        "-s",
        type=int,
        default=1500,
        help="Sample size per metric (default: 1500).",
    )
    parser.add_argument(
        "--direction",
        choices=["higher", "lower"],
        default="higher",
        help="Direction for 'running best': 'higher' (max) or 'lower' (min). Default: higher.",
    )
    parser.add_argument(
        "--normalize",
        action="store_true",
        help="Normalize x-axis to 0-100%% for unequal-length experiments.",
    )
    parser.add_argument(
        "--output",
        "-o",
        default="benchmark.png",
        help="Output image path (default: benchmark.png).",
    )
    parser.add_argument("--title", "-t", default=None, help="Chart title (default: auto).")
    parser.add_argument("--dpi", type=int, default=150, help="Image DPI (default: 150).")
    parser.add_argument("--api-key", default=None, help="SwanLab API key (or use swanlab login).")
    parser.add_argument(
        "--data", default=None, help="Path to a JSON file with pre-fetched benchmark data. Skips API calls."
    )
    parser.add_argument("--host", default=None, help="SwanLab API host URL (for self-hosted).")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    keys = [k.strip() for k in args.keys.split(",") if k.strip()]

    if not keys:
        print("Error: --keys must contain at least one key.", file=sys.stderr)
        return 1

    if not args.paths and not args.data:
        print("Error: provide either PATH arguments or --data <json_file>.", file=sys.stderr)
        return 1

    all_series: Dict[str, List[ExperimentSeries]] = {k: [] for k in keys}

    if args.data:
        # Load pre-fetched benchmark data from JSON
        with open(args.data, "r", encoding="utf-8") as f:
            raw_data = json.load(f)
        # Expected structure: {"experiments": [{"name": "...", "path": "...", "metrics": {<key>: [{"step": s, "value": v}, ...]}}]}
        experiments = raw_data if isinstance(raw_data, list) else raw_data.get("experiments", [raw_data])
        for exp in experiments:
            name = exp.get("name", "unknown")
            path = exp.get("path", "")
            metrics = exp.get("metrics", exp)  # fallback: top-level dict if flat
            for key in keys:
                steps, values = _extract(metrics, key)
                all_series[key].append(ExperimentSeries(name=name, path=path, steps=steps, values=values))
        print(f"Loaded benchmark data from: {args.data}")
    else:
        # Fetch from API
        custom_labels: List[str] = []
        if args.labels:
            custom_labels = [lbl.strip() for lbl in args.labels.split(",")]
            if len(custom_labels) != len(args.paths):
                print(
                    f"Error: --labels has {len(custom_labels)} items but {len(args.paths)} paths given.",
                    file=sys.stderr,
                )
                return 1

        api = _create_api(args.api_key, args.host)

        for i, path in enumerate(args.paths):
            print(f"Fetching [{i + 1}/{len(args.paths)}] {path} ...")
            experiment = api.run(path)
            if not experiment.run_id:
                print(
                    f"Error: Failed to fetch experiment at path '{path}'. Please verify the path and credentials.",
                    file=sys.stderr,
                )
                return 1
            name = custom_labels[i] if custom_labels else experiment.name
            raw = experiment.metrics(keys=keys, sample=args.sample, ignore_timestamp=True)

            for key in keys:
                steps, values = _extract(raw, key)
                all_series[key].append(ExperimentSeries(name=name, path=path, steps=steps, values=values))

    # Print comparison tables
    for key in keys:
        print_comparison_table(all_series[key], key)

    # Plot
    plot_benchmark(
        all_series=all_series,
        keys=keys,
        output=args.output,
        title=args.title,
        direction=args.direction,
        normalize=args.normalize,
        dpi=args.dpi,
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
