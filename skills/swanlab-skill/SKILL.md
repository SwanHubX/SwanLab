---
name: swanlab-skill
metadata:
  version: "0.1.0"
description: >
  Interact with SwanLab — both writing tracking code (init/log/finish/multimedia) and querying
  experiment data via CLI (`swanlab api`). Use this skill when the user wants to write training
  tracking code, log metrics or media, manage experiments, inspect metrics/logs/columns, list
  projects/runs, filter experiments, manage self-hosted users, automate queries via CLI, or
  mentions "swanlab", "experiment tracking", "log metrics", "swanlab api", "swanlab cli".
---

# SwanLab Skill

SwanLab is an AI training experiment tracking platform. This skill covers two usage patterns:

- **Writing tracking code** — use the Python SDK (`swanlab.init`, `swanlab.log`, `swanlab.finish`, media helpers)
- **Reading experiment data** — use the `swanlab api` CLI to query metrics, logs, summaries, media, etc.

---

## Reference Routing

| If the user wants to... | Read this reference |
|---|---|
| Write tracking code (init/log/finish/media) | `references/SDK_QUICKSTART.md` |
| Query data via CLI (metrics/summary/logs/filter/etc.) | `references/CLI_REFERENCE.md` |
| Understand data model / terminology / filter syntax | `references/SWANLAB_CONCEPTS.md` |
| Plot metrics or compare experiments visually | See **Scripts** below |

---

## Run Modes

`swanlab.init(mode=...)` controls where data goes:

| Mode | Local Storage | Cloud Upload | Use Case |
|------|--------------|-------------|----------|
| `online` | Yes (protobuf) | Yes (Transport → HTTP) | Normal cloud usage. Requires login. |
| `local` | Yes (protobuf) | No | Air-gapped / no account needed. |
| `offline` | Yes (protobuf) | No (syncable later via `swanlab sync`) | Save locally, upload to cloud later. |
| `disabled` | No | No | Completely disable all logging. |

Default is `online` if logged in, otherwise the user is prompted interactively (or falls back to `offline`).

---

## Scripts

Two helper scripts are available for visualizing experiment data:

### `scripts/plot_metrics.py` — Single Experiment Line Chart

Trigger when the user wants to **visualize scalar metrics from one experiment** (e.g. "plot my loss curve", "show training metrics chart").

```bash
python scripts/plot_metrics.py username/project_name/run_id --keys loss,acc
python scripts/plot_metrics.py user/proj/run1 -k loss -o loss_chart.png -s 500
python scripts/plot_metrics.py --data metrics.json -k loss,acc -o chart.png   # from saved JSON
```

### `scripts/runs_benchmark.py` — Cross-Experiment Comparison

Trigger when the user wants to **compare the same metric across multiple experiments** (e.g. "compare loss across runs", "benchmark these experiments").

```bash
python scripts/runs_benchmark.py user/proj/run1 user/proj/run2 user/proj/run3 -k loss
python scripts/runs_benchmark.py user/proj/run1 user/proj/run2 -k loss --direction lower
python scripts/runs_benchmark.py user/proj/run1 user/proj/run2 -k loss,acc --normalize
python scripts/runs_benchmark.py --data benchmark_data.json -k loss              # from saved JSON
```

Both scripts require `swanlab login` (or `--api-key` / `--host` flags).

---

## Path Convention

Several CLI commands take a `PATH` argument:

- **Project path**: `username/project_name`
- **Experiment path**: `username/project_name/run_id`

---

## Quick Disambiguation

| User says... | They probably mean... | Route |
|---|---|---|
| "track my training" / "log metrics" | Write tracking code | `SDK_QUICKSTART.md` |
| "log images/audio/text" | Log media data | `SDK_QUICKSTART.md` |
| "my loss curve" / "experiment metrics" | Query scalar data | `CLI_REFERENCE.md > run metrics` |
| "filter experiments" | Query by conditions | `CLI_REFERENCE.md > run filter` |
| "my experiments" / "list runs" | List experiments | `CLI_REFERENCE.md > run list` |
| "compare runs visually" | Cross-experiment chart | `scripts/runs_benchmark.py` |
| "plot metric chart" | Single-experiment chart | `scripts/plot_metrics.py` |
| "experiment config" | Hyperparameters | `CLI_REFERENCE.md > run info` |
| "console output" | Captured logs | `CLI_REFERENCE.md > run logs` |
| "what columns are tracked" | Metric definitions | `CLI_REFERENCE.md > run columns` |
| "check connectivity" / "can I reach swanlab" | Environment check | `swanlab ping` |

---

## Environment Connectivity

Before writing tracking code or running CLI queries, especially in `online` mode, use `swanlab ping` to verify that the current environment can reach the SwanLab server:

```bash
swanlab ping
# Reports: API host, web host, latency, and login status.
# If ping fails, check SWANLAB_API_HOST / network proxy / firewall settings.
```

This is the fastest way to diagnose connectivity issues — run it first when a user reports upload failures, login problems, or unknown mode fallbacks.

---

## Behavioral Constraints

- **Never use `--all` unless the user explicitly asks for it.** This bypasses pagination and fetches everything. Always use default pagination first.
- **Always ask for specific column keys before running `run metrics`, `run medias`, or `run column`.** If unknown, run `run columns PATH` first to discover keys.
- **Always add `--ignore-timestamp` for `run metrics` and `run logs`** unless the user specifically needs timestamps.
- All CLI output is JSON to stdout. Pipe to `jq` for further processing.
