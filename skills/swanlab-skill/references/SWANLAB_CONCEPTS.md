# SwanLab Concepts & Data Model

Reference for understanding SwanLab's entity hierarchy, terminology, and data model. Read this when interpreting user requests or deciding which CLI command to use.

---

## Product Overview

SwanLab is an AI training experiment tracking platform. The Python SDK records metrics, configs, media, console output, and hardware stats during training. Data is persisted locally and optionally synced to the SwanLab cloud (or self-hosted instance).

The CLI (`swanlab api`) provides read-only query access to all tracked data from the command line.

---

## Core Entity Hierarchy

```
User (authenticated account)
+-- Workspace (personal or team namespace, identified by username)
|   +-- Project (groups related experiments)
|   |   +-- Experiment / Run (a single training execution)
|   |   |   +-- Config (input hyperparameters, set at init, immutable)
|   |   |   +-- Column (a metric definition — e.g. "loss", "acc", "image")
|   |   |   |   +-- Scalar Metrics (time-series numeric data)
|   |   |   |   +-- Media Metrics (images, audio, video, molecules, etc.)
|   |   |   +-- Console Logs (captured stdout/stderr output)
|   |   |   +-- Profile (metadata: config, requirements, conda, hardware)
```

### Entity Reference

| Entity | Identified by | Description |
|--------|--------------|-------------|
| **User** | implicit (logged-in) | The authenticated account. Has `username`, `name`, `bio`, `email`, etc. |
| **Workspace** | `username` | A personal or team namespace. Type is `PERSON` or `TEAM`. Every project belongs to a workspace. |
| **Project** | `username/project_name` (path) | Groups related experiments. Has visibility (`PUBLIC` / `PRIVATE`), description, labels. |
| **Experiment (Run)** | `username/project_name/run_id` (path) | The atomic unit — one execution of training code. Created by `swanlab.init()`. Has state, config, metrics, logs. |
| **Column** | `key` within an experiment | A metric definition. Belongs to either `CUSTOM` (user-defined) or `SYSTEM` (auto-collected). Has a data type. |
| **Metric** | queried via column key | The actual time-series data for a column. Three categories: Scalar, Media, Log. |

---

## Key Terminology

### Experiment States

An experiment (run) is always in one of these states:

| State | Meaning |
|-------|---------|
| `RUNNING` | Currently executing |
| `FINISHED` | Completed successfully |
| `CRASHED` | Terminated due to an error |
| `ABORTED` | Manually stopped by the user |
| `OFFLINE` | Created in offline mode, not yet synced |

### Column Classes

| Class | Description |
|-------|-------------|
| `CUSTOM` | User-defined metrics logged via `swanlab.log()` (default) |
| `SYSTEM` | Auto-collected system metrics (CPU, GPU, memory, etc.) |

### Column Data Types

| Type | Category | Description |
|------|----------|-------------|
| `FLOAT` | Scalar | Floating-point numeric values |
| `BOOLEAN` | Scalar | True/false values |
| `STRING` | Scalar | Text values |
| `IMAGE` | Media | Image files (PNG, JPG, etc.) |
| `AUDIO` | Media | Audio files |
| `VIDEO` | Media | Video files |
| `OBJECT3D` | Media | 3D point cloud data (JSON) |
| `MOLECULE` | Media | Biochemical molecule structures |
| `ECHARTS` | Media | Custom chart configurations (JS/TS) |
| `TABLE` | Media | Tabular data |
| `TEXT` | Media | Rich text content |

### Project Visibility

| Visibility | Description |
|------------|-------------|
| `PRIVATE` | Only accessible to workspace members (default) |
| `PUBLIC` | Visible to anyone |

---

## Path Convention

The SwanLab API uses a hierarchical path format to address entities:

```
username                        → Workspace
username/project_name           → Project
username/project_name/run_id    → Experiment
```

- `username`: The workspace's unique identifier.
- `project_name`: Human-readable project name (1-100 chars, `0-9a-zA-Z-_.+`).
- `run_id`: A unique identifier (CUID) assigned to each experiment.

---

## Three Metric Categories

Understanding which metric type to query is critical for answering "where's my data?":

### Scalar Metrics (e.g. loss, accuracy)

Numeric time-series data logged at each training step.

```
Structure per data point:
{
  "index": 0,        // step number or custom x-axis value
  "data": 0.523,     // the metric value (can be int, float, or "NaN"/"INF"/"-INF")
  "timestamp": 1714368000  // Unix timestamp
}
```

**Statistics** are computed server-side per scalar column: `min`, `max`, `avg`, `median`, `latest`.

**Sampling**: By default, `--sample 1500` returns downsampled data for visualization. Use `--all` to get a CSV download URL for the full dataset.

### Media Metrics (e.g. images, audio)

File-based data stored per step. Each media item has a presigned download URL.

```
Structure:
{
  "steps": [0, 1, 2],           // available steps
  "step": 0,                     // current step
  "metrics": [
    {
      "index": 0,
      "items": [
        { "url": "https://..." }  // presigned download URL
      ]
    }
  ]
}
```

**Step-based**: Media is retrieved one step at a time by default (`--step N`). Use `--all` to fetch all steps.

### Console Logs

Captured stdout/stderr output from the training process.

```
Structure per log entry:
{
  "epoch": 0,            // shard index (logs are paginated by shards)
  "level": "INFO",       // DEBUG | INFO | WARN | ERROR
  "message": "Epoch 1/10",
  "tag": "",             // optional tag
  "timestamp": "..."     // ISO timestamp
}
```

**Shard-based pagination**: Logs are fetched by `--offset` (shard index), not page number. Increase offset to get later shards.

---

## Experiment Profile

Each experiment carries a `profile` object containing metadata about the run:

```
{
  "config": { ... },           // Hyperparameters set at swanlab.init(config={...})
  "metadata": { ... },         // Hardware & environment info (auto-collected)
  "requirements": "...",       // pip requirements string
  "conda": "..."               // conda environment YAML
}
```

- **Config** = user inputs. Set once at `swanlab.init()`. Does not change during training.
- **Metadata** = auto-collected system info (Python version, GPU model, OS, etc.).
- **Requirements / Conda** = captured from the active Python environment.

---

## Self-Hosted Instance

SwanLab can be deployed as a self-hosted instance (private deployment). Self-hosted commands are **only valid when the target host is NOT `swanlab.cn`**.

**Host detection rule**: Before using any self-hosted command, check where requests will be sent:
1. If the user passes `--host` explicitly → use that value.
2. Otherwise, check `SWANLAB_API_HOST` / `SWANLAB_WEB_HOST` environment variables.
3. Otherwise, check `.netrc` or SwanLab config (`~/.swanlab`).

If the resolved host contains `swanlab.cn`, self-hosted commands will fail and should not be attempted. Only proceed when the host points to a self-hosted deployment.

Self-hosted-specific features:

| Property | Description |
|----------|-------------|
| `enabled` | Whether this is a self-hosted deployment |
| `expired` | Whether the license has expired |
| `root` | Whether the current user is a root (admin) |
| `plan` | License plan: `free` or `commercial` |
| `seats` | Number of allowed user seats |

Management operations (create user, list users) require root privileges.

---

## Quick Disambiguation

| User says... | They probably mean... | CLI command |
|---|---|---|
| "my experiments" | Experiments in a project | `swanlab api run list -p user/project` |
| "training metrics" | Scalar time-series data | `swanlab api run metrics PATH --keys loss,acc` |
| "my loss curve" | Scalar data for the "loss" key | `swanlab api run metrics PATH --keys loss` |
| "experiment config" | Hyperparameters set at init | `swanlab api run info PATH` (look at `profile.config`) |
| "logged images" | Media metrics | `swanlab api run medias PATH --keys image` |
| "console output" | Captured logs | `swanlab api run logs PATH` |
| "what columns are tracked" | Metric definitions | `swanlab api run columns PATH` |
| "my projects" | Projects in a workspace | `swanlab api project list` |
| "create a new project" | Make a project | `swanlab api project create -n NAME` |
| "who am I" | Current user info | `swanlab api user info` |
| "self-hosted users" | User management | `swanlab api selfhosted list-users` |
| "best loss" | Minimum scalar value | `swanlab api run metrics PATH --keys loss` (check `min` field) |

---

## Key Environment Variables

| Variable | Purpose |
|----------|---------|
| `SWANLAB_API_KEY` | API authentication key |
| `SWANLAB_API_HOST` | API server URL (default `https://api.swanlab.cn`) |
| `SWANLAB_WEB_HOST` | Web dashboard URL (default `https://swanlab.cn`) |
| `SWANLAB_LOG_DIR` | Custom local log directory (default `./swanlog`) |
| `SWANLAB_SAVE_DIR` | Root dir for SwanLab files (default `~/.swanlab`) |

---

## Run Modes

| Mode | Local Storage | Cloud Upload | Use Case |
|------|--------------|-------------|----------|
| `online` | Yes | Yes | Normal cloud usage |
| `local` | Yes | No | Air-gapped / no account |
| `offline` | Yes | No (syncable later via `swanlab sync`) | Cloud later |
| `disabled` | No | No | Disable all logging |
