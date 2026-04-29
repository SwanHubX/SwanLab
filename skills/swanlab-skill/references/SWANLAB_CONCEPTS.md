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

Management operations (create user, list users, list projects, list workspaces, summary) all require root privileges. If a non-root user attempts any of these, the API returns a permission error.

### Instance-Level Management

On a self-hosted deployment, the root (admin) user can query across the entire instance, not just within a single workspace:

| Operation | Root Required | Description |
|-----------|---------------|-------------|
| Create User | Yes | Add a new user to the instance |
| List Users | Yes | All registered users across the instance |
| List Projects | Yes | All projects across all workspaces, with optional filters for creator, workspace, and keyword |
| List Workspaces | Yes | All workspaces (personal + team) on the instance, with optional keyword search |
| Usage Summary | Yes | Aggregate statistics: total users, projects, experiments, storage, etc. |

This differs from the regular `project list` / `workspace info` commands, which are scoped to a single workspace. Instance-level queries have a global view.

### Quick Disambiguation (Self-Hosted)

| User says... | They probably mean... | CLI command |
|---|---|---|
| "all projects on the server" | Instance-wide project listing | `swanlab api selfhosted list-projects` |
| "all workspaces" | Instance-wide workspace listing | `swanlab api selfhosted list-workspaces` |
| "server usage / disk usage" | Usage summary | `swanlab api selfhosted summary` |
| "all users on the server" | Instance-wide user listing | `swanlab api selfhosted list-users` |

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
| "all projects on server" | Instance-wide project listing (self-hosted) | `swanlab api selfhosted list-projects` |
| "all workspaces on server" | Instance-wide workspace listing (self-hosted) | `swanlab api selfhosted list-workspaces` |
| "server usage summary" | System usage stats (self-hosted) | `swanlab api selfhosted summary` |
| "filter experiments" | Query by conditions | `swanlab api run filter -p user/project -f QUERY` |
| "best loss" | Minimum scalar value | `swanlab api run metrics PATH --keys loss` (check `min` field) |

---

## Filter Query (Experiment Filtering)

Experiments can be filtered by submitting a JSON array of filter objects to the runs/shows API. Each filter describes a condition on an experiment field.

### Filter Object Structure

Each filter object must contain four fields:

```json
{
  "key": "state",
  "type": "STABLE",
  "op": "EQ",
  "value": ["FINISHED"]
}
```

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `key` | string | yes | The field name to filter on |
| `type` | string | yes | Sidebar type: `STABLE`, `CONFIG`, or `SCALAR` |
| `op` | string | yes | Filter operator |
| `value` | array | yes | Values to compare against (always a list) |

### Filter Types (`type`)

| Type | Description | Valid Keys |
|------|-------------|------------|
| `STABLE` | Built-in experiment fields (see table below) | `state`, `name`, `description`, `show`, `pin`, `baseline`, `colors`, `cluster`, `job`, `createdAt`, `updatedAt`, `finishedAt`, `pinnedAt`, `labels` |
| `CONFIG` | Hyperparameters from `swanlab.init(config={...})` | Any key from the config dict, e.g. `learning_rate`, `batch_size` |
| `SCALAR` | Scalar metrics from `swanlab.log(...)` | Any metric key, e.g. `train/loss`, `val/acc` |

### STABLE Keys

| Key | Type | Description |
|-----|------|-------------|
| `state` | enum | Experiment state: `RUNNING`, `FINISHED`, `CRASHED`, `ABORTED`, `OFFLINE` |
| `name` | string | Experiment display name |
| `description` | string | Experiment description text |
| `show` | boolean | Whether the experiment is visible in the UI |
| `pin` | boolean | Whether the experiment is pinned to the top |
| `baseline` | boolean | Whether this experiment is marked as a baseline for comparison |
| `colors` | string | Assigned color identifier |
| `cluster` | string | Experiment group/cluster name |
| `job` | string | Distributed job type (e.g. `worker`, `master`) |
| `createdAt` | datetime | Time the experiment was created (UTC+8 on filter, UTC in response) |
| `updatedAt` | datetime | Time the experiment was last updated (UTC+8 on filter, UTC in response) |
| `finishedAt` | datetime | Time the experiment finished (UTC+8 on filter, UTC in response) |
| `pinnedAt` | datetime | Time the experiment was pinned (UTC+8 on filter, UTC in response) |
| `labels` | array | Tags/labels attached to the experiment |

### Filter Operators (`op`)

| Operator | Meaning | Notes |
|----------|---------|-------|
| `EQ` | Equals | Single value |
| `NEQ` | Not equals | Single value |
| `GTE` | Greater than or equal | Numeric / date / string comparison |
| `LTE` | Less than or equal | Numeric / date / string comparison |
| `IN` | In list | Value is an array of options |
| `NOT IN` | Not in list | Value is an array of options |
| `CONTAIN` | Fuzzy contains | Substring / partial match |

### Time Fields & Range Queries

Time-type fields (`createdAt`, `updatedAt`, `finishedAt`, `pinnedAt`) are **range-queried** using `GTE` and `LTE` together to define a time interval.

**Important timezone note:**
- **Filter values** use **UTC+8** (Asia/Shanghai). Input times are interpreted as Beijing time.
- **Response timestamps** are in **UTC**. The returned experiment fields (`createdAt`, `finishedAt`, etc.) are UTC ISO strings.

To filter experiments within a time range, combine `GTE` and `LTE` filters on the same key:

```json
[
  {"key": "createdAt", "type": "STABLE", "op": "GTE", "value": ["2026-01-01T00:00:00"]},
  {"key": "createdAt", "type": "STABLE", "op": "LTE", "value": ["2026-04-29T23:59:59"]}
]
```

- **Value format**: ISO 8601 string (e.g. `"2026-01-01T00:00:00"`). The server interprets filter values as UTC+8 and compares as date objects.
- When reading results, remember the returned timestamps are UTC — they will be 8 hours behind the filter input times.

### Operator Constraints

- Array-type fields (e.g. `labels`) only support `EQ`, `NEQ`, `IN`, `NOT IN`, `CONTAIN`.
- Date-type fields use range queries with `GTE`/`LTE` (see above).
- Numeric fields prefer numeric comparison, falling back to string comparison.

### Examples

Filter by state:
```json
[{"key": "state", "type": "STABLE", "op": "IN", "value": ["FINISHED", "RUNNING"]}]
```

Filter by config value:
```json
[{"key": "learning_rate", "type": "CONFIG", "op": "EQ", "value": ["0.01"]}]
```

Filter by scalar metric:
```json
[{"key": "train/loss", "type": "SCALAR", "op": "LTE", "value": ["0.5"]}]
```

Filter by time range (filter values in UTC+8, response in UTC):
```json
[{"key": "finishedAt", "type": "STABLE", "op": "GTE", "value": ["2026-04-01T00:00:00"]}]
```

Multiple filters combined:
```json
[
  {"key": "state", "type": "STABLE", "op": "EQ", "value": ["FINISHED"]},
  {"key": "learning_rate", "type": "CONFIG", "op": "GTE", "value": ["0.001"]},
  {"key": "train/loss", "type": "SCALAR", "op": "LTE", "value": ["0.5"]}
]
```

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
