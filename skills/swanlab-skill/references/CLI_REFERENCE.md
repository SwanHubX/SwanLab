# SwanLab CLI Reference (`swanlab api`)

Complete reference for all `swanlab api` CLI subcommands. All commands output JSON to stdout. Use `--save` to persist results to a file.

**Prerequisite**: The user must be logged in (`swanlab login`), or supply `--host` / `--api-key` explicitly.

Read `references/SWANLAB_CONCEPTS.md` when you need to understand the data model, entity hierarchy, or disambiguate user requests.

---

## Common Options

Every command accepts these global authentication/override options:

| Option | Short | Type | Description |
|--------|-------|------|-------------|
| `--host` | `-h` | string | SwanLab server host URL. Defaults to the logged-in host. |
| `--api-key` | `-k` | string | API key for authentication. Defaults to the logged-in key. |
| `--save` | | flag | Save JSON output to a file in the current directory. Use `--save <filename>` for a custom name. When used bare, auto-generates `swanlab-YYYYMMDD_HHMMSS-xxxx.json`. |

## Path Convention

Several commands take a `PATH` argument. See `references/SWANLAB_CONCEPTS.md > Path Convention` for the full format and rules. In short:

- **Project path**: `username/project_name`
- **Experiment path**: `username/project_name/run_id`

---

## Commands Reference

### 1. User

#### `swanlab api user info`

Get the currently authenticated user's profile.

```bash
swanlab api user info [--save [FILENAME]] [--host HOST] [--api-key KEY]
```

---

### 2. Workspace

#### `swanlab api workspace info USERNAME`

Get workspace details for a given username.

```bash
swanlab api workspace info USERNAME [--save [FILENAME]] [--host HOST] [--api-key KEY]
```

| Argument | Required | Description |
|----------|----------|-------------|
| `USERNAME` | yes | The workspace username to look up |

---

### 3. Project

#### `swanlab api project info PATH`

Get project details by path.

```bash
swanlab api project info PATH [--save [FILENAME]] [--host HOST] [--api-key KEY]
```

| Argument | Required | Description |
|----------|----------|-------------|
| `PATH` | yes | Project path in `username/project_name` format |

#### `swanlab api project list`

List projects under a workspace. Paginated by default.

```bash
swanlab api project list [OPTIONS]
```

| Option | Short | Default | Description |
|--------|-------|---------|-------------|
| `--page_num` | `-n` | 1 | Page number (>= 1) |
| `--page_size` | `-s` | 20 | Page size. One of: 10, 12, 15, 20, 24, 27, 50, 100 |
| `--workspace` | | logged-in user | Workspace username to list projects from |
| `--all` | | false | Fetch all pages at once |
| `--save` | | off | Save output to file |

#### `swanlab api project create`

Create a new project.

```bash
swanlab api project create -n NAME [OPTIONS]
```

| Option | Short | Default | Description |
|--------|-------|---------|-------------|
| `--name` | `-n` | (required) | Project name. 1-100 chars, only `0-9a-zA-Z-_.+` |
| `--visibility` | `-v` | PRIVATE | `PUBLIC` or `PRIVATE` |
| `--description` | `-d` | none | Project description |
| `--workspace` | `-w` | logged-in user | Target workspace username |
| `--save` | | off | Save output to file |

---

### 4. Experiment (Run)

#### `swanlab api run info PATH`

Get experiment details.

```bash
swanlab api run info PATH [--save [FILENAME]] [--host HOST] [--api-key KEY]
```

| Argument | Required | Description |
|----------|----------|-------------|
| `PATH` | yes | Experiment path in `username/project_name/run_id` format |

#### `swanlab api run list`

List experiments under a project. Paginated by default.

```bash
swanlab api run list -p PROJECT_PATH [OPTIONS]
```

| Option | Short | Default | Description |
|--------|-------|---------|-------------|
| `--project_path` | `-p` | (required) | Project path in `username/project_name` format |
| `--page_num` | `-n` | 1 | Page number (>= 1) |
| `--page_size` | `-s` | 20 | Page size. One of: 10, 12, 15, 20, 24, 27, 50, 100 |
| `--all` | | false | Fetch all pages |
| `--save` | | off | Save output to file |

#### `swanlab api run filter -p PROJECT_PATH -f FILTER_QUERY`

Filter experiments under a project by a structured query. Returns matching experiments without pagination.

See `references/SWANLAB_CONCEPTS.md > Filter Query` for the full filter object structure, supported types, operators, and constraints.

```bash
swanlab api run filter -p PROJECT_PATH -f FILTER_QUERY [--save [FILENAME]] [--host HOST] [--api-key KEY]
```

| Option | Short | Default | Description |
|--------|-------|---------|-------------|
| `--project_path` | `-p` | (required) | Project path in `username/project_name` format |
| `--filter_query` | `-f` | (required) | Filter query as an inline JSON string or path to a `.json` file |
| `--save` | | off | Save output to file |

**Quick examples:**

```bash
# Filter by state
swanlab api run filter -p user/project -f '[{"key":"state","type":"STABLE","op":"EQ","value":["FINISHED"]}]'

# From a JSON file
swanlab api run filter -p user/project -f ./filters.json
```

For full filter syntax, operators, and advanced examples (config, scalar, time range, multi-condition), see `SWANLAB_CONCEPTS.md > Filter Query`.

#### `swanlab api run columns PATH`

List columns (metric definitions) under an experiment.

```bash
swanlab api run columns PATH [OPTIONS]
```

| Argument | Required | Description |
|----------|----------|-------------|
| `PATH` | yes | Experiment path in `username/project_name/run_id` format |

| Option | Short | Default | Description |
|--------|-------|---------|-------------|
| `--page_num` | `-n` | 1 | Page number (>= 1) |
| `--page_size` | `-s` | 20 | Page size. One of: 10, 12, 15, 20, 24, 27, 50, 100 |
| `--search` | | none | Fuzzy search keyword (matches column **name**, case-insensitive) |
| `--class` | | CUSTOM | Column class: `CUSTOM` or `SYSTEM` (see `SWANLAB_CONCEPTS.md > Column Classes`) |
| `--type` | | all | Data type filter (see `SWANLAB_CONCEPTS.md > Column Data Types` for all valid values) |
| `--all` | | false | Fetch all pages |
| `--save` | | off | Save output to file |

> **Notes on columns**:
> - Columns are **paginated** — use `--page_num` / `--page_size` to iterate in order to fetch all pages.
> - To query system metrics (CPU, GPU, memory, etc.), you must explicitly pass `--class SYSTEM`. The default is `CUSTOM` (user-defined metrics only).
> - **Resumed experiments have no system columns.** When an experiment is resumed via `swanlab.init(resume=...)`, system metrics (hardware monitoring) are not re-collected, so `--class SYSTEM` will return an empty list.

#### `swanlab api run column PATH --key KEY`

Get a single column by key name.

```bash
swanlab api run column PATH --key KEY [OPTIONS]
```

| Argument | Required | Description |
|----------|----------|-------------|
| `PATH` | yes | Experiment path in `username/project_name/run_id` format |

| Option | Short | Default | Description |
|--------|-------|---------|-------------|
| `--key` | | (required) | Column key name |
| `--class` | | CUSTOM | `CUSTOM` or `SYSTEM` (see `SWANLAB_CONCEPTS.md > Column Classes`) |
| `--type` | | all | Data type filter (see `SWANLAB_CONCEPTS.md > Column Data Types`) |
| `--save` | | off | Save output to file |

#### `swanlab api run metrics PATH --keys KEYS`

Get scalar metric data for specified keys.
Keys is a comma-separated list of column keys, e.g. `loss,acc`.
Returns an object mapping each key to its data points (step, value). Use `--all` to fetch all data points without sampling limit. See `SWANLAB_CONCEPTS.md > Scalar Metrics` for data structure details.

```bash
swanlab api run metrics PATH --keys KEYS [OPTIONS]
```

| Argument | Required | Description |
|----------|----------|-------------|
| `PATH` | yes | Experiment path in `username/project_name/run_id` format |

| Option | Short | Default | Description |
|--------|-------|---------|-------------|
| `--keys` | | (required) | Comma-separated metric keys, e.g. `loss,acc` |
| `--sample` | `-s` | 1500 | Sample size for scalars (>= 1). Max 1500; use `--all` for full export. |
| `--ignore-timestamp` | | false | Strip timestamps from metric data |
| `--all` | | false | Fetch all data points without sampling limit |
| `--range-type` | | `step` | Range filter axis: `step` or `timestamp`. Defaults to `step` when any range option is provided. |
| `--range-start` | | none | Range start value (inclusive). For `step` type: step number (int >= 0). For `timestamp` type: Unix timestamp in **milliseconds** (int >= 0). |
| `--range-end` | | none | Range end value (inclusive). Same type as `--range-start`. |
| `--range-head` | | none | Return only the first N data points (int >= 1). Mutually exclusive with `--range-tail`. |
| `--range-tail` | | none | Return only the last N data points (int >= 1). Mutually exclusive with `--range-head`. |
| `--range-last` | | none | Last N milliseconds of data (int >= 1). Mutually exclusive with `--range-start`/`--range-end`. Can combine with `--range-head`/`--range-tail`. |
| `--save` | | off | Save output to file |

**Range query constraints:**
- `--range-head` and `--range-tail` are mutually exclusive.
- `--range-last` is mutually exclusive with `--range-start`/`--range-end`.
- `--range-start` must be ≤ `--range-end`.
- `--range-head`/`--range-tail` can be combined with `--range-last` or `--range-start`/`--range-end`.
- Range query is only supported for SCALAR metrics. It downloads CSV data and applies client-side filtering.
- When using `--range-type timestamp`, each CSV row must have a timestamp column. Rows missing timestamps are skipped with a warning.

**Quick examples:**

```bash
# Get steps 100–200
swanlab api run metrics PATH --keys loss --range-start 100 --range-end 200

# Get first 50 data points
swanlab api run metrics PATH --keys loss --range-head 50

# Get data by timestamp range (Unix milliseconds)
swanlab api run metrics PATH --keys loss --range-type timestamp --range-start 1714368000000 --range-end 1714454400000

# Get data from the last 5 minutes
swanlab api run metrics PATH --keys loss --range-last 300000

# Get last 30 data points
swanlab api run metrics PATH --keys loss --range-tail 30
```

> **Tip**: For large metric data, use `--save` to write JSON to file, then plot with `scripts/plot_metrics.py --data file.json -k loss`. For quick stats, prefer `run summary` over full metrics.

#### `swanlab api run summary PATH`

Get scalar metric summaries (statistics) for an experiment. Returns per-key aggregates: last step/value, min/max, avg, median, stdDev. See `SWANLAB_CONCEPTS.md > Scalar Summary` for response structure details.

```bash
swanlab api run summary PATH [--keys KEYS] [--save [FILENAME]] [--host HOST] [--api-key KEY]
```

| Argument | Required | Description |
|----------|----------|-------------|
| `PATH` | yes | Experiment path in `username/project_name/run_id` format |

| Option | Short | Default | Description |
|--------|-------|---------|-------------|
| `--keys` | | all | Comma-separated scalar keys, e.g. `loss,acc`. Omit to query all scalar keys. |
| `--save` | | off | Save output to file |

#### `swanlab api run medias PATH --keys KEYS`

Get media data for specified keys. See `SWANLAB_CONCEPTS.md > Media Metrics` for the response structure and how presigned URLs work.

```bash
swanlab api run medias PATH --keys KEYS [OPTIONS]
```

| Argument | Required | Description |
|----------|----------|-------------|
| `PATH` | yes | Experiment path in `username/project_name/run_id` format |

| Option | Short | Default | Description |
|--------|-------|---------|-------------|
| `--keys` | | (required) | Comma-separated media keys, e.g. `image,audio` |
| `--step` | `-s` | 0 | Step number to retrieve |
| `--all` | | false | Fetch all steps |
| `--save` | | off | Save output to file |

#### `swanlab api run logs PATH`

Get console logs captured during a run. See `SWANLAB_CONCEPTS.md > Console Logs` for shard-based pagination details.

```bash
swanlab api run logs PATH [OPTIONS]
```

| Argument | Required | Description |
|----------|----------|-------------|
| `PATH` | yes | Experiment path in `username/project_name/run_id` format |

| Option | Short | Default | Description |
|--------|-------|---------|-------------|
| `--offset` | `-o` | 0 | Log offset (shard index) |
| `--level` | `-l` | INFO | Log level: `DEBUG`, `INFO`, `WARN`, `ERROR` |
| `--ignore-timestamp` | | false | Strip timestamps from log entries |
| `--save` | | off | Save output to file |

#### `swanlab api run export-logs PATH`

Export experiment console logs as a downloadable `.log` file. Returns a presigned download URL. See `SWANLAB_CONCEPTS.md > Log Export` for conceptual details.

```bash
swanlab api run export-logs PATH [OPTIONS]
```

| Argument | Required | Description |
|----------|----------|-------------|
| `PATH` | yes | Experiment path in `username/project_name/run_id` format |

| Option | Short | Default | Description |
|--------|-------|---------|-------------|
| `--start` | | 0 | Start row index (0-based) |
| `--rows` | `-r` | 500000 | Number of rows to export (1–500000) |
| `--save` | | off | Save the export URL as JSON to file |

---

### 5. Self-Hosted

These commands are only valid for **self-hosted SwanLab instances**. All self-hosted management commands require **root (admin) privileges**. See `SWANLAB_CONCEPTS.md > Self-Hosted Instance` for the host detection rule — do not use these commands if the resolved host contains `swanlab.cn`.

#### `swanlab api self-hosted info`

Show self-hosted instance info.

```bash
swanlab api self-hosted info [--save [FILENAME]] [--host HOST] [--api-key KEY]
```

#### `swanlab api self-hosted create-user`

Create a user in the self-hosted instance (root only).

```bash
swanlab api self-hosted create-user -u USERNAME -p PASSWORD
```

| Option | Short | Description |
|--------|-------|-------------|
| `--username` | `-u` | Username to create (required) |
| `--password` | `-p` | Password for the new user (required) |

#### `swanlab api self-hosted list-users`

List users in the self-hosted instance (root only).

```bash
swanlab api self-hosted list-users [OPTIONS]
```

| Option | Short | Default | Description |
|--------|-------|---------|-------------|
| `--page_num` | `-n` | 1 | Page number |
| `--page_size` | `-s` | 20 | Page size |
| `--all` | | false | Fetch all users |

#### `swanlab api self-hosted list-projects`

List all projects across the self-hosted instance (root only). Supports filtering by creator, workspace, and keyword search.

```bash
swanlab api self-hosted list-projects [OPTIONS]
```

| Option | Short | Default | Description |
|--------|-------|---------|-------------|
| `--page_num` | `-n` | 1 | Page number |
| `--page_size` | `-s` | 20 | Page size |
| `--all` | | false | Fetch all projects |
| `--search` | | none | Search keyword |
| `--creator` | | none | Filter by creator username |
| `--workspace` | | none | Filter by workspace username |

#### `swanlab api self-hosted list-workspaces`

List all workspaces in the self-hosted instance (root only).

```bash
swanlab api self-hosted list-workspaces [OPTIONS]
```

| Option | Short | Default | Description |
|--------|-------|---------|-------------|
| `--page_num` | `-n` | 1 | Page number |
| `--page_size` | `-s` | 20 | Page size |
| `--all` | | false | Fetch all workspaces |
| `--search` | | none | Search keyword |

#### `swanlab api self-hosted summary`

Show system usage summary for the self-hosted instance (root only). Includes aggregate statistics such as total users, projects, experiments, and storage usage. See `SWANLAB_CONCEPTS.md > Self-Hosted Instance > Usage Summary` for details.

```bash
swanlab api self-hosted summary [--save [FILENAME]] [--host HOST] [--api-key KEY]
```

---

## Behavioral Constraints

- **Use `--all` only when the user explicitly asks for it.** This flag bypasses pagination and fetches the entire dataset in one go, which puts heavy load on the database. For paginated list commands, always use default pagination (`--page_num` / `--page_size`). Only add `--all` when the user says something like "fetch all", "get everything", or "I want the complete list".
- **Always ask the user for specific column keys before running `run metrics`, `run medias`, or `run column`.** These commands accept a `--keys` or `--key` parameter. Querying without a specific key forces the server to scan and return data for all columns, which is expensive. If the user doesn't know the key names, first run `run columns PATH` to list available columns, then use the returned keys for targeted queries. Always discover key names from `run columns` output before querying specific metrics.
- **Always add `--ignore-timestamp` for `run metrics` and `run logs`.** Unless the user specifically asks to keep timestamps, include this flag by default. It produces cleaner, more readable output by removing Unix timestamps from every data point.
- All output is JSON to stdout. Pipe to `jq` or similar tools for further processing.
- `--save` without a filename auto-generates `swanlab-YYYYMMDD_HHMMSS-xxxx.json` in the current directory.
- If not logged in, supply both `--host` and `--api-key` on every invocation. See `SWANLAB_CONCEPTS.md > Key Environment Variables` for the full variable list.
