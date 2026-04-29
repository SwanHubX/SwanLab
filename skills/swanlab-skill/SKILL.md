---
name: swanlab-skill
description: >
  Interact with SwanLab via the `swanlab api` CLI subcommands to query experiments, projects,
  workspaces, users, columns, metrics, media, logs, and self-hosted instance management.
  Use this skill whenever the user wants to retrieve SwanLab data from the command line,
  inspect experiment metrics or logs, list projects/runs/columns, manage self-hosted users,
  or automate any SwanLab API query via CLI. Also trigger when the user mentions "swanlab api",
  "swanlab cli", or asks to fetch/inspect/query SwanLab experiment data, training metrics,
  or project information from the terminal — even if they don't explicitly say "CLI" or "api".
---

# SwanLab Skill — CLI-Based Interaction

This skill enables querying the SwanLab platform via the `swanlab api` CLI group.
All commands output JSON to stdout. Use `--save` to persist results to a file.

**Prerequisite**: The user must be logged in (`swanlab login`), or supply `--host` / `--api-key` explicitly.

Read `references/SWANLAB_CONCEPTS.md` when you need to understand the data model, entity hierarchy, or disambiguate user requests.

## Common Options

Every command accepts these global authentication/override options (applied via decorator, not shown in subcommand help):

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

# Filter by config value
swanlab api run filter -p user/project -f '[{"key":"learning_rate","type":"CONFIG","op":"GTE","value":["0.001"]}]'

# Filter by time range — filter values are UTC+8, response timestamps are UTC
swanlab api run filter -p user/project -f '[
  {"key":"createdAt","type":"STABLE","op":"GTE","value":["2026-01-01T00:00:00"]},
  {"key":"createdAt","type":"STABLE","op":"LTE","value":["2026-04-29T23:59:59"]}
]'

# From a JSON file
swanlab api run filter -p user/project -f ./filters.json

# Multiple filters combined
swanlab api run filter -p user/project -f '[
  {"key":"state","type":"STABLE","op":"IN","value":["FINISHED","RUNNING"]},
  {"key":"learning_rate","type":"CONFIG","op":"EQ","value":["0.01"]}
]'
```

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
| `--class` | | CUSTOM | Column class: `CUSTOM` or `SYSTEM` (see `references/SWANLAB_CONCEPTS.md > Column Classes`) |
| `--type` | | all | Data type filter (see `references/SWANLAB_CONCEPTS.md > Column Data Types` for all valid values) |
| `--all` | | false | Fetch all pages |
| `--save` | | off | Save output to file |

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
| `--class` | | CUSTOM | `CUSTOM` or `SYSTEM` (see `references/SWANLAB_CONCEPTS.md > Column Classes`) |
| `--type` | | all | Data type filter (see `references/SWANLAB_CONCEPTS.md > Column Data Types`) |
| `--save` | | off | Save output to file |

#### `swanlab api run metrics PATH --keys KEYS`

Get scalar metric data for specified keys.
Keys is a comma-separated list of column keys, e.g. `loss,acc`.
Returns an object mapping each key to its data points (step, value). For large datasets, use `--all` to get a CSV download URL instead of inlined data. See `references/SWANLAB_CONCEPTS.md > Scalar Metrics` for data structure details.


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
| `--all` | | false | Fetch all data points (returns CSV download URL instead of inlined data) |
| `--save` | | off | Save output to file |

#### `swanlab api run medias PATH --keys KEYS`

Get media data for specified keys. See `references/SWANLAB_CONCEPTS.md > Media Metrics` for the response structure and how presigned URLs work.

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

Get console logs captured during a run. See `references/SWANLAB_CONCEPTS.md > Console Logs` for shard-based pagination details.

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

---

### 5. Self-Hosted

These commands are only valid for **self-hosted SwanLab instances**. All self-hosted management commands (create-user, list-users, list-projects, list-workspaces, summary) require **root (admin) privileges**. See `references/SWANLAB_CONCEPTS.md > Self-Hosted Instance` for the host detection rule — do not use these commands if the resolved host contains `swanlab.cn`.

#### `swanlab api selfhosted info`

Show self-hosted instance info.

```bash
swanlab api selfhosted info [--save [FILENAME]] [--host HOST] [--api-key KEY]
```

#### `swanlab api selfhosted create-user`

Create a user in the self-hosted instance (root only).

```bash
swanlab api selfhosted create-user -u USERNAME -p PASSWORD
```

| Option | Short | Description |
|--------|-------|-------------|
| `--username` | `-u` | Username to create (required) |
| `--password` | `-p` | Password for the new user (required) |

#### `swanlab api selfhosted list-users`

List users in the self-hosted instance (root only).

```bash
swanlab api selfhosted list-users [OPTIONS]
```

| Option | Short | Default | Description |
|--------|-------|---------|-------------|
| `--page_num` | `-n` | 1 | Page number |
| `--page_size` | `-s` | 20 | Page size |
| `--all` | | false | Fetch all users |

#### `swanlab api selfhosted list-projects`

List all projects across the self-hosted instance (root only). Supports filtering by creator, workspace, and keyword search.

```bash
swanlab api selfhosted list-projects [OPTIONS]
```

| Option | Short | Default | Description |
|--------|-------|---------|-------------|
| `--page_num` | `-n` | 1 | Page number |
| `--page_size` | `-s` | 20 | Page size |
| `--all` | | false | Fetch all projects |
| `--search` | | none | Search keyword |
| `--creator` | | none | Filter by creator username |
| `--workspace` | | none | Filter by workspace username |

#### `swanlab api selfhosted list-workspaces`

List all workspaces in the self-hosted instance (root only).

```bash
swanlab api selfhosted list-workspaces [OPTIONS]
```

| Option | Short | Default | Description |
|--------|-------|---------|-------------|
| `--page_num` | `-n` | 1 | Page number |
| `--page_size` | `-s` | 20 | Page size |
| `--all` | | false | Fetch all workspaces |
| `--search` | | none | Search keyword |

#### `swanlab api selfhosted summary`

Show system usage summary for the self-hosted instance (root only). Includes aggregate statistics such as total users, projects, experiments, and storage usage. See `references/SWANLAB_CONCEPTS.md > Self-Hosted Instance > Usage Summary` for details.

```bash
swanlab api selfhosted summary [--save [FILENAME]] [--host HOST] [--api-key KEY]
```

---

## Behavioral Constraints

- **Never use `--all` unless the user explicitly asks for it.** This flag bypasses pagination and fetches the entire dataset in one go, which puts heavy load on the database. For paginated list commands, always use default pagination (`--page_num` / `--page_size`). Only add `--all` if the user says something like "fetch all", "get everything", or "I want the complete list".
- **Always ask the user for specific column keys before running `run metrics`, `run medias`, or `run column`.** These commands accept a `--keys` or `--key` parameter. Querying without a specific key forces the server to scan and return data for all columns, which is expensive. If the user doesn't know the key names, first run `run columns PATH` to list available columns, then use the returned keys to make targeted queries. Never guess or fabricate key names — always discover them first.
- **Always add `--ignore-timestamp` for `run metrics` and `run logs`.** Unless the user specifically asks to keep timestamps, include this flag by default. It produces cleaner, more readable output by removing Unix timestamps from every data point.
- All output is JSON to stdout. Pipe to `jq` or similar tools for further processing.
- `--save` without a filename auto-generates `swanlab-YYYYMMDD_HHMMSS-xxxx.json` in the current directory.
- If not logged in, supply both `--host` and `--api-key` on every invocation. See `references/SWANLAB_CONCEPTS.md > Key Environment Variables` for the full variable list.
