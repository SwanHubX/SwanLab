# SwanLab SDK Quick Start

Reference for writing SwanLab tracking code. Based on the current SDK API surface.
All data **reading** operations should use the `swanlab api` CLI — see `CLI_REFERENCE.md`.

---

## 1. Install & Login

```bash
pip install swanlab
```

Login (one-time):

```bash
swanlab login
```

Or login programmatically:

```python
import swanlab
swanlab.login(api_key="your-api-key", save=True)
```

---

## 2. Create an Experiment — `swanlab.init()`

```python
run = swanlab.init(
    project="my-project",        # project name (defaults to cwd name)
    name="resnet-exp-1",         # experiment name (auto-generated if omitted)
    description="baseline run",  # experiment description
    config={                     # hyperparameters / metadata (dict, JSON path, or YAML path)
        "learning_rate": 0.01,
        "epochs": 10,
        "batch_size": 32,
    },
)
```

### Key Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `project` | `str` | cwd name | Project name |
| `name` | `str` | auto-generated | Experiment display name |
| `description` | `str` | none | Experiment description |
| `config` | `dict` / `str` / `Path` | none | Hyperparameters dict or path to JSON/YAML file |
| `mode` | `str` | `"online"` | Run mode — see SKILL.md for details |
| `workspace` | `str` | current user | Target workspace username |
| `tags` | `list[str]` | none | Tags for categorization |
| `group` | `str` | none | Group name for distributed/related experiments — see §11 |
| `job_type` | `str` | none | Role label, e.g. `"train"`, `"worker"` — see §11 |
| `resume` | `str` / `bool` | `"never"` | Resume strategy — see §9 |
| `id` | `str` | auto | Custom run ID (1–512 chars, no `<>:"/\|?*#%` or control chars). Required when `resume="must"`. |
| `log_dir` | `str` | `"./swanlog"` | Local log directory |
| `callbacks` | `list` | none | Callback objects for lifecycle hooks |
| `reinit` | `bool` | `False` | If True, finish current run before starting new one |

### Access Config at Runtime

`run.config` is a mutable mapping with attribute access. Read hyperparameters after init:

```python
lr = run.config.learning_rate
epochs = run.config["epochs"]
```

---

## 3. Log Metrics — `swanlab.log()`

Log a dict of scalar values. `step` auto-increments if omitted.

```python
swanlab.log({"loss": 0.5, "accuracy": 0.92})
swanlab.log({"loss": 0.3, "accuracy": 0.95}, step=10)  # explicit step
```

- Keys support `/` for nested grouping: `swanlab.log({"train/loss": 0.5, "val/loss": 0.6})`
- Nested dicts are auto-flattened: `{"train": {"loss": 0.5}}` → `"train.loss"`
- Values can be `int`, `float`, or `str`
- `step` must be a non-negative integer

---

## 4. Log Media — Two Styles

### Style A: Pass media objects to `swanlab.log()`

```python
swanlab.log({
    "sample_image": swanlab.Image(tensor_or_path),
    "sample_audio": swanlab.Audio(numpy_array, sample_rate=16000),
    "sample_text":  swanlab.Text("generated output"),
    "sample_video": swanlab.Video("video.mp4"),
    "chart":        swanlab.ECharts(pyecharts_chart),
    "point_cloud":  swanlab.Object3D(numpy_array),
    "mol":          swanlab.Molecule("CCO"),  # SMILES string
})
```

### Style B: Use convenience functions

```python
swanlab.log_image(key="image", data=pil_image, caption="epoch 5")
swanlab.log_audio(key="audio", data="audio.wav", sample_rate=44100)
swanlab.log_text(key="output", data="hello world")
swanlab.log_video(key="demo", data="demo.mp4")
swanlab.log_echarts(key="chart", data=pyecharts_obj)
swanlab.log_object3d(key="pc", data=numpy_array)
swanlab.log_molecule(key="mol", data="CCO")
```

### Media Input Types

| Media Class | Accepted Input Types |
|-------------|---------------------|
| `Image` | PIL Image, numpy array, torch Tensor, matplotlib Figure, file path (`str`) |
| `Audio` | numpy array, file path (`str`) |
| `Video` | file path (`str`) |
| `Text` | string (`str`) |
| `ECharts` | pyecharts chart object |
| `Object3D` | numpy array, dict, file path (`str`) |
| `Molecule` | SMILES string, file path (`str`), RDKit Mol object |

All media `log_*` functions accept optional `caption` (string) and `step` (int) parameters.

---

## 5. Finish an Experiment — `swanlab.finish()`

```python
swanlab.finish()                          # state="success" (default)
swanlab.finish(state="crashed", error="OOM")  # report a crash
```

- `state`: `"success"` / `"crashed"` / `"aborted"`
- `error`: optional error message (required when `state="crashed"`)
- `async_log_timeout`: optional timeout (seconds) for pending async_log tasks
- SwanLab auto-calls `finish()` at program exit if not called explicitly

---

## 6. Other SDK Utilities

### Async Logging — `swanlab.async_log()`

Run a function in the background, auto-log its return dict when done:

```python
future = swanlab.async_log(compute_expensive_metric, x=data, mode="threading")
# finish() waits for all pending async_log tasks
```

| `mode` | Execution |
|--------|-----------|
| `"threading"` | Background thread (default) |
| `"spawn"` | New process |
| `"fork"` | Forked process |
| `"asyncio"` | Async event loop |

---

## 7. Complete Example

```python
import swanlab
import random

# Initialize
run = swanlab.init(
    project="my-project",
    config={"learning_rate": 0.01, "epochs": 10},
)

# Training loop
for epoch in range(2, run.config.epochs):
    acc = 1 - 2**-epoch - random.random() / epoch
    loss = 2**-epoch + random.random() / epoch
    swanlab.log({"accuracy": acc, "loss": loss})

# Finish
swanlab.finish()
```

---

## 8. External Tool Conversion

SwanLab can import data from other tracking tools:

```python
from swanlab import sync_wandb, sync_tensorboardX, sync_tensorboard_torch, sync_mlflow
```

See the `converter/` module for detailed usage.

---

## 9. Resume / Continue Training — `resume` + `id`

Resume an existing experiment to append new data (instead of creating a new one).

```python
# First run — remember the ID
run = swanlab.init(project="my-project")
run_id = run.id       # e.g. "abc123"
swanlab.log({"loss": 2.0})
run.finish()

# Later — resume that same experiment
run = swanlab.init(project="my-project", resume="must", id=run_id)
swanlab.log({"loss": 0.5})   # appended after the existing data
run.finish()
```

### Resume Strategies

| `resume` | Existing experiment with same `id` | No matching experiment |
|-----------|-----------------------------------|------------------------|
| `"must"` | Resume it | **Error** |
| `"allow"` / `True` | Resume it | Create a new one |
| `"never"` / `False` | **Error** (if `id` is set) | Create a new one |

### Finding the Experiment ID

- In the UI: experiment → **Environment** tab → **Log Directory** field shows the run ID
- In the URL: `https://swanlab.cn/@<user>/<project>/runs/<exp_id>/...`
- Via CLI: `swanlab api run info user/project/run_id`

### Environment Variable Alternative

When modifying `swanlab.init()` is inconvenient (e.g. framework integration):

```bash
export SWANLAB_RESUME=must
export SWANLAB_RUN_ID=<exp_id>
```

---

## 10. Sync Offline Data — `swanlab sync`

Upload locally stored experiment data to the cloud. Works with `mode="offline"` and `mode="local"` experiments.

```bash
# Sync a single experiment
swanlab sync ./swanlog/run-abc123

# Sync to a different project / workspace
swanlab sync ./swanlog/run-abc123 -p other-project -w team-workspace

# Sync onto an existing experiment (incremental — only new data is uploaded)
swanlab sync ./swanlog/run-abc123 --id <exp_id>

# Batch sync all runs
swanlab sync ./swanlog/run-*

# Self-hosted target
swanlab sync ./swanlog/run-abc123 --host https://your-swanserver
```

Key options: `-p` (project), `-w` (workspace), `--id` (existing experiment ID), `-k` (API key), `--host`.

---

## 11. Distributed Training — `group` + `job_type`

SwanLab supports two patterns for distributed training:

### Pattern A: Log only from rank 0

Simplest approach — only the main process initializes and logs:

```python
import os, swanlab

local_rank = int(os.environ.get("LOCAL_RANK", 0))

if local_rank == 0:
    swanlab.init(project="ddp-training", name="resnet-ddp")

for epoch in range(10):
    # ... training ...
    if local_rank == 0:
        swanlab.log({"loss": loss.item()})

if local_rank == 0:
    swanlab.finish()
```

### Pattern B: Each process logs independently

Each process creates its own experiment, grouped with `group` and distinguished by `job_type`:

```python
import os, swanlab

rank = int(os.environ.get("RANK", 0))

swanlab.init(
    project="ddp-training",
    group="exp-20260604",                    # same group for all ranks
    job_type="main" if rank == 0 else f"worker-{rank}",
)

for epoch in range(10):
    # ... training ...
    swanlab.log({"loss": loss.item()})

swanlab.finish()
```

### Environment Variables for Distributed

| Variable | Purpose |
|----------|---------|
| `SWANLAB_GROUP` | Group name for associating experiments |
| `SWANLAB_JOB_TYPE` | Role label, e.g. `"train"`, `"worker"` |
| `SWANLAB_NAME` | Experiment name |
| `SWANLAB_DESCRIPTION` | Experiment description |
