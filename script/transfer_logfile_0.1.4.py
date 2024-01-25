#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-25 12:37:08
@File: test/test_wandb.py
@IDE: vscode
@Description:
    Swanlab update script for upgrading from v0.1.4 to v0.1.5, transitioning from file system to database.
    This script depends on swanlab v0.1.5, so please upgrade to swanlab v0.1.5 and backup the original log data before using this script—even though this script will not delete the original log data.

    The following is the usage method:
        python transfer_logfile_0.1.4.py -i /path/to/swanlog_old -o /path/to/swanlog_new
"""
import os
from swanlab.env import ROOT
from swanlab.db import connect, Experiment, Project, Chart, Tag, Source, Namespace, Display
import json
import yaml
import argparse


def run(old_path: str, new_dir: str = "swanlog.new"):
    if not os.path.exists(old_path):
        raise FileNotFoundError(f"{old_path} not found")
    if not os.path.isdir(old_path):
        raise FileNotFoundError(f"{old_path} is not a directory")
    # Parent directory of the old log path
    old_parent_dir = os.path.dirname(old_path)
    # New log path
    new_path = os.path.join(old_parent_dir, new_dir)
    # Create new log path
    if os.path.exists(new_path):
        raise FileExistsError(f"{new_path} already exists")
        # shutil.rmtree(new_path)
        # print(f"{new_path} already exists, remove it")
    os.mkdir(new_path)
    print("Begin to migrate swanlog...")
    # Inject environment variable
    os.environ[ROOT] = new_path
    # Next, create the database file, set up tables, and migrate data
    db = connect(autocreate=True)
    # Copy all folders from the old directory to the new directory
    for name in os.listdir(old_path):
        # Check if it's a folder
        if os.path.isdir(os.path.join(old_path, name)):
            os.system(f"cp -r {os.path.join(old_path, name)} {new_path}")
    # Project information in the old directory
    old_projects = json.loads(open(os.path.join(old_path, "project.json"), "r").read())
    name = old_projects.get("name", "")
    Project.init(name=name)
    # Iterate through all experiments in project.json
    for experiment in old_projects.get("experiments", []):
        # Create experiment
        exp: Experiment = Experiment.create(experiment["name"], experiment["name"], experiment.get("description", ""))
        print(f"migrate {experiment['name']}...")
        # Create experiment configuration in the corresponding folder files/config.yaml
        if not os.path.exists(os.path.join(new_path, experiment["name"], "files")):
            os.mkdir(os.path.join(new_path, experiment["name"], "files"))
        with open(os.path.join(new_path, experiment["name"], "files", "config.yaml"), "w") as f:
            # Iterate through the experiment's configuration fields, change the original value to the value field, and add a desc field
            config = experiment.get("config", {})
            for key in config.keys():
                config[key] = {"value": config[key], "desc": ""}
            yaml.dump(config, f)
            print(f"create config.yaml in {experiment['name']}")
        # Iterate through all charts of the experiment, located in the experiment folder's chart.json
        if not os.path.exists(os.path.join(new_path, experiment["name"], "chart.json")):
            print(f"chart.json not found in {experiment['name']}")
            continue
        # Change experiment status
        exp.status = experiment.get("status", -1)
        # Save
        exp.save()
        print(f"update status of {experiment['name']} to {experiment.get('status', -1)}")
        # Read chart.json
        charts = json.loads(open(os.path.join(new_path, experiment["name"], "chart.json"), "r").read())
        # Generate chart, tag, and source
        chart_dict = {}
        for chart in charts.get("charts", []):
            # Create chart
            c = Chart.create(
                chart["source"][0],
                experiment_id=exp.id,
                type=chart.get("type", "default"),
                reference=chart.get("reference", "default"),
                config=chart.get("config", None),
            )
            chart_dict[chart["chart_id"]] = c
            # Get error
            error = chart.get("error", None)
            # Create tag
            tag = Tag.create(exp.id, chart["source"][0], type="number")
            Source.create(tag.id, c.id, error)
        # Generate namespace
        for namespace in charts.get("namespaces", []):
            n = Namespace.create(namespace["namespace"], exp.id)
            # Iterate through charts
            for chart_id in namespace.get("charts", []):
                # Create display
                Display.create(chart_dict[chart_id].id, n.id)


if __name__ == "__main__":
    parse = argparse.ArgumentParser()
    # Old log path, to avoid ambiguity, please provide an absolute path, e.g., /home/username/swanlog
    parse.add_argument("-i", "--old_path", type=str, default="path/to/swanlog_old")
    # This script will create a new log folder named swanlog.new in the upper directory of the old log path
    # Please ensure that the upper directory of the old log path has write permission, and there is no swanlog.new folder under the old log path
    parse.add_argument("-o", "--new_dir", type=str, default="path/to/swanlog_new")
    args = parse.parse_args()

    _old_path = args.old_path
    _new_dir = args.new_dir

    if not os.path.isabs(_old_path):
        _old_path = os.path.abspath(_old_path)

    if not os.path.isabs(_new_dir):
        _new_path = os.path.abspath(_new_dir)

    run(_old_path, _new_path)
