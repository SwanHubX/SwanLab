"""
@author: cunyue
@file: test_set_run_store.py
@time: 2025/7/21 17:19
@description: 测试 set_run_store 功能
"""

import nanoid
import pytest

from swanlab.proto.v0 import Project, Experiment
from swanlab.sync.sync_utils import set_run_store
from tutils.setup import UseMockRunState


@pytest.fixture
def temp_project() -> Project:
    """
    mock Project object for testing.
    """
    return Project.model_validate(
        {
            "name": "temp-project",
            "workspace": "temp-workspace",
            "public": False,
        }
    )


@pytest.fixture
def temp_experiment() -> Experiment:
    """
    mock Experiment object for testing.
    """
    return Experiment.model_validate(
        {
            "id": nanoid.generate(alphabet="0123456789abcdefghijklmnopqrstuvwxyz", size=21),
            "name": "temp-experiment",
            "colors": ("red", "blue"),
            "description": "Temp experiment for testing.",
            "tags": ["test", "sync"],
        }
    )


class TestSetRunStore:
    """
    测试设置运行存储的相关信息。
    """

    def test_default(self, temp_project, temp_experiment):
        with UseMockRunState() as mock:
            run_store = mock.store
            set_run_store(run_store, temp_project, temp_experiment)
            assert run_store.run_id is None, "Run ID should be None by default."
            assert run_store.resume == "never", "Resume should be set to 'never' by default."
            assert run_store.run_name == "temp-experiment", "Run name should match the experiment name."
            assert run_store.visibility is False, "Visibility should be False by default."
            assert run_store.workspace == "temp-workspace", "Workspace should match the project workspace."
            assert run_store.project == "temp-project", "Project should match the project project."
            assert run_store.run_colors == ("red", "blue"), "Run colors should match the experiment colors."
            assert run_store.tags == ["test", "sync"], "Run tags should match the experiment tags."
            assert (
                run_store.description == "Temp experiment for testing."
            ), "Run description should match the experiment description."

    def test_set_project(self, temp_project, temp_experiment):
        with UseMockRunState() as mock:
            run_store = mock.store
            set_run_store(run_store, temp_project, temp_experiment, project="custom-project")
            assert run_store.project == "custom-project", "Project should be set to 'custom-project'."

    def test_set_workspace(self, temp_project, temp_experiment):
        with UseMockRunState() as mock:
            run_store = mock.store
            set_run_store(run_store, temp_project, temp_experiment, workspace="custom-workspace")
            assert run_store.workspace == "custom-workspace", "Workspace should be set to 'custom-workspace'."

    def test_set_id_new(self, temp_project, temp_experiment):
        with UseMockRunState() as mock:
            run_store = mock.store
            set_run_store(run_store, temp_project, temp_experiment, id="new")
            assert run_store.run_id is None, "Run ID should be equal to the experiment ID."
            assert run_store.resume == "never", "Resume should be 'never' for 'new'."

    def test_set_id_auto(self, temp_project, temp_experiment):
        with UseMockRunState() as mock:
            run_store = mock.store
            set_run_store(run_store, temp_project, temp_experiment, id="auto")
            assert run_store.run_id == temp_experiment.id, "Run ID should match the experiment ID."
            assert run_store.resume == "allow", "Resume should be 'allow' for 'auto'."

    def test_set_id_custom(self, temp_project, temp_experiment):
        with UseMockRunState() as mock:
            run_store = mock.store
            custom_id = nanoid.generate(alphabet="0123456789abcdefghijklmnopqrstuvwxyz", size=21)
            set_run_store(run_store, temp_project, temp_experiment, id=custom_id)
            assert run_store.run_id == custom_id, "Run ID should be set to the custom ID."
            assert run_store.resume == "must", "Resume should be 'must' for a specific ID."
