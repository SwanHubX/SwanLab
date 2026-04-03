import os
import pathlib

from swanlab.core_python.save import collect_save_files


def test_collect_save_files_expands_directories_and_preserves_base_path(tmp_path):
    base_path = tmp_path / "playground"
    nested_dir = base_path / "hahah" / "nested"
    nested_dir.mkdir(parents=True)

    root_file = base_path / "hahah" / "root.txt"
    nested_file = nested_dir / "child.txt"
    root_file.write_text("root")
    nested_file.write_text("child")

    files = collect_save_files(
        pathlib.PurePath(os.path.abspath(base_path / "hahah" / "*")),
        pathlib.PurePath(os.path.abspath(base_path)),
        tmp_path / "run_files",
    )

    assert [file.name for file in files] == [
        "hahah/nested/child.txt",
        "hahah/root.txt",
    ]
    assert [file.source_path for file in files] == [
        str(nested_file.absolute()),
        str(root_file.absolute()),
    ]
    assert [file.target_path for file in files] == [
        str((tmp_path / "run_files" / "hahah" / "nested" / "child.txt").absolute()),
        str((tmp_path / "run_files" / "hahah" / "root.txt").absolute()),
    ]
