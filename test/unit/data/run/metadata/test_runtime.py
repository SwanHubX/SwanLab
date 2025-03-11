from swanlab.data.run.metadata.runtime import parse_git_url


def test_parse_git_url():
    # ssh
    assert parse_git_url("git@github.com:swanhubx/swanlab.git") == "https://github.com/swanhubx/swanlab"
    assert parse_git_url("git@localhost:8000/swanhubx/swanlab.git") == "https://localhost:8000/swanhubx/swanlab"
    # https
    assert parse_git_url("https://github.com/swanhubx/swanlab.git") == "https://github.com/swanhubx/swanlab"
    assert parse_git_url("https://localhost:8000/swanhubx/swanlab.git") == "https://localhost:8000/swanhubx/swanlab"
    # no .git
    assert parse_git_url("git@github.com:swanhubx/swanlab") == "https://github.com/swanhubx/swanlab"
    assert parse_git_url("git@localhost:8000/swanhubx/swanlab") == "https://localhost:8000/swanhubx/swanlab"
    assert parse_git_url("https://github.com/swanhubx/swanlab") == "https://github.com/swanhubx/swanlab"
    assert parse_git_url("https://localhost:8000/swanhubx/swanlab") == "https://localhost:8000/swanhubx/swanlab"
    