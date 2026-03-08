.PHONY:  init sync format test

init:
	uv sync --all-extras
	pre-commit install

sync:
	uv sync --all-extras

format:
	-uvx ruff check --select I --fix . --quiet
	-uvx ruff format . --quiet

unit:
	uv run pytest tests/unit