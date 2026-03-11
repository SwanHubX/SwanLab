.PHONY:  init sync format test proto

init:
	uv sync --all-extras
	pre-commit install

sync:
	uv sync --all-extras

proto:
	uv run scripts/generate_protos.py

format:
	-uvx ruff check --select I --fix . --quiet
	-uvx ruff format . --quiet

unit:
	uv run pytest tests/unit