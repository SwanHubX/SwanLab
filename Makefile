.PHONY: format test

sync:
	-uv sync --all-extras

format:
	-uvx ruff check --select I --fix . --quiet
	-uvx ruff format . --quiet

test:
	-uv run pytest