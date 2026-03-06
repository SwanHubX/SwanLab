.PHONY: format

format:
	-uvx ruff check --select I --fix . --quiet
	-uvx ruff format . --quiet
	