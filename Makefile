.PHONY:  init sync format proto unit bench clean build

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

bench:
	uv run pytest tests/benchmark

build:
	@if [ -n "$(VERSION)" ]; then \
		python -c "import json; data=json.load(open('swanlab/package.json')); data['version']='$(VERSION)'; json.dump(data,open('swanlab/package.json','w'),indent=2)" && \
		echo "Updated swanlab/package.json version to $(VERSION)"; \
	else \
		echo "VERSION not set, using default version in swanlab/package.json"; \
	fi
	uv build

publish:
	@if [ -z "$$PYPI_TOKEN" ]; then \
		echo "PYPI_TOKEN not set"; \
		exit 1; \
	fi
	UV_PUBLISH_TOKEN="$$PYPI_TOKEN" uv publish

clean:
	@bash scripts/clean_pycache.sh .