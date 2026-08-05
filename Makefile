SKILL_NAME := pr-review-lab
SKILL_SOURCE := docs/skills/$(SKILL_NAME)
AGENTS_SKILL_DIR := .agents/skills
AGENTS_SKILL_LINK := $(AGENTS_SKILL_DIR)/$(SKILL_NAME)
# symlink target relative to AGENTS_SKILL_DIR, so the link stays valid if the repo moves
SKILL_LINK_TARGET := ../../$(SKILL_SOURCE)

.PHONY:  init sync format proto unit bench clean build publish link-skills unlink-skills relink-skills

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

link-skills:
	@if [ ! -d "$(SKILL_SOURCE)" ]; then \
		echo "Missing skill source: $(SKILL_SOURCE)"; \
		exit 1; \
	fi
	@mkdir -p "$(AGENTS_SKILL_DIR)"
	@if [ -e "$(AGENTS_SKILL_LINK)" ] && [ ! -L "$(AGENTS_SKILL_LINK)" ]; then \
		echo "Refusing to replace non-symlink: $(AGENTS_SKILL_LINK)"; \
		exit 1; \
	fi
	@ln -sfn "$(SKILL_LINK_TARGET)" "$(AGENTS_SKILL_LINK)"

unlink-skills:
	@if [ -L "$(AGENTS_SKILL_LINK)" ]; then \
		rm "$(AGENTS_SKILL_LINK)"; \
	fi

relink-skills: unlink-skills link-skills
