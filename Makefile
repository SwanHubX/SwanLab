AGENTS_SKILL_DIR := .agents/skills
SKILLS_DIR := docs/skills

.PHONY: init sync format proto unit bench clean build publish backport link-skills unlink-skills relink-skills core-lint core-fmt core-test core-build core-tidy

# ----------------------------------
# SKILL (docs/skills)
# ----------------------------------

link-skills:
	@if [ ! -d "$(SKILLS_DIR)" ]; then \
		echo "Missing skills source: $(SKILLS_DIR)"; \
		exit 1; \
	fi
	@mkdir -p "$(AGENTS_SKILL_DIR)"
	@for skill_dir in $(SKILLS_DIR)/*/; do \
		name=$$(basename "$$skill_dir"); \
		link="$(AGENTS_SKILL_DIR)/$$name"; \
		if [ -e "$$link" ] && [ ! -L "$$link" ]; then \
			echo "Refusing to replace non-symlink: $$link"; \
			exit 1; \
		fi; \
		ln -sfn "../../$(SKILLS_DIR)/$$name" "$$link"; \
		echo "Linked $$name"; \
	done

unlink-skills:
	@for link in $(AGENTS_SKILL_DIR)/*; do \
		if [ -L "$$link" ]; then \
			rm "$$link"; \
		fi; \
	done

relink-skills: unlink-skills link-skills

# ----------------------------------
# Protobuf generation (scripts/generate_protos.py)
# ----------------------------------

proto:
	uv run scripts/generate_protos.py

# ----------------------------------
# Backport (scripts/backport.sh)
# ----------------------------------

backport:
	@if [ -z "$(PR)" ] || [ -z "$(BRANCH)" ]; then \
		echo "Usage: make backport PR=<pr-number> BRANCH=<release-branch>"; \
		echo "Example: make backport PR=1743 BRANCH=release/v0.8"; \
		exit 1; \
	fi
	bash scripts/backport.sh PR=$(PR) BRANCH=$(BRANCH)

clean:
	@bash scripts/clean_pycache.sh .

# ----------------------------------
# Python package (swanlab/)
# ----------------------------------

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

# ----------------------------------
# Go core (core/)
# ----------------------------------

core-lint:
	cd core && golangci-lint run

core-fmt:
	cd core && golangci-lint fmt

core-test:
	cd core && go test $$(go list ./... | grep -v /proto/)

core-build:
	cd core && go build ./...

core-tidy:
	cd core && go mod tidy
