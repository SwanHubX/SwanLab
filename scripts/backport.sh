#!/usr/bin/env bash
# Backport a merged PR to a release branch:
# cherry-pick its merge commit onto origin/<BRANCH> and open a backport PR.
#
# Usage:
#   make backport PR=<pr-number> BRANCH=<release-branch>
#   # equivalent: bash scripts/backport.sh PR=1743 BRANCH=release/v0.8
set -euo pipefail

usage() {
  echo "Usage: bash scripts/backport.sh PR=<pr-number> BRANCH=<release-branch>"
  echo "Example: bash scripts/backport.sh PR=1743 BRANCH=release/v0.8"
  exit 1
}

PR=""
BRANCH=""
for arg in "$@"; do
  case "$arg" in
    PR=*) PR="${arg#PR=}" ;;
    BRANCH=*) BRANCH="${arg#BRANCH=}" ;;
    *) echo "error: unknown argument '$arg'"; usage ;;
  esac
done

[[ "$PR" =~ ^[0-9]+$ ]] || { echo "error: PR must be a PR number"; usage; }
[ -n "$BRANCH" ] || { echo "error: BRANCH is required"; usage; }
command -v gh >/dev/null 2>&1 || { echo "error: gh CLI not found, install it from https://cli.github.com/"; exit 1; }

# cherry-pick 不能与未提交的改动混在一起，要求工作区干净
[ -z "$(git status --porcelain)" ] || { echo "error: working tree not clean, commit or stash first"; exit 1; }

STATE=$(gh pr view "$PR" --json state --jq .state)
[ "$STATE" = "MERGED" ] || { echo "error: PR #$PR is not merged yet (state: $STATE)"; exit 1; }

SHA=$(gh pr view "$PR" --json mergeCommit --jq .mergeCommit.oid)
[ -n "$SHA" ] && [ "$SHA" != "null" ] || { echo "error: no merge commit found for PR #$PR"; exit 1; }
TITLE=$(gh pr view "$PR" --json title --jq .title)

git fetch origin "$BRANCH"

# squash 合入的提交 message 带有 "(#<PR>)"；目标分支上已存在则说明此前合入或 backport 过，拦截重复 backport
if git log -1 --format=%H --fixed-strings --grep="(#$PR)" "origin/$BRANCH" | grep -q .; then
  echo "error: origin/$BRANCH already contains a commit referencing (#$PR), duplicate backport?"
  exit 1
fi

BACKPORT_BRANCH="backport/$BRANCH/pr-$PR"
if git show-ref --verify --quiet "refs/heads/$BACKPORT_BRANCH"; then
  echo "error: local branch '$BACKPORT_BRANCH' already exists, delete it first: git branch -D $BACKPORT_BRANCH"
  exit 1
fi

echo "==> Backporting #$PR ($SHA) to $BRANCH"
git switch --create "$BACKPORT_BRANCH" "origin/$BRANCH"

# squash 合入产生单父提交；真正的 merge commit 需要 -m 1 才能取到 PR 侧改动
if git rev-parse --verify --quiet "$SHA^2" >/dev/null 2>&1; then
  pick=(cherry-pick -m 1 "$SHA")
else
  pick=(cherry-pick "$SHA")
fi

# 失败时保留冲突现场，由使用者手动解决后续步骤
if ! git "${pick[@]}"; then
  cat <<EOF

error: cherry-pick failed, repository left on '$BACKPORT_BRANCH'.

Resolve manually (on conflict):
  1. fix conflicts, then: git add <files> && git cherry-pick --continue
  2. git push -u origin $BACKPORT_BRANCH
  3. gh pr create --base $BRANCH --head $BACKPORT_BRANCH \\
       --title "[backport] $TITLE" --body "Backport of #$PR"

Or abort everything:
  git cherry-pick --abort; git branch -D $BACKPORT_BRANCH
EOF
  exit 1
fi

git push -u origin "$BACKPORT_BRANCH"
gh pr create --base "$BRANCH" --head "$BACKPORT_BRANCH" \
  --title "[backport] $TITLE" \
  --body "Backport of #$PR"

echo "==> Done. You are now on '$BACKPORT_BRANCH'."
