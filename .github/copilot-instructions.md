# tools-iuc Repository Instructions

This repository hosts Galaxy tool wrappers maintained by the Intergalactic
Utilities Commission (IUC).

## Galaxy tool development

Before creating or modifying anything under `tools/`, `data_managers/`,
`tool_collections/`, `suites/`, `deprecated/`, or `macros/`, load the
**tool-dev** skill (`.agents/skills/tool-dev/SKILL.md`). It is the primary
Galaxy-specific reference — do not reinvent its guidance:

- `SKILL.md` — wrapper patterns, macros, Cheetah gotchas, inputs/outputs,
  testing, `.shed.yml`, and the IUC PR review checklist
- `references/testing.md` — planemo testing, assertions, failure analysis
- `references/tool-placement.md` — where a new tool belongs

The skill defers to the [IUC standards](https://galaxy-iuc-standards.readthedocs.io/)
for the rules. The target repository's conventions (sibling tools, this CI
setup) take precedence — when you deviate from the skill, say so and name the
sibling tool or convention you followed.

## Repository facts

- One flat directory per suite: `tools/<tool_name>/` holds all wrappers,
  `macros.xml`, `.shed.yml`, and `test-data/`.
- Macros shared across suites live in `macros/` at the repo root.
- CI runs `planemo lint`, `planemo test`, `flake8` (Python) and `styler` (R);
  see `.github/workflows/`.
- Read [CONTRIBUTING.md](../../CONTRIBUTING.md) and
  [docs/guide_for_reviewers.md](../../docs/guide_for_reviewers.md).
- If AI assisted a contribution, identify the AI harness and model in the PR
  description.

## Attribution

The `tool-dev` skill is vendored from
[galaxyproject/galaxy-skills](https://github.com/galaxyproject/galaxy-skills)
(`dev-skills/tool-dev`), MIT licensed © Galaxy Skills Contributors. Do not
edit the vendored files in this repository; propose changes upstream and
re-sync. `skills-lock.json` records sync state and is intentionally
git-ignored.
