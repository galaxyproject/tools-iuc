---
applyTo: "tools/**, data_managers/**, tool_collections/**, suites/**, deprecated/**, macros/**"
---

# Galaxy Tool Wrapper Review Guidance

Load the `tool-dev` skill (`.agents/skills/tool-dev/SKILL.md`) for detailed
conventions. These instructions cover review priorities and, critically, what
NOT to flag. Do not approve changes automatically and do not suppress genuine
findings — the goal is fewer false positives, not fewer comments.

## Real correctness issues — always report

- Cheetah truthiness on booleans/selects: `#if $param:` is always true for a
  select whose value is the string `False`. Requires
  `#if str($param) == 'True':` (exact casing match).
- Falsy zero: `#if $int_param:` skips the flag when the value is 0. Use
  `#if str($int_param):`.
- Cheetah (`#if`, `#for`) inside an `<xml>` macro — does not work; must be a
  `<token>`.
- Loop-variable shadowing in `<repeat>` blocks: params must be accessed
  through the loop variable (`$item.x`), not the outer name.
- Output `<filter>` checking the wrong value (copy-paste between outputs).
- Missing `.strip()` on optional text params (whitespace-only input treated
  as present).
- Unanchored `discover_datasets` regex patterns (must end with `$`).
- Orphaned parameters: in `<inputs>` but not used in `<command>` (or vice
  versa).
- `expect_num_outputs` not matching the outputs actually produced by the test
  case (filtered-out outputs don't count).
- Version bumps: `@TOOL_VERSION@` changed but `@VERSION_SUFFIX@` not reset to
  `0`; version-token-only diffs should be checked against upstream release
  notes before being treated as mechanical.
- Review the whole tool directory on version bumps — shared macros affect
  sibling wrappers absent from the diff.

## Existing conventions — do NOT flag these as problems

- `detect_errors="exit_code"` on script/configfile-based wrappers (Python/R)
  is intentional; `aggressive` false-fails on benign library warnings.
- Older tools using RST help, bare `name=` without `argument=`, or older
  profile versions predate current guidance — only raise these for new tools
  or when the touched lines are affected. Suggest modernization at most as a
  non-blocking note.
- `.lint_skip` / `.tt_skip` entries are deliberate; don't suggest deleting
  them without addressing the underlying issue.
- Test data over 1 MB may predate the guideline; flag only for new test data.
- Single-tool directories without `macros.xml` are fine.
- Suite `remote_repository_url` must use `tree/main` (default branch is
  `main`, not `master`).

## Review style for this repository

- Anchor findings to the IUC standards or the `tool-dev` skill; explain the
  why in one sentence.
- Don't demand fixes for style that `planemo format` / `planemo lint` would
  settle — suggest running them instead.
- Distinguish clearly between required changes and optional suggestions, as
  docs/guide_for_reviewers.md requires.
- If a change is tool XML, note whether `planemo lint` and `planemo test`
  results were provided; missing test coverage for a new code path is a
  legitimate finding.
