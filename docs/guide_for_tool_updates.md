# Galaxy Tool Update Guide

This guide describes how to update an existing Galaxy tool wrapper or tool
suite in this repository to a newer upstream software version. It complements
the [new tool review guide](guide_for_reviewers.md): this document describes the
update workflow, while the reviewer guide provides the corresponding review
checks.

An update is more than changing a version token. Upstream releases may change
command-line arguments, defaults, output formats, exit behavior, dependencies,
or licensing. Treat the version change, wrapper adaptation, tests, and review
evidence as one change.

## 1. Establish the update scope

Before editing files:

1. Identify the requested upstream version and the wrapper directory.
2. Find every tool XML file and macro file in the tool suite.
3. Identify the main Conda requirement and how the suite defines
   `@TOOL_VERSION@` and `@VERSION_SUFFIX@`.
4. Confirm that the requested package version exists in the best-practice Conda
   channels and that a suitable BioContainer is available when applicable.
5. Read the upstream release notes for every version between the current and
   requested versions.

Keep an explicitly requested target version. Do not silently update to a newer
release merely because it is the newest version available from Conda.

Stop and report the blocker when the requested package is unavailable, the
upstream release is incompatible with the wrapper's supported interface, or the
required behavior cannot be verified.

## 2. Perform the mechanical update

[Planemo's `autoupdate` command](https://planemo.readthedocs.io/en/latest/autoupdate.html)
can update wrappers that follow the standard version-token layout. Preview its
changes first:

```sh
planemo autoupdate --dry-run tools/<tool>/<wrapper>.xml
```

Then run it without `--dry-run` when the proposed target is correct:

```sh
planemo autoupdate tools/<tool>/<wrapper>.xml
```

For a tool suite, update each intended wrapper or use `--recursive` after
checking that every wrapper in the directory belongs in the same update:

```sh
planemo autoupdate --dry-run --recursive tools/<tool>
planemo autoupdate --recursive tools/<tool>
```

If the requested version is not the newest version selected by `autoupdate`,
edit the version tokens and requirements deliberately instead.

After an upstream software update:

- Set `@TOOL_VERSION@` to the upstream/package version used by the wrapper.
- Keep the main requirement synchronized with `@TOOL_VERSION@`.
- Reset `@VERSION_SUFFIX@` to `0`.
- Ensure every affected tool's `version` attribute resolves to the intended
  `<upstream-version>+galaxy0` value.

For a wrapper-only change with no upstream software update, keep
`@TOOL_VERSION@` unchanged and increment `@VERSION_SUFFIX@` instead.

## 3. Adapt the wrapper to the upstream release

Compare the updated software's help and release notes with the wrapper. Check at
least:

- Added, removed, renamed, or behavior-changing command-line arguments.
- Changed defaults, accepted values, validation constraints, or required
  combinations.
- Changed input and output formats, filenames, headers, ordering, or precision.
- Changed exit codes, warnings, logging, or error-detection behavior.
- Added or changed runtime dependencies, data files, environment variables, or
  licenses.
- Help text, citations, and version reporting that mention old behavior.

Make the smallest coherent update. Avoid unrelated cleanup unless it is needed
for the new release or removes a concrete lint or maintenance problem. In
particular, changing the Galaxy tool `profile` can change runtime semantics; do
so intentionally and test the consequences.

## 4. Update tests deliberately

Run the existing tests before changing their expectations when practical, then
use failures to distinguish upstream behavior changes from wrapper defects.

- Add or adjust tests for new or changed behavior exposed by the wrapper.
- Preserve coverage of unchanged behavior.
- Prefer semantic assertions over broad tolerances or complete expected-output
  replacement.
- When output contains nondeterministic content, loosen only the unstable part
  and retain assertions for stable invariants.
- Do not accept regenerated test data without inspecting and explaining every
  meaningful difference.
- Keep test data small and remove files that are no longer used.

Use `planemo autoupdate --test` only when its selected target version and test
scope are correct. Treat `--update_test_data` as a proposal generator, not as
approval of the new outputs.

## 5. Validate the complete update

Run linting and tests against the complete tool directory so shared macros and
every affected wrapper are covered:

```sh
planemo lint --biocontainers tools/<tool>
planemo test --install_galaxy tools/<tool>
```

Run any additional formatter or language-specific checks required by files in
the directory. Review the final diff after testing and confirm that it contains
only intentional changes.

## 6. Prepare the pull request

In the pull request description, record:

- The old and new upstream versions.
- A link to the relevant upstream release notes.
- Important wrapper adaptations and changed defaults or outputs.
- Why test expectations or test data changed.
- The lint and test commands run, including any tests that could not be run.

Update every wrapper in a suite that shares the changed version unless a partial
update is intentional. The `skip-version-check` label is only for the documented
case where a suite intentionally updates a subset of its tools; it is not a way
to bypass an unexplained versioning failure.
