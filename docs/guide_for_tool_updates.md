# Galaxy Tool Update Guide

This guide describes how to update an existing Galaxy tool wrapper or tool
suite in this repository to a newer upstream software version. It complements
the [new tool review guide](guide_for_reviewers.md): this document describes the
update workflow, while the reviewer guide provides the corresponding review
checks.

The [IUC Standards and Best
Practices](https://galaxy-iuc-standards.readthedocs.io/en/latest/) remain the
canonical source for general wrapper requirements. In particular, follow the
standards for [tool XML and
versioning](https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/tool_xml.html),
[dependencies](https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/package_xml.html),
and [ToolShed
readiness](https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/integration_checklist.html).
This guide only describes the additional work needed when updating an existing
tool in tools-iuc.

An update is more than changing a version token. Upstream releases may change
command-line arguments, defaults, output formats, exit behavior, dependencies,
or licensing. Treat the version change, wrapper adaptation, tests, and review
evidence as one change.

## 1. Establish the update scope

Before editing files:

1. Identify the target upstream version and the wrapper directory.
2. Find every tool XML file and macro file in the tool suite.
3. Determine which wrappers share the version and requirement macros and must
   move together.
4. Confirm that the target package version exists in the best-practice Conda
   channels. Do not require a BioContainer for the new version before merge;
   that image may only be built after the update lands on the main branch.
5. Read the upstream release notes for every version between the current and
   target versions.

Keep the pull request description and the diff in agreement about the target
version. If the target changes while the pull request is open, say so instead
of silently retargeting it.

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

If the target version is not the newest version selected by `autoupdate`,
edit the version tokens and requirements deliberately instead.

Apply the [IUC tool-version
rules](https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/tool_xml.html#tool-versions)
to every affected wrapper. Check the complete suite after automation runs;
shared macros can change wrappers that do not otherwise appear in the diff.

## 3. Adapt the wrapper to the upstream release

Compare the updated software's help and all intervening release notes with the
wrapper. Identify changes that affect the wrapped interface, such as command
line options and defaults, input or output formats, dependencies, or licensing.
Apply the corresponding IUC standards instead of restating them here.

Make the smallest coherent update. Avoid unrelated cleanup unless it is needed
for the new release or to meet a current standard. When updating the tool
`profile` to satisfy the [IUC profile
standard](https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/tool_xml.html#tool-profile),
test any resulting change in Galaxy semantics.

## 4. Update tests deliberately

Run the existing tests before changing their expectations when practical, then
use failures to distinguish upstream behavior changes from wrapper defects.

Follow the [IUC test
standards](https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/tool_xml.html#tests).
Add or adjust coverage for release changes exposed by the wrapper, and preserve
coverage of unchanged behavior. Inspect every regenerated test-data difference
and tie it to an upstream or wrapper behavior change.

Use `planemo autoupdate --test` only when its selected target version and test
scope are correct. Treat `--update_test_data` as a proposal generator, not as
approval of the new outputs.

## 5. Validate the complete update

Run linting and tests against the complete tool directory so shared macros and
every affected wrapper are covered:

```sh
planemo lint tools/<tool>
planemo test tools/<tool>
```

Use the Galaxy instance, profile, or other test-engine options appropriate for
your local environment. CI tests with `--biocontainers` (except for the paths
listed in `.tt_biocontainer_skip`), so use it locally too when you can and test
what CI tests. Pull request CI runs the repository lint and test jobs; it does
not require the new BioContainer to exist. The main-branch lint job adds the
BioContainer check after the update has merged.

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
