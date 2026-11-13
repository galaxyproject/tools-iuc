# mothur (archived wrappers)

Archived copies of mothur wrappers whose commands mothur no longer provides. The rest of
the mothur suite lives in `tools/mothur`.

## parse.list

`parse.list` was removed from mothur in version 1.43.0 (1.42.1 is the last release that
still knows the command); the wrapper is pinned to `mothur=1.39.5`, the version it was
written and last tested against.

The wrapper and the `macros.xml` next to it are frozen at their state before the 1.48.6
upgrade of the suite — `macros.xml` is a private copy of the suite macros from that time so
that the archived wrapper keeps its own dependency pin, its own version token and its own
test expectations.

This directory is not built, linted or tested by CI: the pull-request workflows ignore
`deprecated/**`. To run the tests locally, put mothur 1.39.5 on the `PATH` and call

    planemo test --no_dependency_resolution deprecated/tools/mothur/parse.list.xml
