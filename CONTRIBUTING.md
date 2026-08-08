# Contributing

This document describes how to contribute to this repository. Pull
requests containing bug fixes, updates, and extensions to the existing
tools and tool suites in this repository will be considered for
inclusion.

## How to Contribute

* Make sure you have a [GitHub account](https://github.com/signup/free)
* Make sure you have git [installed](https://docs.github.com/en/get-started/git-basics/set-up-git)
* Fork the repository on [GitHub](https://github.com/galaxyproject/tools-iuc/fork)
* Make the desired modifications - consider using a [feature branch](https://github.com/Kunena/Kunena-Forum/wiki/Create-a-new-branch-with-git-and-manage-branches).
* Make sure you have added the necessary tests for your changes and they pass.
* Make sure submitted tools meet IUC [Best Practices](https://galaxy-iuc-standards.readthedocs.io/)
* Open a [pull request](https://docs.github.com/en/pull-requests/reference/pull-requests)
  with these changes.

## What to contribute

* Wrappers for new [OSI-approved licensed](https://opensource.org/licenses) tools
* Visualization Plugins
* Updates for tools
* Enhancements for tools (e.g. supporting new parameters)
* Bug fixes
* Documentation improvements
* New test cases

### Abandoned Tools

* For tools of general interest, the IUC is usually willing to adopt tools that
  you (the developer) are abandoning.
* If there are tools that you find useful but seem to be abandoned and not
  updated, you're welcome to create an issue recommending that the IUC adopt
  that tool.

## What not to contribute

* Things already wrapped and currently maintained by other users
* Wrappers without tests
* New datatypes. These should be added directly to the [Galaxy codebase](https://github.com/galaxyproject/galaxy).

## Tests

All contributed tools should include test cases. They need not
necessarily cover all uses of the program, but should ensure that it is
generally working. The Galaxy Hub has a
[page](https://galaxyproject.org/admin/tools/writing-tests/) on writing
tests.

The IUC strongly recommends testing with [planemo](https://github.com/galaxyproject/planemo/), which provides a simple command line utility for testing functionality

```console
$ planemo test --install_galaxy my_tool.xml
```

## Requirements for Contributions

Before a PR will be accepted, the IUC has [some requirements](https://galaxy-iuc-standards.readthedocs.io/) on the
submitted code (which we will be happy to help you achieve if you need the
assistance).

* Continuous integration tests must pass: 
    * The tests must pass (`planemo test --install_galaxy my_tool.xml`)
    * The tools must pass linting by planemo (`planemo lint my_tool.xml`)
    * Any Python or R script must pass linting with `flake8` and `styler`, respectively
* If there's a relevant paper for the tool, it should be cited in a [citation](https://docs.galaxyproject.org/en/latest/dev/schema.html#tool-citations) block
* The tool must be licensed to allow use by anyone. The OSI maintains a [list of appropriate licenses](https://opensource.org/licenses)
* At least one approving review by an IUC member
