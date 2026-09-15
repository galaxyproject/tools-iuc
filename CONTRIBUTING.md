# Contributing

This document describes how to contribute to this repository.
Pull requests containing new tools or bug fixes, updates, and extensions to the existing
tools and tool suites in this repository will be considered for inclusion.

## How to Contribute

* Make sure you have a [GitHub account](https://github.com/signup/free)
* Make sure you have git [installed](https://docs.github.com/en/get-started/git-basics/set-up-git)
* Fork the repository on [GitHub](https://github.com/galaxyproject/tools-iuc/fork)
* Make the desired modifications - consider using a [feature branch](https://github.com/Kunena/Kunena-Forum/wiki/Create-a-new-branch-with-git-and-manage-branches).
* Make sure you have added the necessary tests for your changes and they pass.
* Make sure submitted tools meet IUC [Best Practices](https://galaxy-iuc-standards.readthedocs.io/)
* Open a [pull request](https://docs.github.com/en/pull-requests/reference/pull-requests)
  with these changes.


## What to Contribute

* Wrappers for new [OSI-approved licensed](https://opensource.org/licenses) tools
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

### Eligibility Criteria for New Tools

In general, contributions of new tools from all scientific domains are welcome if they are in line with a few guiding principles:

* New tool contributions should, as far as possible, wrap software that is developed and packaged upstream.

  It is possible to run custom scripts from a tool wrapper, but these scripts should normally not be at the core of your contribution.
  This repo is intended as a home for Galaxy tool *wrappers*, not for the wrapped tools themselves.

* New tool contributions should not promote trivial forks of maintained upstream tools.

  It is not our aim to collect highly similar wrappers for all kinds of flavors of the same software.
  In particular, when one of these flavors is clearly identifiable as the original software and is still actively maintained,
  any fork or rewrite not approved by the original author will have to demonstrate *significant* improvements over the original to justifiy its inclusion here.

The above criteria are not rigid rules, but exist to guide reviewers when they are in doubt about the value of a proposed new contribution.

If considering a new tool contribution, which might fail to meet one of the criteria, consider opening an issue to gather feedback from IUC members first.


## What Not to Contribute

* Things already wrapped and currently maintained by other users
* Any new tools clearly not meeting the eligibility criteria above
* New datatypes. These should be added directly to the [Galaxy codebase](https://github.com/galaxyproject/galaxy).


## How to Get Your Contribution Accepted

### Adherence to IUC Standards and Best Practices

For your contribution to be accepted, your code has to adhere to the [IUC standards and best practices](https://galaxy-iuc-standards.readthedocs.io/).

### Tests

All contributed tools should include test cases.
They need not necessarily cover all uses of the program, but should ensure that it is generally working.
The Galaxy Hub has a [page](https://galaxyproject.org/admin/tools/writing-tests/) on writing tests.

The IUC strongly recommends testing with [planemo](https://github.com/galaxyproject/planemo/), which provides simple subcommands for checking and testing tool wrappers.
In particular, you may want to use [planemo lint](https://planemo.readthedocs.io/en/latest/commands.html#lint-command), [planemo test](https://planemo.readthedocs.io/en/latest/commands.html#test-command) and [planemo serve](https://planemo.readthedocs.io/en/latest/commands.html#serve-command).

The IUC's continuous integration pipeline will run `planemo lint` and `planemo test` on any tool wrappers in your contribution, too, so it is a good idea to use the same tooling locally.
CI runs will also include linting of any Python and/or R scripts with `flake8` and `styler`, respectively.

Your contribution will have to pass all linting and all tests before it can get accepted.

### Other Requirements

* If there's a relevant paper for the tool, it should be cited in a [citation](https://docs.galaxyproject.org/en/latest/dev/schema.html#tool-citations) block
* The tool must be licensed to allow use by anyone. The OSI maintains a [list of appropriate licenses](https://opensource.org/licenses)
* At least one approving review by an IUC member
