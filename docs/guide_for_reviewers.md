# New Tool Review Guide

This document describes a checklist suitable as a guide for reviewers of pull requests against the IUC's tools repo. It is a guide only and reviewer discretion is still advised.

**This document is a work in progress!**

The comprehensive list is aimed at new tools. Obviously for tool updates just use the appropriate section of the checklist on the PR diffs.

This checklist is based on the IUC's [Best Practices](https://galaxy-iuc-standards.readthedocs.io/en/latest/index.html) document.

## Repository

* [ ] Is this tool appropriate for the IUC repo (see [CONTRIBUTING.md](https://github.com/galaxyproject/tools-iuc/blob/master/CONTRIBUTING.md))
* [ ] Does the tool already exist in the toolshed?
    * [ ] Is this tool warranted?
    * [ ] Does the IUC account have access to the current repo?
* [ ] If the repository contains more than one tool, should it be separated or made a tool collection?
* [ ] Is there a `.shed.yml` file?
* [ ] Is there a tool `.xml` file?
* [ ] No `tool_dependencies.xml` file. (Has been deprecated and is no longer "Best Practice")


## Files

### .shed.yml

* [ ] Is there a correctly formatted `.shed.yml` file?
* [ ] Are the categories appropriate?
* [ ] Does the name match the `.xml` file or in the case of multiple tools the folder name?
* [ ] Is the owner set to `iuc`? If not, is the owner set to a current wrapper?
* [ ] Is there a `description`?
* [ ] Are there `homepage_url` & `remote_repository_url` fields and do they point somewhere sensible?

### tool.xml

**Linting**

* [ ] Does the tool pass planemo/travis lint with no warnings or errors?
* [ ] 4 spaces indentation?

**Order of XML Elements**

* [ ] Are the xml elements in the suggested order as in [Best Practices Coding Style](http://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/tool_xml.html#coding-style)?

**&lt;tool&gt; (name and id etc.)**

* [ ] id and name sensible and not previously used?
* [ ] Version follow [PEP 440](https://www.python.org/dev/peps/pep-0440/) with `+galaxyN`?
    - [ ] Alphanumeric only and no `-`
* [ ] Is there a `@TOOL_VERSION@` macro token used? (Should there be?)
* [ ] If there is a `profile`, is it appropriate?

**&lt;description&gt;**

* [ ] Is there a description tag?
* [ ] Is the description of suitable length?

**&lt;macros&gt;**

* [ ] If there is more than one tool present (tool collection), is there a `macros.xml` file?
* [ ] Is it appropriate?

**&lt;edam_topics&gt; & &lt;edam_operations&gt;**

* Link to the [EDAM browser](https://bioportal.bioontology.org/ontologies/EDAM?p=classes)

**&lt;[parallelism]&gt;**

**&lt;requirements&gt;**

* [ ] Are there corresponding conda packages in best practice channels?
* [ ] Are they versioned correctly with `@TOOL_VERSION@`? (or multiple packages/docker containers with correctly described versions)
* [ ] If using R/Perl/Python/Ruby packages etc versions specified in correctly formatted `*_environment` tags?

**&lt;~~code~~&gt;**

* [ ] this tagset has been deprecated and should not be used.  [galaxyproject/galaxy#2712](https://github.com/galaxyproject/galaxy/issues/2712)

**&lt;stdio&gt;**

* [ ] Is `<stdio>` tagset warranted (or should it use `detect_errors` in the `<command>` tag)

**&lt;version_command&gt;**

* [ ] Is there a version command?
* [ ] `<![CDATA[ ... ]]>` tags?

**&lt;command&gt;**

* [ ] `<![CDATA[ ... ]]>` tags?
* [ ] Is there `detect_errors` in the `<command ... >` tag? Or is it handled by a more comprehensive `<stdio>` tagset?
* [ ] No `interpreter` field in `<command ... >` tag - This is deprecated.
* [ ] Text parameters, input and output files `'single quoted'`?
* [ ] Is the Cheetah indented and readable?
* [ ] Are multiple commands joined with `&&`?
* [ ] Are any extra temporary files (such as indices etc.) created in the CWD?

**&lt;environment_variables&gt;**



**&lt;configfiles&gt;**

**&lt;inputs&gt; and &lt;parameters&gt;**

*General*
* [ ] Order of parameter attributes as per [Best Practices Coding Style](http://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/tool_xml.html#coding-style)
* [ ] Do the `argument` attributes include the long form of the underlying tool parameters?

*Boolean*
* [ ] Are the `truevalue` and `falsevalue` set with the underlying tool parameter?

*Dynamic Options*
* [ ] Do conditional parameters use a `select` and not a `boolean`
* [ ] Are the advanced parameters of the tool hidden with appropriate `<section>` tags?

**&lt;request_param_translation&gt;**

**&lt;outputs&gt;**

* [ ] Are optional output files selected using `<filter>` tags (with corresponding `<select>` tags in the `<inputs>` section)

**&lt;tests&gt;**

* [ ] Is most of the functionality of the tool tested?
* [ ] Are the test datasets included in the `test-data` directory?
* [ ] Is the test data of suitable size? (i.e. small..)
* [ ] Do the tests pass?
* [ ] Is the output filtering tested using the `expect_num_outputs` attribute?
* [ ] Are there unused files in `test-data/`
* [ ] In the case where the tool uses built in reference data, is there a `tool-data/tool_data_table_conf.xml.test` file?

**&lt;help&gt;**

* [ ] Is it book-ended with `<![CDATA[ ... ]]>` tags?
* [ ] Is it correctly formatted [restructuredText](http://docutils.sourceforge.net/docs/ref/rst/restructuredtext.html)?
* [ ] Are any images in the `./static/images` directory?

**&lt;citations&gt;**

* [ ] Is there a citation
    - [ ] Is it in bibtex or doi format? (doi preferred)

### Data Tables

For tools that use the built in reference data and indices:

* [ ] Is there a `tool-data/*.loc.sample` file?
* [ ] Is there a `tool-data/tool_data_table_conf.xml.sample` file?
