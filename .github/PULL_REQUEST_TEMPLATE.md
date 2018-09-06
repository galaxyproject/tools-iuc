FOR CONTRIBUTOR:
* [ ] - I have read the [CONTRIBUTING.md](https://github.com/galaxyproject/tools-iuc/blob/master/CONTRIBUTING.md) document and this tool is appropriate for the tools-iuc repo.
* [ ] - License permits unrestricted use (educational + commercial)
* [ ] - This PR adds a new tool or tool collection
* [ ] - This PR updates an existing tool or tool collection
* [ ] - This PR does something else (explain below)

FOR REVIEWER:
* [ ] .shed.yml file ok
    - [ ] Toolshed user `iuc` has access to associated toolshed repo(s)
* [ ] Indentation is correct (4 spaces)
* [ ] Tool version/build ok
* [ ] `<command/>`
  - [ ] Text parameters, input and output files `'single quoted'`
  - [ ] Use of `<![CDATA[ ... ]]>` tags
  - [ ] Parameters of type `text` or having `optional="true"` attribute are checked with `if str($param)` before being used
* [ ] Data parameters have a `format` attribute containing datatypes recognised by Galaxy
* [ ] Tests
  - [ ] Parameters are reasonably covered
  - [ ] Test files are appropriate
* [ ] Help
  - [ ] Valid restructuredText and uses `<![CDATA[ ... ]]>` tags
* [ ] Complies with other best practice in [Best Practices Doc](http://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/tool_xml.html)
