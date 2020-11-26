## Notes for setting up the SpyBOAT Galaxy tool

### Toolshed(s)

Register at:

- https://testtoolshed.g2.bx.psu.edu/
- https://toolshed.g2.bx.psu.edu/


### commands

First put toolshed account details into `.planemo.yml`

- planemo lint [tool.xml]
- planemo shed_init --name spyboat
- planemo shed_create --shed_target testtoolshed
- planemo shed_update --shed_target toolshed


### Hosting repo

- https://github.com/galaxyproject/tools-iuc

#### Guidelines

- https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/tool_xml.html