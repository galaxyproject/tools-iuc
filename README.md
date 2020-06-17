[![Build Status](https://travis-ci.org/galaxyproject/tools-iuc.svg?branch=master)](https://travis-ci.org/galaxyproject/tools-iuc)
[![Gitter](https://badges.gitter.im/galaxyproject/tools-iuc.svg)](https://gitter.im/galaxy-iuc/iuc?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge) 

Galaxy Tools maintained by IUC
==============================

![](iuc_logo.png)

This repo contains a subset of Galaxy repositories used in the Tool Shed (https://toolshed.g2.bx.psu.edu/).

These repositories are maintained and developed by the [Intergalactic Utilities Commission](https://galaxyproject.org/iuc/) ([Github Org](https://github.com/galaxy-iuc/))

Pull Requests with dependencies specified as conda-package will be automatically tested and verified using [planemo](https://github.com/galaxyproject/planemo). If the tests pass and the pull request is merged, the tool will be automatically uploaded to the [Test](http://testtoolshed.g2.bx.psu.edu/)- and [Main Tool Shed](http://toolshed.g2.bx.psu.edu/).

Please note, if you don’t want to run the tests or the automatic upload, add `[ci skip]` to the git commit message.
Commits that have [ci skip] anywhere in the commit messages are ignored by Travis CI.


Other repositories with Galaxy tools:
 * [Björn Grüning repo](https://github.com/bgruening/galaxytools)
 * [Galaxy devteam repo](https://github.com/galaxyproject/tools-devteam)
 * Peter Cock's repos:
   * [blast repo](https://github.com/peterjc/galaxy_blast)
   * [pico repo](https://github.com/peterjc/pico_galaxy)
   * [mira repo](https://github.com/peterjc/galaxy_mira)
 * [ENCODE tools](https://github.com/modENCODE-DCC/Galaxy)
 * [Biopython repo](https://github.com/biopython/galaxy_packages)
 * [Galaxy Proteomics repo](https://github.com/galaxyproteomics/tools-galaxyp)
 * [Colibread Galaxy Tools](https://github.com/genouest/tools-colibread)
 * [Greg von Kuster's repo](https://github.com/gregvonkuster/galaxy-csg)
 * [EI repo](https://github.com/TGAC/earlham-galaxytools)
 * [AAFC-MBB Canada repo](https://github.com/AAFC-MBB/Galaxy/tree/master/wrappers)
 * [Mark Einon's repo](https://gitlab.com/einonm/galaxy-tools)
 * [National Microbiology Laboratory's repo](https://github.com/phac-nml/galaxy_tools)

 
