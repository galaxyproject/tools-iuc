
=================
IMPORTANT NOTE REGARDING SYSTEM CONFIGURATION
=================

All of the Galaxy wrappers contained herein call the respective bbtools' shell wrapper, which calls the underlying java-based tool. Unlike a C-based program, java will grab a pre-determined amount of memory at the very beginning of the execution.

Some of the algorithms (e.g. bbnorm) utilise a hash table, and potential collusions can decrease the numeric accuracy of the output. This problem is expected to become more pronounced if the fraction of the memory occupied w.r.t. allocated memory becomes high, i.e. when the available memory is low and/or the input file is big. If the tool generates a warning to stderr, and will be caught by the Galaxy wrapper resulting in a failed job. However, `count min sketch <https://en.wikipedia.org/wiki/Count%E2%80%93min_sketch>` does not run out of memory, this is a gradual effect, and will NOT trigger a fatal error unless the load reaches this critically high level. You can read more about the implications of this at the `BBtools manual <https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbnorm-guide/>`

If you are administering a heteregenous computing environment with multiple nodes of very different quantities of physically available RAM, it is recommended to define a global cap on the RAM to be used to avoid introducing run-to-run bias by exporting an environmental variable, by something like:
export _JAVA_OPTIONS="-Xmx2048m -Xms256m"

The tool currently considers the following limits, in the given priority order:
1) _JAVA_OPTIONS
2) JAVA_TOOL_OPTIONS
3) GALAXY_MEMORY_MB
4) 4 GB

