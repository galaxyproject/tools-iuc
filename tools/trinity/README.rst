Trinity memory usage
====================

As described on the official Trinity website
`FAQ <http://trinityrnaseq.sourceforge.net/trinity_faq.html#ques_comp_resources_required>`_,
trinity requires a large amount of memory to perform the assembly: "roughly
~1G of RAM per 1M reads to be assembled"

By default, this tool is configured to limit the memory consumption to 30G.
You might need to lower this limit if the machine(s) executing the jobs have less memory available.
If you have a lot of reads to assemble and a machine with enough memory, you can increase it.
In both cases, you can edit the TRINITY_MEM_OPTIONS in the file:

<tool_dependency_dir>/environment_settings/TRINITY_MEM_OPTIONS/iuc/trinity/<hash_string>/env.sh

to lower the maximum memory usage to 2G for example:

TRINITY_MEM_OPTIONS='--max_memory 2G --bflyHeapSpaceMax 2G'
