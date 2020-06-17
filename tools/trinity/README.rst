Trinity memory usage
====================

As described on the official Trinity website
`FAQ <http://trinityrnaseq.sourceforge.net/trinity_faq.html#ques_comp_resources_required>`_,
trinity requires a large amount of memory to perform the assembly: "roughly
~1G of RAM per 1M reads to be assembled"

By default, this tool is configured to limit the memory consumption to 1G.
You might need to lower this limit if the machine(s) executing the jobs have less memory available.
If you have a lot of reads to assemble and a machine with enough memory, you can increase it.

This tool uses the GALAXY_MEMORY_MB environmental variable to limit the its memory usage. If you use a supported cluster type, you just need to limit the memory in the destination section of the job_conf.xml file to the amount of available memory (in MB), for example with a Slurm cluster::

    <destination id="cluster_high_mem" runner="slurm
        <param id="nativeSpecification">--mem=1024</param>
    </destination>

If you cluster type does not support memory limit detection, you can also set manually the GALAXY_MEMORY_MB environmental variable in the destination section of the job_conf.xml file to the amount of available memory (in MB)::

    <?xml version="1.0"?>
    <!-- A sample job config that explicitly configures job running the way it is configured by default (if there is no explicit config). -->
    <job_conf>
        <plugins>
            <plugin id="local" type="runner" load="galaxy.jobs.runners.local:LocalJobRunner" workers="4"/>
        </plugins>
        <handlers>
            <handler id="main"/>
        </handlers>
        <destinations>
            <destination id="local" runner="local">
                <env id="GALAXY_MEMORY_MB">1024</env>
            </destination>
        </destinations>
    </job_conf>


