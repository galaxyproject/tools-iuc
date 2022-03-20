DISCO memory usage
==================


By default, this tool is configured to limit the memory consumption to 4G.
You might need to lower this limit if the machine(s) executing the jobs have less memory available.
If you have a lot of reads to assemble and a machine with enough memory, you can increase it.

In both cases, you can set the DISCO_MAX_MEMORY environmental variable in the destination section of the job_conf.xml file::

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
                <env id="DISCO_MAX_MEMORY">4</env>
            </destination>
        </destinations>
    </job_conf>
