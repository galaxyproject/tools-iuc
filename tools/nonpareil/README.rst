Nonpareil memory usage
=======================

By default, this tool is configured to limit the memory consumption to 1024 Mib.

Ideally this value should be larger than the sequences to analyze (discarding non-sequence elements like headers or quality). This is particularly important when running in multiple cores. This value is approximated. Maximum value in this version: 4194303.

You can set the NONPAREIL_MAX_MEMORY environmental variable in the destination section of the job_conf.xml file:

```
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
            <env id="NONPAREIL_MAX_MEMORY">1024</env>
        </destination>
    </destinations>
</job_conf>
```



