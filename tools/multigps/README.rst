MultiGPS wrapper for Galaxy
================================

* http://mahonylab.org/software/multigps/

MultiGPS performs significant EM optimization of binding events along the genome and across experimental
conditions, and it integrates motif-finding via MEME.  The tool loads all data into memory, so the potential
exists for time and memory intensive analyses if running over many conditions or large datasets.

Setting the memory allocation in Galaxy for this tool is handled using the <env id="_JAVA_OPTIONS"> tag for
a selected job runner in the job_conf.xml file.

License
-------

MIT

