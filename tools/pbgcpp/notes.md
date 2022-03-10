This is working

```bash

singularity exec \
    -B /media,${PWD},${TMPDIR},/tmp --nv -H $(mktemp -d) --pwd ${PWD} --containall --cleanenv --writable-tmpfs \
        docker://quay.io/biocontainers/genomicconsensus@sha256:de72299d4fb4f2bd25abdb0309527c0ad5b39e3e6b1216f76456324a642962ab \
    variantCaller \
        --numWorkers 32 \
        --referenceFilename test-data/All4mer.V2.01_Insert.fa \
        --outputFilename test-data/output.fa \
        test-data/out.aligned_subreads.bam
```

But I'm getting a strange error when I run planemo.

```
Job in error state.. tool_id: genomicconsensus_arrow, exit_code: 255, stderr: FATAL:   container creation failed: mount /tmp/tmpqd7x7e8m/job_working_directory/000/5->/tmp/tmpqd7x7e8m/job_working_directory/000/5 error: while mounting /tmp/tmpqd7x7e8m/job_working_directory/000/5: destination /tmp/tmpqd7x7e8m/job_working_directory/000/5 doesn't exist in container
```

I don't know how it's running Singularity but it's not working.

This is how Singularity gets run:

```bash
cd working; \
SINGULARITYENV_GALAXY_SLOTS=$GALAXY_SLOTS SINGULARITYENV_GALAXY_MEMORY_MB=$GALAXY_MEMORY_MB SINGULARITYENV_GALAXY_MEMORY_MB_PER_SLOT=$GALAXY_MEMORY_MB_PER_SLOT SINGULARITYENV_HOME=$HOME SINGULARITYENV__GALAXY_JOB_HOME_DIR=$_GALAXY_JOB_HOME_DIR SINGULARITYENV__GALAXY_JOB_TMP_DIR=$_GALAXY_JOB_TMP_DIR SINGULARITYENV_TMPDIR=$TMPDIR SINGULARITYENV_TMP=$TMP SINGULARITYENV_TEMP=$TEMP \
        singularity -s exec \
            -B /home/tom/.planemo/galaxy:/home/tom/.planemo/galaxy \
            -B /home/tom/.planemo/planemo_tmp_07y8xlaa:/home/tom/.planemo/planemo_tmp_07y8xlaa \
            -B /tmp/tmpqd7x7e8m/job_working_directory/000/5:/tmp/tmpqd7x7e8m/job_working_directory/000/5 \
            -B /tmp/tmpqd7x7e8m/job_working_directory/000/5/outputs:/tmp/tmpqd7x7e8m/job_working_directory/000/5/outputs \
            -B "$_GALAXY_JOB_TMP_DIR:$_GALAXY_JOB_TMP_DIR" \
            -B "$_GALAXY_JOB_HOME_DIR:$_GALAXY_JOB_HOME_DIR" \
            -B /tmp/tmpqd7x7e8m/job_working_directory/000/5/working:/tmp/tmpqd7x7e8m/job_working_directory/000/5/working \
            -B /tmp/tmpqd7x7e8m/files:/tmp/tmpqd7x7e8m/files \
            -B /home/tom/.planemo/planemo_tmp_07y8xlaa/test-data:/home/tom/.planemo/planemo_tmp_07y8xlaa/test-data \
            -B /home/tom/.planemo/galaxy/tool-data:/home/tom/.planemo/galaxy/tool-data \
            -B /home/tom/.planemo/galaxy/tool-data:/home/tom/.planemo/galaxy/tool-data \
            --home $HOME:$HOME \
            docker://quay.io/biocontainers/genomicconsensus@sha256:de72299d4fb4f2bd25abdb0309527c0ad5b39e3e6b1216f76456324a642962ab \
            /bin/sh \
            /tmp/tmpqd7x7e8m/job_working_directory/000/5/tool_script.sh \
            > ../outputs/tool_stdout \
            2> ../outputs/tool_stderr; \
            return_code=$?; 


```


```bash
cd working; \
SINGULARITYENV_GALAXY_SLOTS=$GALAXY_SLOTS SINGULARITYENV_GALAXY_MEMORY_MB=$GALAXY_MEMORY_MB SINGULARITYENV_GALAXY_MEMORY_MB_PER_SLOT=$GALAXY_MEMORY_MB_PER_SLOT SINGULARITYENV_HOME=$HOME SINGULARITYENV__GALAXY_JOB_HOME_DIR=$_GALAXY_JOB_HOME_DIR SINGULARITYENV__GALAXY_JOB_TMP_DIR=$_GALAXY_JOB_TMP_DIR SINGULARITYENV_TMPDIR=$TMPDIR SINGULARITYENV_TMP=$TMP SINGULARITYENV_TEMP=$TEMP \
        singularity -s exec \
            -B /home/tom/.planemo/galaxy:/home/tom/.planemo/galaxy \
            -B /tmp \
            -B /tmp/tmpqd7x7e8m/job_working_directory/000/5:/tmp/tmpqd7x7e8m/job_working_directory/000/5 \
            -B /tmp/tmpqd7x7e8m/files:/tmp/tmpqd7x7e8m/files \
            -B /home/tom/.planemo/galaxy/tool-data:/home/tom/.planemo/galaxy/tool-data \
            -B /home/tom/.planemo/galaxy/tool-data:/home/tom/.planemo/galaxy/tool-data \
            --home $HOME:$HOME \
            docker://quay.io/biocontainers/genomicconsensus@sha256:de72299d4fb4f2bd25abdb0309527c0ad5b39e3e6b1216f76456324a642962ab \
            /bin/sh \
            /tmp/tmpqd7x7e8m/job_working_directory/000/5/tool_script.sh \
            > ../outputs/tool_stdout \
            2> ../outputs/tool_stderr; \
            return_code=$?; 

```

This is the (modified) tool script

```bash
#!/bin/sh

# The following block can be used by the job system
# to ensure this script is runnable before actually attempting
# to run it.
if [ -n "$ABC_TEST_JOB_SCRIPT_INTEGRITY_XYZ" ]; then
    exit 42
fi
set -e

# Check if container was created by installing conda packages,
# and if so, source scripts to populate environment variables
# that would be set by activating the conda environment.
if [ -d /usr/local/etc/conda/activate.d ]; then
  export CONDA_PREFIX=/usr/local
  for f in /usr/local/etc/conda/activate.d/*.sh; do
    case "$f" in
      "/usr/local/etc/conda/activate.d/activate-"*) :;;
      *) . "$f" ;;
    esac;
  done
fi
cp '/tmp/tmpqd7x7e8m/files/1/a/6/dataset_1a6e460a-0062-41eb-81e5-f609722d693a.dat' 'reference.fasta' &&  cp '/tmp/tmpqd7x7e8m/files/a/2/6/dataset_a26ff453-da74-4c0c-9874-19900c336f32.dat' 'reference.fasta.fai' &&  cp '/tmp/tmpqd7x7e8m/files/a/8/e/dataset_a8e7ef63-b3e2-4ef2-a39d-d894c6d51f5c.dat' 'input.bam' &&  cp '/tmp/tmpqd7x7e8m/files/9/9/1/dataset_991127ce-95b3-4b49-82d2-5f295f48404a.dat' 'input.bam.pbi' &&   variantCaller --numWorkers ${GALAXY_SLOTS:-4} --referenceFilename 'reference.fasta' --outputFilename output.fa 'input.bam'

```

Running this from `/tmp/tmpqd7x7e8m/job_working_directory/000/5/working`

```bash
singularity exec    -B /tmp --nv -H $(mktemp -d) --pwd ${PWD} --containall --cleanenv --writable-tmpfs           docker://quay.io/biocontainers/genomicconsensus@sha256:de72299d4fb4f2bd25abdb0309527c0ad5b39e3e6b1216f76456324a642962ab      variantCaller                --numWorkers 32                 --referenceFilename reference.fasta          --outputFilename output.fa              input.bam
```

Gets me the corrupt BAM error

```python
[ERROR] A BGZF (e.g. a BAM file) block should start with '\x1f\x8b\x08\x04', not 'PBI\x01'; handle.tell() now says 4
Traceback (most recent call last):
  File "/usr/local/lib/python2.7/site-packages/pbcommand/cli/core.py", line 138, in _pacbio_main_runner
    return_code = exe_main_func(*args, **kwargs)
  File "/usr/local/lib/python2.7/site-packages/GenomicConsensus/main.py", line 340, in args_runner
    return tr.main()
  File "/usr/local/lib/python2.7/site-packages/GenomicConsensus/main.py", line 258, in main
    with AlignmentSet(options.inputFilename) as peekFile:
  File "/usr/local/lib/python2.7/site-packages/pbcore/io/dataset/DataSetIO.py", line 2746, in __init__
    super(AlignmentSet, self).__init__(*files, **kwargs)
  File "/usr/local/lib/python2.7/site-packages/pbcore/io/dataset/DataSetIO.py", line 1994, in __init__
    super(ReadSet, self).__init__(*files, **kwargs)
  File "/usr/local/lib/python2.7/site-packages/pbcore/io/dataset/DataSetIO.py", line 461, in __init__
    self.updateCounts()
  File "/usr/local/lib/python2.7/site-packages/pbcore/io/dataset/DataSetIO.py", line 2574, in updateCounts
    self.assertIndexed()
  File "/usr/local/lib/python2.7/site-packages/pbcore/io/dataset/DataSetIO.py", line 2404, in assertIndexed
    self._assertIndexed((IndexedBamReader, CmpH5Reader))
  File "/usr/local/lib/python2.7/site-packages/pbcore/io/dataset/DataSetIO.py", line 1951, in _assertIndexed
    self._openFiles()
  File "/usr/local/lib/python2.7/site-packages/pbcore/io/dataset/DataSetIO.py", line 2075, in _openFiles
    resource = IndexedBamReader(location)
  File "/usr/local/lib/python2.7/site-packages/pbcore/io/align/BamIO.py", line 379, in __init__
    self.pbi = PacBioBamIndex(pbiFname)
  File "/usr/local/lib/python2.7/site-packages/pbcore/io/align/PacBioBamIndex.py", line 185, in __init__
    with BgzfReader(pbiFilename) as f:
  File "/usr/local/lib/python2.7/site-packages/pbcore/io/align/_bgzf.py", line 557, in __init__
    self._load_block(handle.tell())
  File "/usr/local/lib/python2.7/site-packages/pbcore/io/align/_bgzf.py", line 584, in _load_block
    block_size, self._buffer = _load_bgzf_block(handle, self._text)
  File "/usr/local/lib/python2.7/site-packages/pbcore/io/align/_bgzf.py", line 422, in _load_bgzf_block
    % (_bgzf_magic, magic, handle.tell()))
ValueError: A BGZF (e.g. a BAM file) block should start with '\x1f\x8b\x08\x04', not 'PBI\x01'; handle.tell() now says 4
```

```bash
singularity exec    \
    -B ${PWD},/tmp --nv -H $(mktemp -d) --pwd ${PWD} --containall --cleanenv --writable-tmpfs           \
    docker://quay.io/biocontainers/genomicconsensus@sha256:de72299d4fb4f2bd25abdb0309527c0ad5b39e3e6b1216f76456324a642962ab \
   variantCaller                \
   --numWorkers 32                 \
   --referenceFilename test-data/All4mer.V2.01_Insert.fa          \
   --outputFilename output.fa              \
   test-data/out.aligned_subreads.bam


```

## Switch to pbgcpp

```bash
# align the reads
singularity exec \
    docker://quay.io/biocontainers/pbmm2:1.8.0--hdfd78af_0 \
    pbmm2 align -j 1 \
    bnd.bam \
    bnd-ref.fasta \
    CRAMTMP/pbmm2.bam \
    --sort -j 1 -J 1 --short-sa-cigar

# try out with pbgcpp
singularity exec    \
    -B ${PWD},/tmp --nv -H $(mktemp -d) --pwd ${PWD} --containall --cleanenv --writable-tmpfs           \
    pbgcpp_2.0.2.sif \
        gcpp \
            --num-threads 32 \
            --reference singularity-test-data/bnd-ref.fasta \
            --output output.fa \
            singularity-test-data/pbmm2.bam





