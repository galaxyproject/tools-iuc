# FoldMason Galaxy tools

The Galaxy wrappers are MIT licensed, as declared by the `license` attributes
and the tools-iuc repository license. FoldMason itself is separately licensed
under [GPL-3.0](https://github.com/steineggerlab/foldmason/blob/main/LICENSE.md).
The wrapper license does not change the license of the installed software.

## ProstT5 inference and GPU configuration

Structure input requires neither a model nor a GPU. Sequence input requires a
FoldMason build with ProstT5 support and compatible GGUF weights registered in
the shared `huggingface` data table; see `tool-data/huggingface.loc.sample` for
installation instructions. The original build
`quay.io/biocontainers/foldmason:4.dd3c235--h5021889_0` lacks ProstT5 support.
The sequence tests check model selection and the missing-model/unsupported-build
error paths; they do not test successful inference.

ProstT5 inference uses the CPU by default (`--gpu 0`). An administrator can set
`GALAXY_FOLDMASON_GPU=1` on a GPU job destination to pass `--gpu 1` instead.
For example, add this to the destination's Galaxy job configuration:

```yaml
env:
  - name: GALAXY_FOLDMASON_GPU
    value: "1"
# For Docker destinations, also pass the value into the container:
params:
  docker_env_GALAXY_FOLDMASON_GPU: "1"
```

That destination also needs a CUDA-enabled FoldMason build and GPU access for
the job/container. Configure GPU allocation and `CUDA_VISIBLE_DEVICES` through
the job runner or scheduler. This wrapper variable enables inference on GPUs;
it does not allocate GPU resources.

## Structure identifiers

Both tools strip `.pdb`/`.cif` (and an optional `.gz` suffix) from input names
and replace characters outside letters, digits, underscores, dots and hyphens
with underscores. These normalized names must be unique. A separate marker
per normalized identifier detects collisions even when the staged files have
different extensions, such as `same.pdb` and `same.cif`. The job fails with an
explanation before running FoldMason if names collide.
