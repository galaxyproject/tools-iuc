# LIANA+ for Galaxy

Galaxy wrappers for [LIANA+](https://liana.readthedocs.io/en/stable/) **2.0.0**, covering cell-cell communication and spatial analysis. Each tool form provides method-specific parameters and help.

## Tools

| Tool | Use |
| --- | --- |
| [Methods](liana_methods.xml) | Ligand-receptor inference, rank aggregation and spatial statistics |
| [MISTy](misty.xml) | Fit generic or ligand-receptor spatial multi-view models |
| [Multi-view](multi.xml) | Build views or tensors, convert results and run NMF |
| [Plotting](plot.xml) | Visualize interactions, spatial curves and MISTy results |
| [Utilities](utils.xml) | Prepare spatial inputs, extract results and summarize curves |
| [Resources](resource.xml) | Retrieve, filter or translate interaction resources |

## Running an analysis

1. Upload an AnnData (`.h5ad`) file with expression, gene identifiers and observation annotations. For single-cell ligand-receptor inference, use nonnegative, log-normalized expression. Scaled or centered expression is unsuitable.
2. Select the expression source. **Use adata.raw** requires an initialized `.raw`; disable it to use `.X` or a named layer. There is no automatic fallback from `.raw` to `.X`.
3. Select the observation column defining the cell groups, and an interaction resource matching the organism and gene identifiers. Custom ligand-receptor resources use a tab-separated table with `ligand` and `receptor` columns; translation options are available in Resources.
4. Run Methods and pass its annotated output to a compatible Plotting method, keeping the result key consistent. Single-cell methods, including rank aggregation, compare groups without restricting interactions by spatial proximity.
5. Spatial methods need coordinates or a neighbor graph as specified in their help: defaults are `adata.obsm['spatial']` and `adata.obsp['spatial_connectivities']`. Utilities → `spatial_neighbors` builds a graph from coordinates. Cross-PCF and LRIC also output a curves table for `lric_lineplot`.
6. For MISTy, fit a model before plotting its `.h5mu` result. Target-metric plots require a statistic; interaction plots require one view name present in the result.

Built-in and history resources are available without a cache. To populate cached resource selectors, a Galaxy administrator installs and runs the [LIANA Resources Data Manager](../../data_managers/data_manager_liana_resources/README.md). Users then select an installed resource in the tool form.

## Upgrading existing workflows to 2.0.0

Upgrade each workflow step explicitly and review its settings:

- Recheck the expression source and grouping column.
- Cross-PCF and LRIC use analytical nulls; angular sampling and minimum-expressing-cell counts were removed. For LRIC, choose expression-proportion filters explicitly; counts cannot be converted automatically. Annulus width must be a positive integer multiple of the radius step.
- MISTy free-text sorting keys were removed. Recheck the required MuData, statistic and view selections.
- Metalinks annotation filters are available under **Fetch with LIANA**. Reselect that source and re-enter any previous biospecimen filter. Cached resources support the shared edge filters; their presets are already filtered.

## Development

From this directory, run the standard wrapper checks:

```console
planemo lint .
planemo test .
```

See the [LIANA documentation](https://liana.readthedocs.io/en/stable/) and [API reference](https://liana.readthedocs.io/en/stable/api.html) for the underlying methods.
