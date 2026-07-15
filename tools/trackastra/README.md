# Trackastra track

Galaxy wrapper for the upstream `trackastra track` CLI. It accepts paired image/mask time series as ordered TIFF collections or multi-page TIFF datasets and returns an edge table, CTC track table, and tracked-mask collection.

Models are selected from the `trackastra_models` Tool Data Table or supplied from history as `model.pt`, `config.yaml`, and `train_config.yaml`.

## Upgrade note

Version `0.5.3+galaxy2` additionally fixes the Galaxy 26.0 functional tests:

- the bundled history-model checkpoint uses PyTorch's non-ZIP serialization so Galaxy does not auto-extract it during test-data upload;
- tracked-mask tests validate collection membership and non-empty TIFF files without invoking Galaxy's ZSTD image decoder.

It retains the command and collection fixes introduced in `+galaxy1`.
