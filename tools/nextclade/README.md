### Updating the generated_macros.xml

Nextclade outputs differs depending on which dataset is used for clade assignment. To deal with this, the `datasets_to_macros.py` script runs nextclade to get a list of datasets and then runs nextclade with sample data and collects the expected columns from each output dataset. This information is used to generation macros, thus:

```
# requires nextclade to be in the path
./datasets_to_macros.py generated_macros.xml
```

Should be run before updating the nextclade tool. Note that there are a few special cases in here:

1. There are two sets of Influenza datasets with the same name, but different reference sequences. These are special cased in the code to distinguish them.

2. Some information about the SARS-CoV-2 and MPXV is used for testing, thus tokens are created specifically containing info about these datasets.

The `generated_macros.xml` now also includes the tool version, to ensure that the generated data matches the version of nextclade used to generate it.