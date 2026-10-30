## Updating the wrapper

1. Update the VERSION token in macros.xml. The value must be a tag from
   https://github.com/Gaius-Augustus/Tiberius/tags.
2. Regenerate the `options_macros.xml` file by running `python3
   update_options_macros.py macros.xml options_macros.xml`.
3. Run tests etc.
