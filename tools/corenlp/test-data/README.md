# Test Data for Stanford CoreNLP

## Running Tests Locally

The tests require the Stanford CoreNLP English language models. Download them before running tests:

```bash
cd test-data
curl -L -o stanford-corenlp-4.5.10-models-english.jar \
  https://repo1.maven.org/maven2/edu/stanford/nlp/stanford-corenlp/4.5.10/stanford-corenlp-4.5.10-models-english.jar
```

Then run tests from the tool directory:

```bash
cd ..
planemo test --docker
```

## ToolShed Submission

The model JAR file (424MB) is too large to include in the repository. For ToolShed automated testing:

1. The `tool_data_table_conf.xml.test` file points to `test-data/corenlp_models.loc`
2. The `.loc` file references `stanford-corenlp-4.5.10-models-english.jar` using `${__HERE__}`
3. Download the model JAR to `test-data/` before running `planemo shed_test`

### Automated Testing on ToolShed

When submitting to ToolShed, automated tests may not work due to the large model file. Consider:

1. Including a download script in the repository documentation
2. Using manual verification for tool releases
3. Setting up a test data cache if available on the test infrastructure

## Files

- `corenlp_models.loc` - Test data table for language models
- `input.txt` - Sample text for basic tests
- `2.txt` - Multi-sentence sample
- `sa-input.txt` - Sample for sentiment analysis
- `*.conll`, `*.conllu`, `*.json` - Expected output files
