# Galaxy Tool Testing Reference

Standalone reference for planemo testing. Used by tool-dev and other skills (nf-to-galaxy, galaxy-integration).

For the full tool development guide, see `../SKILL.md`. For the complete assertion specification, see the [Galaxy Tool XSD Schema](https://docs.galaxyproject.org/en/latest/dev/schema.html).

**Standard:** the IUC rules for tests are in [tool_xml: Tests](https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/tool_xml.html#tests) and the "Testing Your Tool" part of the [ToolShed Readiness Checklist](https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/integration_checklist.html). This file covers how to write and debug the tests, not what is required.

---

## Running Tests

### Full Test Suite

```bash
# With biocontainers (standard for CLI-wrapping tools)
planemo test --biocontainers tools/{tool_name}/

# Without containers (API tools with fixtures)
planemo test tools/{tool_name}/
```

### Specific Test

```bash
# Run a specific test by 0-based index
planemo test --test_index 0 tools/{tool_name}/{tool}.xml

# Single tool file
planemo test --biocontainers tools/{tool_name}/{tool}.xml
```

### With a Local Galaxy Instance

```bash
planemo test --galaxy_root /path/to/galaxy tools/{tool_name}/
```

---

## Test Output Files

After running tests, planemo produces:
- `tool_test_output.json` — machine-readable detailed results
- `tool_test_output.html` — human-readable report

---

## Writing Tests

### Basic Test

```xml
<tests>
    <test expect_num_outputs="1">
        <param name="input" value="test_input.txt"/>
        <output name="output">
            <assert_contents>
                <has_text text="expected output"/>
            </assert_contents>
        </output>
    </test>
</tests>
```

### With Sections and Conditionals

```xml
<test expect_num_outputs="2">
    <conditional name="output_mode">
        <param name="mode" value="advanced"/>
    </conditional>
    <section name="options">
        <param name="threads" value="4"/>
    </section>
    <output name="output" file="expected_output.txt"/>
</test>
```

**Avoid pipe syntax** for nested inputs. Use explicit nesting:

```xml
<!-- BAD: pipe syntax -->
<param name="filtering|min_score" value="0.5"/>

<!-- GOOD: explicit nesting -->
<section name="filtering">
    <param name="min_score" value="0.5"/>
</section>
```

### With Repeat Elements

```xml
<test expect_num_outputs="1">
    <section name="filters">
        <repeat name="filter_list">
            <param name="filter_value" value="value1"/>
        </repeat>
        <repeat name="filter_list">
            <param name="filter_value" value="value2"/>
        </repeat>
    </section>
    <output name="output">
        <assert_contents>
            <has_text text="value1"/>
        </assert_contents>
    </output>
</test>
```

### Collection Output

```xml
<test expect_num_outputs="1">
    <param name="input" value="test_input.fa"/>
    <output_collection name="results" type="list" count="5">
        <element name="file1">
            <assert_contents>
                <has_text text="expected"/>
            </assert_contents>
        </element>
    </output_collection>
</test>
```

Use `min` instead of exact `count` when output count may vary:

```xml
<output_collection type="list" min="10">
```

---

## Assertion Types

### Content Assertions (`<assert_contents>`)

```xml
<assert_contents>
    <!-- Text matching -->
    <has_text text="must contain this"/>
    <has_text text="must not contain" negate="true"/>
    <has_line line="exact line match"/>
    <!-- Anchored at both ends: the expression must match a whole line (see Gotchas) -->
    <has_line_matching expression="regex.*pattern"/>

    <!-- Structure -->
    <has_n_lines n="10"/>
    <has_n_lines min="5" max="20"/>
    <has_n_columns n="4"/>
    <has_n_columns n="4" sep="\t"/>

    <!-- Size -->
    <has_size value="1000" delta="100"/>
    <has_size min="500"/>
</assert_contents>
```

### File Comparison Modes

```xml
<!-- Exact match (default) -->
<output name="output" file="expected.txt"/>

<!-- Partial match — output must contain expected file's content -->
<output name="output" file="expected_subset.txt" compare="contains"/>

<!-- Allow N lines to differ (timestamps, floats) -->
<output name="output" file="expected.txt" lines_diff="2"/>

<!-- Size-based comparison for binary files -->
<output name="output" file="expected.bam" ftype="bam" compare="sim_size" delta="100"/>
```

### Stream and Command Assertions

```xml
<!-- Assert on stderr content -->
<assert_stderr>
    <has_text text="WARNING: low coverage"/>
</assert_stderr>

<!-- Assert on stdout content -->
<assert_stdout>
    <has_text text="Processing complete"/>
</assert_stdout>

<!-- Assert on the constructed command line -->
<assert_command>
    <has_text text="--expected-flag"/>
</assert_command>
```

### Testing Compressed Outputs

Use `decompress="true"` to test content inside compressed files. It is an attribute of the `<output>` (or collection `<element>`) being checked, not of the individual assertion:

```xml
<output name="output" ftype="fastqsanger.gz" decompress="true">
    <assert_contents>
        <has_text text="expected"/>
    </assert_contents>
</output>
```

It applies to `assert_contents` and to every `compare` mode except `sim_size`, and by default handles gzip, bz2, and zip.

---

## Common Test Failures

### Data Changes Over Time

**Problem**: Upstream database changes, exact counts fail.

**Solution**: Use `min=` instead of exact values.

```xml
<!-- Fragile: exact count -->
<has_n_lines n="142"/>
<output_collection type="list" count="12">

<!-- Robust: minimum -->
<has_n_lines min="140"/>
<output_collection type="list" min="10">
```

### Variable Output

**Problem**: Output varies between runs (timestamps, random seeds, floating point).

**Solution**: Test for required content, not exact match.

```xml
<assert_contents>
    <has_text text="required header"/>
    <has_n_columns n="4"/>
</assert_contents>
```

Or use `lines_diff` for mostly-stable output:

```xml
<output name="output" file="expected.tsv" lines_diff="2"/>
```

### Floating Point

**Problem**: Floating point precision differences across platforms.

**Solution**: Use size-based assertions with delta.

```xml
<has_size value="1000" delta="50"/>
```

### Compressed Output

**Problem**: Testing content of compressed files.

**Solution**: Use `decompress="true"` on the `<output>` or `<element>`. Inside an `<output_collection>`:

```xml
<element name="sample1" decompress="true">
    <assert_contents>
        <has_text text="expected"/>
    </assert_contents>
</element>
```

---

## Tool-Data (Data Table) Tests

**Standard:** [ToolShed Readiness Checklist](https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/integration_checklist.html), "If your tool uses tool-data".

A tool that reads a data table (`<options from_data_table="mytool_indexes">`) cannot be tested until the table exists in the test Galaxy. Two files do that:

`tool_data_table_conf.xml.test`, next to the tool XML. `${__HERE__}` is the tool directory:

```xml
<tables>
    <table name="mytool_indexes" comment_char="#">
        <columns>value, dbkey, name, path</columns>
        <file path="${__HERE__}/test-data/mytool_indexes.loc"/>
    </table>
</tables>
```

`test-data/mytool_indexes.loc`, tab-separated, one row per test index, with the index itself committed under `test-data/`:

```
test_index	hg38_test	Test index	${__HERE__}/mytool_index
```

The test selects the row by its `value` column:

```xml
<test expect_num_outputs="1">
    <param name="reference" value="test_index"/>
    <param name="reads" value="reads.fastq.gz"/>
    <output name="output" file="expected.bam" ftype="bam" compare="sim_size"/>
</test>
```

`planemo test` picks up `tool_data_table_conf.xml.test` from the tool directory by itself. The same tool also ships `tool_data_table_conf.xml.sample` and `tool-data/mytool_indexes.loc.sample` for real installs (see `../SKILL.md` Section 5). Copy a sibling tool's four files rather than writing them from memory; `tools/bowtie2` in tools-iuc is a working example.

---

## Hostile-Input Test

**Standard:** [security item 2](https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/security.html#security-free-text): test punctuation, whitespace, quotes, newlines, empty values, and boundary lengths, not only the default.

Every free-text parameter gets two kinds of test.

**Rejected.** A value outside the validator's language never reaches the shell. The test expects failure and has no output checks:

```xml
<test expect_failure="true">
    <param name="input" value="input.tsv"/>
    <param name="sample_name" value="ok; touch pwned"/>
</test>
```

Useful hostile values: `a;b`, `$(id)`, `` `id` ``, `a'b`, `a"b`, `a b`, `../x`, `-x` (looks like a flag), and a trailing-junk value such as `good;rm` that only an end-anchored validator rejects. If this test *passes without* `expect_failure`, the validator is not anchored.

**Accepted and intact.** The awkward values the validator does allow (spaces, dots, the longest reasonable value, the empty string for an optional param) arrive as one argument, unchanged:

```xml
<test expect_num_outputs="1">
    <param name="input" value="input.tsv"/>
    <param name="sample_name" value="Sample 1.a"/>
    <assert_command>
        <has_text text="--name 'Sample 1.a'"/>
    </assert_command>
    <output name="output">
        <assert_contents>
            <has_text text="Sample 1.a"/>
        </assert_contents>
    </output>
</test>
```

> **Verify:** a validator rejection happens before a job exists. Current Galaxy counts that as the expected failure for `expect_failure="true"` (checked in the test interactor of a 26.1 development checkout); this has not been checked on older Galaxy or planemo versions. If the rejected-value test errors instead of passing, say so in the hand-back rather than deleting the test.

---

## Gotchas

Field findings from real wrapper work. None of these is in the IUC standards, and most are not in the Galaxy or planemo docs either. Each entry says when to delete it.

### A test `<collection>` only binds to a `data_collection` param

A `type="data" multiple="true"` param accepts a list collection in the UI, but a test-side `<collection>` block cannot feed it. Give it a comma-separated list of files instead:

```xml
<!-- multiple="true" data param -->
<param name="inputs" value="a.bam,b.bam"/>

<!-- data_collection param -->
<param name="paired_input">
    <collection type="paired">
        <element name="forward" value="r1.fastq.gz"/>
        <element name="reverse" value="r2.fastq.gz"/>
    </collection>
</param>
```

*Delete when the planemo or Galaxy test docs state which param types a test `<collection>` can bind to.*

### Unsniffable inputs need an explicit `ftype`

Galaxy guesses the datatype of each uploaded test file. Files it cannot sniff (Salmon `quant.sf`, custom tabular or text formats) can end up with a type the param does not accept. Set it:

```xml
<param name="quant" value="quant.sf" ftype="tabular"/>
```

*Delete when the planemo docs cover `ftype` on test inputs.*

### `has_line_matching` is anchored at both ends; `has_text` is literal

`has_line_matching` wraps the expression as `^expression$` per line, so `expression="chr1"` does not match the line `chr1\t100`. Write `chr1\t.*`, or use `has_text_matching` for an unanchored search. `has_text` is a literal substring match: regex metacharacters in `text=` are not interpreted.

*Delete when the schema docs for these assertions say so explicitly.*

### `planemo test` fails building the Galaxy virtualenv

Symptom: the run dies before any test, building packages such as `h5py` or `pyo3`, on a machine whose default `python3` is newer than the Galaxy release supports (seen with Python 3.14 and Galaxy 25.1). Pin the interpreter planemo uses for Galaxy:

```bash
planemo test --galaxy_python_version 3.12 --biocontainers tools/mytool/
```

This is troubleshooting for a local environment problem, not part of the standard command. Do not add it to documented or CI commands.

*Delete when planemo picks a supported interpreter by itself or documents this failure.*

### No mulled image on quay for a multi-requirement tool

Symptom: `planemo test --biocontainers` cannot find a container for a tool with two or more requirements. The combination has not been registered and built yet (see `../SKILL.md` Section 9). Build it locally, then re-run the test:

```bash
planemo mull tools/mytool/mytool.xml
planemo test --biocontainers tools/mytool/mytool.xml
```

Say in the hand-back that the image was built locally, so the reviewer knows the public image is still pending.

> **Verify:** `planemo mull` is taken from planemo's own help text ("build containers for specified tools ... Galaxy can automatically discover this container"); it was not run end to end for this skill. The field report used `mulled-build build` directly; its exact arguments were not recorded, so they are not shown here.

*Delete when the planemo docs describe testing a tool whose mulled container is not published yet.*

### Python `csv.writer` writes CRLF

`csv.writer` defaults to `\r\n` line endings. In a wrapper script that produces tabular output, this breaks `$`-anchored assertions and line-based assertions, and confuses downstream tabular tools. Pass `lineterminator="\n"`:

```python
writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
```

*Delete when this is covered in the Galaxy tool-development docs or the IUC standards.*

### Generator expressions do not compile in Cheetah

`#if any(str($b.code) for $b in $barcodes)` fails to compile (seen on Galaxy 25.1). The failure shows up as a tool that will not load or a command that cannot be built, not as a test assertion. Use an explicit loop:

```
#set $has_code = False
#for $b in $barcodes
    #if str($b.code).strip()
        #set $has_code = True
    #end if
#end for
```

*Delete when the Galaxy Cheetah documentation lists unsupported Python constructs.*

### `find ... -exec ... \;` is mangled in a command block

Reported on Galaxy 25.1; the exact mechanism was not investigated. Replace it with a shell glob or loop, which also keeps file enumeration at runtime where remote execution needs it:

```
mkdir -p out && mv work/*.fastq out/
for f in work/*.fastq; do gzip "\$f"; done
```

*Delete when the Galaxy command-block docs cover escaping of `\;`.*

### `version_command` needs CDATA when it redirects

`2>&1` contains `&`, which is not valid in XML text. `<command>` is always in CDATA, so the problem only shows up in `version_command`:

```xml
<version_command><![CDATA[mytool --version 2>&1 | head -n 1]]></version_command>
```

*Delete when the IUC standards or schema docs show CDATA on `version_command`.*

---

## Analyzing Test Failures

### Quick Summary from tool_test_output.json

```python
import json

with open("tool_test_output.json") as f:
    data = json.load(f)

print(f"Summary: {data.get('summary', {})}")

for t in data.get("tests", []):
    status = t.get("data", {}).get("status")
    if status != "success":
        print(f"\n=== {t.get('id')} ===")
        print(f"Status: {status}")
        print(f"Problems: {t.get('data', {}).get('output_problems', [])}")
```

### Detailed Job Info

```python
for t in data.get("tests", []):
    if t.get("data", {}).get("status") != "success":
        job = t.get("data", {}).get("job", {})
        if job.get("stderr"):
            print(f"stderr: {job['stderr'][-500:]}")
        if job.get("stdout"):
            print(f"stdout: {job['stdout'][-500:]}")
```

---

## Test Numbering

Tests are 0-indexed:
- Test 0 = First `<test>` in file
- Test 1 = Second `<test>` in file
- etc.

When `tool_test_output.json` reports `test-5`, that's the 6th test.

---

## Test Structure Rules

- Give each test a **unique purpose** (e.g., "defaults", "compression on", "filtering active")
- Always include `expect_num_outputs` to verify file counts
- **Test data: aim for under 1 MB in total.** Shrink inputs first. For large outputs, use assertions or `compare="contains"` instead of a golden file. For a large input that cannot be made smaller, deposit it in the [Galaxy Tool Developers Zenodo community](https://zenodo.org/communities/galaxy_tool_dev_sig) and reference it with `<param name="input" location="https://zenodo.org/records/RECORD/files/FILE_NAME"/>`
- Include tests for output filters to verify filtering actually occurs
- Include tests for error conditions with expected failure
- Include a hostile-input test for every free-text parameter
- Include a test with compressed input and one with uncompressed input when the tool accepts both

---

## Related

- **Full tool development guide**: `../SKILL.md`
- **Tool placement guide**: `tool-placement.md`
