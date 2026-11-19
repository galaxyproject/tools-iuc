---
name: tool-dev
description: How to author, test, and update Galaxy tool wrappers so they meet the IUC standards, with pre-PR passes for security, ToolShed readiness, remote execution, and review flags. Use when authoring or modifying Galaxy tool wrappers for tools-iuc; defers to the IUC standards for the rules.
---

# Galaxy Tool Development Reference

The [IUC standards][std] say what a correct wrapper is. This skill says how an agent gets there and proves it: worked examples, gotchas the standards do not cover, decisions to make, and the checks to run and record before hand-back. It does not restate the standards. Where a section opens with a **Standard:** line, read that page for the rule; the text below it is examples and practice only.

## When to Use This Skill

- Creating a new Galaxy tool wrapper from scratch
- Wrapping a CLI bioinformatics tool or an external API
- Modifying an existing tool for IUC submission
- Updating an existing tool to a new upstream version
- Debugging planemo lint or test failures
- Preparing a tools-iuc PR
- Reviewing Galaxy XML conventions

Data managers are out of scope here. Follow [Data Manager Best Practices][std-dm] (naming `data_manager_$name`, stable identifiers with no wall-clock values, tests with stable expected values).

## Authority

Order of precedence:

1. **The target repository's convention beats everything below.** If sibling tools in the repo you are contributing to do something differently (help format, `detect_errors`, inline tokens vs `macros.xml`, profile), follow the repo and **say so in the hand-back**, naming the convention and the sibling tool you took it from.
2. **The IUC standards beat this skill.** If the two disagree, the standards win; report the disagreement so the skill can be fixed.
3. **This skill** covers how to apply the standards, plus agent practice they omit.

| Standards page | Covers |
|---|---|
| [tool_xml][std-tool-xml] | Versions, IDs, profile, bio.tools, required files, creator, params, tests, command, help, coding style |
| [security][std-security] | Tool Security Checklist (ten items) -- drives pass 10a |
| [integration_checklist][std-checklist] | ToolShed Readiness Checklist -- drives pass 10b |
| [pulsar][std-pulsar] | Remote Execution Checklist (six items) -- drives pass 10c |
| [package_xml][std-package] | Requirements, Bioconda/conda-forge, multi-requirement (mulled) containers |
| [licensing][std-licensing] | OSI rule, licensing layers, use restrictions |
| [shed_yml][std-shed] | `.shed.yml` fields, naming, auto-generated suites, include/exclude |
| [repositories][std-repos] | Repository layout, README, suites |
| [data_managers][std-dm] | Data managers (out of scope for this skill) |
| [management][std-management] | Publishing to the ToolShed (IUC CI does this; do not upload by hand) |

Last aligned against galaxy-iuc/standards 3f6039ac (2026-10-21) plus [PR #95][std-pr95]. PR #95 is a proposal (Markdown help preferred for new tools, when `exit_code` is right, the `@PROFILE@` token, the profile ceiling); where this skill cites it, the published standards may not say it yet.

Other references:

- [Galaxy Tool XSD Schema](https://github.com/galaxyproject/galaxy/blob/dev/lib/galaxy/tool_util/xsd/galaxy.xsd) — validate XML against this
- [Galaxy Core Test Tools](https://github.com/galaxyproject/galaxy/tree/dev/test/functional/tools) — reference implementations and edge cases
- [Planemo Documentation](https://planemo.readthedocs.io/) — lint, test, serve, shed workflows
- `references/testing.md` — standalone planemo testing reference (also used by other skills)
- `references/tool-placement.md` — where to create tools decision guide
- `references/tool-help-diagrams.md` -- SVG diagrams for tool help (RST help only)
- `tool-selection-diagram/SKILL.md` — generate tool selection flowchart diagrams for multi-tool suites

---

## 1. Directory Structure

**Standard:** [repositories][std-repos] for layout, [shed_yml][std-shed] for suites.

All tools in a suite live in **one flat directory** under `tools/<tool_name>/`. No subdirectories for individual tools.

```
tools/mytool/
├── macros.xml                    # Only when two or more tools share something (see Section 3)
├── mytool_align.xml              # Tool wrapper A
├── mytool_filter.xml             # Tool wrapper B
├── .shed.yml                     # Tool Shed metadata
├── static/
│   └── images/
│       └── mytool_tool_selection.png   # Optional: help image (forces RST help, see Section 2)
├── tool-data/                    # Only when the tool uses a data table (see Section 5)
│   └── mytool_indexes.loc.sample
├── tool_data_table_conf.xml.sample     # Only with a data table
├── tool_data_table_conf.xml.test       # Only with a data table
└── test-data/
    ├── input.fastq.gz            # Shared test inputs
    ├── reference.fasta
    ├── expected_align.bam         # Golden file outputs
    └── expected_filter.bed
```

Add `.py` wrapper scripts only when the upstream CLI can't produce the output Galaxy needs (format conversion, multi-step pipelines, etc.) — most tools don't need them.

**Key rules:**
- Tool IDs contain only `[a-z0-9_-]`; separate words with underscores or dashes and match what sibling tools in the repo do ([Tool ids][std-ids])
- Test data goes in a single shared `test-data/` directory
- A shared `macros.xml` is per suite, not per tool. A single-tool directory can keep its tokens inline instead (Section 3)
- One `.shed.yml` per suite

> **Verify:** [repositories][std-repos] says every unrestricted ToolShed repository should contain a `README.rst`. Many tools-iuc tool directories have none, and it is not clear how strictly reviewers apply it. Check sibling tools before adding one.

---

## 2. XML Wrapper Structure

**Standard:** [tool_xml][std-tool-xml]. Element and attribute order: [Coding Style][std-style].

### Element Order (Strict)

Galaxy XML elements must appear in the standards' order. `planemo lint` checks it. The elements an agent normally writes, in order:

```xml
<tool id="..." name="..." version="..." profile="...">
    <description>...</description>
    <macros>
        <import>macros.xml</import>   <!-- or inline <token>s for a single tool -->
    </macros>
    <xrefs>
        <xref type="bio.tools">mytool</xref>
    </xrefs>
    <expand macro="requirements"/>
    <required_files>                  <!-- only when the wrapper ships scripts or config files -->
        <include path="mytool_convert.py"/>
    </required_files>

    <version_command>mytool --version</version_command>

    <command detect_errors="aggressive"><![CDATA[
        ...
    ]]></command>

    <inputs>
        ...
    </inputs>

    <outputs>
        ...
    </outputs>

    <tests>
        ...
    </tests>

    <help format="markdown"><![CDATA[
        ...
    ]]></help>

    <expand macro="citations"/>
    <creator>                         <!-- optional, see below -->
        ...
    </creator>
</tool>
```

`<configfiles>` and `<environment_variables>`, when used, go between `<command>` and `<inputs>`. Section 13 has the rest.

### Tool Element Attributes

```xml
<tool id="mytool_function"
      name="MyTool Function"
      version="@TOOL_VERSION@+galaxy@VERSION_SUFFIX@"
      profile="@PROFILE@">
```

- `id`: [Tool ids][std-ids]. Prefix with the suite name for multi-tool suites (e.g., `bedtools_intersect`).
- `name`: Human-readable, title case. How users and admins find the tool — be specific, not generic.
- `version`: [Tool versions][std-versions]. Use `@TOOL_VERSION@+galaxy@VERSION_SUFFIX@` built from tokens (inline or in `macros.xml`, see Section 3), and reuse `@TOOL_VERSION@` in the `<requirement>`. Exception: a tool that is not a wrapper for one underlying tool omits the `+galaxy` suffix and uses any PEP 440 version. When to bump is in Section 11.
- `profile`: [Tool profile][std-profile]. Use a recent profile, but not one newer than the major public servers run: Galaxy refuses to load a tool whose profile is newer than the server, so a newer profile is not a goal in itself ([PR #95][std-pr95]). Do not copy a number from this skill. Look at what sibling tools in the target repo use, and go higher only when the wrapper needs a feature from a later profile (for example `format="markdown"` help needs 24.2). In a suite, put it in a `@PROFILE@` token so all tools move together; a single tool can set it inline.

### bio.tools Cross-References

**Standard:** [Tool cross-references][std-biotools].

Link to the upstream tool's bio.tools entry with `<xrefs>`. Reference only the tool being wrapped, not its dependencies. Create the bio.tools entry if none exists (needs at least one EDAM Topic and Operation, specific terms rather than root terms). Improve EDAM annotation in the bio.tools entry, not with `edam_topics`/`edam_operations` elements in the XML -- the standards call that not recommended.

```xml
<xrefs>
    <xref type="bio.tools">mytool</xref>
</xrefs>
```

### `<required_files>`

**Standard:** [Required files][std-required-files]; how to get it right: [pulsar item 1][pulsar-1].

Every file the wrapper ships and reaches through `$__tool_directory__` needs an `<include>`: the entry script, helpers it imports or sources, config files, files copied in with `cp`. Not macros. This applies to any wrapper-shipped file, not only Python. The check is one include per `$__tool_directory__` reference across `<command>`, `<version_command>`, and `<configfiles>`, with the same relative path, and no dead includes.

Files under `test-data/` and `tool-data/` are excluded by default, so a script or fixture the job needs at runtime cannot live there.

### `<creator>`

**Standard:** [Creator][std-creator].

Optional. Reviewers invite it ("you can add yourself as creator") but do not require it. It goes last, after `<citations>`. Do not invent identities: add a person only when the human you are working for asks, with their real ORCID or URL, and an email only with their explicit consent.

```xml
<creator>
    <person givenName="Ada" familyName="Lovelace" identifier="https://orcid.org/0000-0000-0000-0000"/>
    <organization name="Example Lab" url="https://example.org"/>
</creator>
```

### Command Block

**Standard:** [Command tag][std-command], [Exit Code Detection][std-exit].

Tripwire: every text parameter, input path, and output path goes in **single** quotes. Double quotes are not equivalent. See Section 4 and [security item 1][sec-1].

**`detect_errors`:** `aggressive` is the default. It fails the job on a non-zero exit code or when `error:` or `exception:` appears in the tool's output (case-insensitive; the standards say stderr, and Galaxy's `stdio.py` checks stdout as well). A line that only says `WARNING` is not a trigger. Use `exit_code` instead in two cases, and put a comment above `<command>` saying which one applies ([PR #95][std-pr95]):

- the tool prints `error:` or `exception:` text on stderr during normal, successful operation;
- reporting errors is the tool's job (validators, linters, QC reporters).

"The tool is chatty" is not a reason by itself. Check the actual stderr of a successful run for the two trigger strings before switching, and check what sibling tools for the same upstream do.

```xml
<!-- exit_code: mytool prints "0 error: none found" on stderr when a run succeeds -->
<command detect_errors="exit_code"><![CDATA[
```

Wrap the command in CDATA and chain multiple commands with `&&`. Most tools call the upstream binary directly:

```xml
<command detect_errors="aggressive"><![CDATA[
    mytool align
        --input '$input_fastq'
        --reference '$reference'
        --output '$output_bam'
        --threads \${GALAXY_SLOTS:-1}
        #if str($min_quality):
            --min-quality $min_quality
        #end if
        @CMD_OPTIONAL_FLAG@
]]></command>
```

For multi-step commands, chain with `&&`:

```xml
<command detect_errors="aggressive"><![CDATA[
    ln -s '$input_fasta' input.fa &&
    mytool index input.fa &&
    mytool align input.fa --output '$output_bam'
]]></command>
```

When the CLI can't produce the output Galaxy needs (format conversion, multi-step pipelines), use a Python wrapper script and declare it in `<required_files>`:

```xml
<command detect_errors="aggressive"><![CDATA[
    python '$__tool_directory__/mytool_convert.py'
        --input '$input_file'
        --output '$output_file'
]]></command>
```

**`version_command`:** use the tool's real version command (`mytool --version`). Reviewers always ask for the real command. `echo @TOOL_VERSION@` is a fallback only for tools with no way to print a version. Wrap it in CDATA when it contains `&` or `<`, for example when the version goes to stderr:

```xml
<version_command><![CDATA[mytool --version 2>&1 | head -n 1]]></version_command>
```

**Parameter parity rule:** Every param in `<inputs>` must appear in `<command>`, and every flag in `<command>` must trace back to an `<inputs>` param or a macro token. Orphaned params are a lint warning and a review flag.

**Index generation:** [Generating Indices][std-indices]. Symlink inputs into the working directory and index there; inputs are read-only.

**Composite datatypes:** inputs use `$input.extra_files_path`, outputs use `$output.files_path` ([Command tag][std-command]).

### Output Paths and `from_work_dir`

**Standard:** [pulsar items 5 and 6][pulsar-5]. Outputs are real files inside the job directory. Never `ln -s` an output to an input or to anything outside the job tree, and never write a discovered output to an absolute path.

When a CLI tool writes to a fixed filename or prefix, use a staging directory so Galaxy can find the output predictably:

```xml
<command detect_errors="aggressive"><![CDATA[
    mkdir -p staging &&
    mytool --output-prefix staging/result '$input_file' &&
    mv staging/result.tsv '$output_file'
]]></command>
```

Or use `from_work_dir` on the output:

```xml
<outputs>
    <data name="output_file" format="tabular" from_work_dir="result.tsv"
          label="${tool.name} on ${on_string}"/>
</outputs>
```

HTML, SVG, and JavaScript outputs are an administrator trust decision, not a free choice. See [security item 10][sec-10] before emitting one.

### Conditional Outputs (Filters)

Use `<filter>` to create outputs only when certain params are set:

```xml
<outputs>
    <data name="output_log" format="txt" label="${tool.name} log">
        <filter>output_log == True</filter>
    </data>
</outputs>
```

Set `expect_num_outputs` to the number of outputs actually produced by each test case. Outputs whose `<filter>` evaluates to False are not produced and should not be counted.

### Dynamic Output Discovery

For tools producing variable numbers of output files, use `discover_datasets`. The `directory` must be inside the job working directory:

```xml
<collection name="split_output" type="list" label="Split files">
    <discover_datasets pattern="__name_and_ext__" directory="output_dir"/>
</collection>
```

### Help Section

**Standard:** [Help tag][std-help], which describes reStructuredText only at the pinned commit; Markdown as the preferred format for new tools is proposed in [PR #95][std-pr95]. Reviewers already offer it ("you can use markdown if this is any easier") and nobody asks for conversion either way.

Pick the format by this table, then write in that format's idioms only:

| Situation | Format |
|---|---|
| New tool, no help images | Markdown: `<help format="markdown">`. Galaxy supports it from 24.2, so the tool's profile must be >= 24.2 |
| New tool that ships help images | RST. Galaxy rewrites `.. image::` paths only for RST help (still true in 26.2), so images do not resolve in Markdown help |
| Editing an existing RST tool | Keep RST. Convert only if you convert the whole help section, and only when the tool has no images. Existing RST tools need not convert |
| Target repo or brief asks for RST | RST, and note it in the hand-back |

Markdown idioms: `##` / `###` headings, single backticks for inline code, fenced blocks for examples, pipe tables, and a blank line before every list (reviewers flag a missing one). RST idioms: `**bold**` pseudo-headers or underlined titles, `-----` rules between sections, double backticks for code, `.. image::` for images. Do not mix them: `-----` rules and double backticks in a `format="markdown"` section are the usual mistake.

Either way: wrap in CDATA, say what the tool does, what each input and output is, and keep it concise. Citations go in `<citations>`, not the help text. For images and diagrams see [Including Images][std-images] and `references/tool-help-diagrams.md`.

---

## 3. macros.xml Patterns

**Standard:** [Tool versions][std-versions] (tokens), [integration_checklist][std-checklist] (a macros file is for avoiding duplication).

### Tokens Are Required, the File Is Not

The standards require the version *tokens*, not a `macros.xml` *file*. Reviewers ask for shared things to go in macros, and question a macro file when nothing is reused.

| Tool directory | Where tokens and macros go |
|---|---|
| One tool, nothing shared | Inline `<macros><token name="@TOOL_VERSION@">...</token>...</macros>` in the tool XML |
| Two or more tools, or anything reused | Shared `macros.xml`, imported with `<import>macros.xml</import>` |

Follow the target repo if its single-tool directories do it differently.

### Tokens vs XML Macros

Get this right from the start.

**`<token>`** — For simple text substitution and Cheetah snippets. Expanded inline with `@NAME@` syntax. Use for version strings, Cheetah command fragments.

**`<xml>`** — For structured XML element trees. Expanded with `<expand macro="name"/>`. Use for parameters, requirements, citations.

**The convention:** anything containing Cheetah template logic (`#if`, `#for`) goes in a `<token>`, not an `<xml>`. Reviewers ask for this. (An `<xml>` macro body is parsed as XML, so a bare Cheetah fragment there is fragile at best; this skill has not verified exactly how it fails.)

```xml
<macros>
    <!-- TOKENS: version strings and Cheetah snippets -->
    <token name="@TOOL_VERSION@">1.2.3</token>
    <token name="@VERSION_SUFFIX@">0</token>
    <!-- Example value only. Choose it by the profile rule in Section 2. -->
    <token name="@PROFILE@">24.2</token>

    <!-- GOOD: Cheetah snippet as token -->
    <token name="@CMD_OPTIONAL_FLAG@"><![CDATA[
        #if str($optional_param).strip()
            --optional-flag '$optional_param'
        #end if
    ]]></token>

    <!-- BAD: Cheetah in xml macro -- reviewers will ask for a token -->
    <xml name="cmd_optional_flag">
        #if str($optional_param).strip()
            --optional-flag '$optional_param'
        #end if
    </xml>

    <!-- XML MACROS: structured element trees -->
    <xml name="requirements">
        <requirements>
            <requirement type="package" version="@TOOL_VERSION@">mytool</requirement>
            <yield/>  <!-- tools can inject additional requirements -->
        </requirements>
    </xml>

    <xml name="organism_param">
        <param name="organism" type="select" label="Organism">
            <option value="human" selected="true">Human (hg38)</option>
            <option value="mouse">Mouse (mm10)</option>
        </param>
    </xml>

    <xml name="citations">
        <citations>
            <citation type="doi">10.1234/example</citation>
        </citations>
    </xml>
</macros>
```

### Named Yields and Token Parameterization

For complex macros, use **named yields** to inject content into specific slots:

```xml
<!-- In macros.xml -->
<xml name="complex_inputs">
    <param name="shared_param" type="text" label="Common input"/>
    <yield name="extra_params"/>
    <param name="shared_flag" type="boolean" label="Common flag"/>
</xml>

<!-- In tool XML — content goes into the named yield slot -->
<expand macro="complex_inputs">
    <yield name="extra_params">
        <param name="special_option" type="select" label="Tool-specific"/>
    </yield>
</expand>
```

Use **token parameterization** on xml macros to pass values into the macro at expansion time. Token parameters are replaced with the value passed at expansion — use them for attribute values like labels, defaults, and formats:

```xml
<!-- In macros.xml -->
<xml name="score_param" tokens="default_score,score_help">
    <param argument="--min-score" type="float" value="@DEFAULT_SCORE@" min="0.0" max="1.0"
           label="Minimum score" help="@SCORE_HELP@"/>
</xml>

<!-- In tool XML — each tool can set its own default and help text -->
<expand macro="score_param" default_score="0.5" score_help="Filter results below this threshold"/>
```

### API Tool Macros

The following macros are only needed when wrapping external APIs (not CLI tools).

**Standard:** [security item 8][sec-8]. Secrets never appear in tool XML, test fixtures, command-line arguments, logs, or a templated `<configfile>`. The tool reads them from the environment.

**Credentials** -- Galaxy's vault-backed secrets system (`<credentials>` arrived in Galaxy 25.1, so it needs a server at 25.1 or later). Place `<credentials>` inside `<requirements>` in the macro, not as a separate macro:

```xml
<xml name="requirements">
    <requirements>
        <requirement type="package" version="@TOOL_VERSION@">mytool</requirement>
        <credentials name="mytool" version="1.0" label="MyTool API" description="API key for MyTool service">
            <secret name="api_key" inject_as_env="MYTOOL_API_KEY" label="API Key" description="Your MyTool API key"/>
        </credentials>
    </requirements>
</xml>
```

Tools then just use `<expand macro="requirements"/>` -- no yield or separate credentials expand needed. Do not expand `\$MYTOOL_API_KEY` onto the command line; let the script read the environment.

**Test fixture param** — hidden param for fixture-based CI testing (bypasses API calls):

```xml
<xml name="test_fixture_param">
    <param name="test_fixture" type="hidden" value=""/>
</xml>

<token name="@CMD_TEST_FIXTURE@"><![CDATA[
    #if str($test_fixture).strip()
        --test-fixture '$__tool_directory__/fixtures/$test_fixture'
    #end if
]]></token>
```

Three things this pattern needs:

- Fixtures live in a directory the job can reach (`fixtures/` here), declared in `<required_files>` (`<include path="fixtures/" type="prefix"/>`). Files under `test-data/` are excluded from staging by default, so a fixture read through `$__tool_directory__/test-data/...` is missing under remote execution.
- Recorded fixtures are scrubbed before commit (Section 8).
- The script accepts only a bare filename for `--test-fixture` and rejects path separators.

> **Verify:** whether a `type="hidden"` param can be set by a user through the API. If it can, `$test_fixture` is a user-controlled value concatenated into a tool-directory path, which is why the script-side filename check above is not optional. This skill has not confirmed it either way.

---

## 4. Cheetah Templating

**Standard:** [security items 1, 4, 5, 6][std-security]; [Command tag][std-command]. Cheetah should be indented and mainly PEP 8 conformant ([Coding Style][std-style]).

### Dollar Sign Escaping

This causes more bugs than anything else in Galaxy XML. The rules:

| Context | Syntax | Example |
|---------|--------|---------|
| Text param, dataset path, output path, collection value | `'$param'` | `--input '$input_file'` |
| `integer` / `float`, boolean `truevalue`/`falsevalue`, fixed-option `select` only | `$param` | `--count $max_count` |
| Shell environment variable | `\${VAR}` | `--workers \${GALAXY_SLOTS:-1}` |
| Cheetah loop variable | `'$item'` | `'$item'` inside `#for` |
| Literal dollar in shell | `\$` | `awk '{print \$1}'` |

**Critical:** Galaxy parameters are single-quoted (`'$param'`). Double quotes (`"$param"`) are not a substitute: the shell still expands inside them ([security item 1][sec-1]). An unquoted parameter is acceptable only when its type or finite value set makes it safe, as in row two; a `select` with dynamic options or a `text` param never qualifies. Environment variables use backslash-escaped dollar with braces (`\${VAR:-default}`).

### Nested Contexts and Configfiles

**Standard:** [security item 1][sec-1] (nested interpreters), [item 6][sec-6] (data, not code).

Single quotes protect a value from the outer shell only. A parameter placed inside `bash -c '...'`, an `awk` or `sed` program, `Rscript -e`, `python -c`, or a generated script is parsed again by that inner language, and the outer quoting does nothing there. Do not interpolate parameters into program text. Hand them over as data:

```xml
<command detect_errors="aggressive"><![CDATA[
    Rscript '$__tool_directory__/plot.R' '$params_json' '$output'
]]></command>
<configfiles>
    <inputs name="params_json" data_style="paths"/>
</configfiles>
```

`<configfiles><inputs>` writes every parameter to a JSON file that the script parses, so no user value is ever program text. A templated `<configfile>` is acceptable only when its content is consumed strictly as data (a key/value table, a config format). A `<configfile>` that renders R, Python, or shell with `$param` inside it is generated code and has the same problem as `bash -c`. Never offer `eval`, a script-body param, or a free-text "extra options" param. Secrets never go in a configfile.

### Keep Cheetah Side-Effect Free

**Standard:** [security item 5][sec-5], [pulsar items 2 and 3][pulsar-2].

Cheetah in `<command>` and `<configfile>` runs on the Galaxy server before the job starts, with Galaxy's privileges and the server's filesystem. In a template: no `#import os` / `subprocess` / `pathlib` / `glob`, no `open(`, no `os.path.*` on data paths, no `$__app__`. Use dataset metadata (`$input.get_size()`, `$input.metadata.*`, `$input.is_of_type(...)`, `$input.element_identifier`). When metadata is not enough, do the work at runtime: a shell glob in the command, or plain script text in a configfile. `#import re` and `#import json` for pure string work are fine.

### Conditional Parameters

For text params, always use `str()` and `.strip()` — never test the raw param directly:

```
#if str($optional_text).strip()
    --flag '$optional_text'
#end if
```

Convert with `str()` **before** calling any method (`.split()`, `.strip()`, joins, formatting), then single-quote each resulting shell argument: method calls on the wrapped value can bypass sanitization ([security item 4][sec-4]). The aim is a value sanitized once for the context that consumes it, not `str()` everywhere.

For select params, use string comparison:

```
#if str($mode) == "advanced"
    --advanced-flag
#end if
```

### Boolean Parameters

**Standard:** [Booleans][std-booleans].

Put the CLI flag in `truevalue`/`falsevalue` (Section 5) and emit the parameter directly. No `#if` is needed:

```
mytool $gzip '$input'
```

`str($flag)` is whatever `truevalue`/`falsevalue` say, so it is `"true"`/`"false"` only when those attributes are left at their defaults. Comparing `str($flag) == "true"` on a boolean with `truevalue="--gzip"` is always false. When you need to branch on a boolean, test it directly:

```
## GOOD: works whatever truevalue is
#if $my_flag
    --enable-feature --feature-level 3
#end if

## BAD once truevalue/falsevalue are set: the string is the flag, never "true"
#if str($my_flag) == "true"
```

The standards use both idioms (the licensing page has `#if not $non_commercial_use` and a string comparison on a default-valued boolean), and neither is declared canonical. Follow sibling tools in the target repo.

> **Verify:** `#if $flag` is truthy by the checked state in current Galaxy (the parameter wrapper's `__bool__` returns `bool(value)`). This skill has not checked it against old profiles.

### Numeric Zero Gotcha

Integer and float params: `0` is a valid value but falsy in Cheetah/Python. Use `str()` to avoid skipping zero:

```
## BAD: skips the flag when value is 0
#if $min_score
    --min-score $min_score
#end if

## GOOD: treats 0 as a valid value
#if str($min_score):
    --min-score $min_score
#end if
```

### For Loops (Multi-Select Parameters)

```
--output-types
#for $fmt in $output_types
    '$fmt'
#end for
```

### Conditionals (Tool Sections)

```xml
<conditional name="output_mode">
    <param name="mode" type="select" label="Output mode">
        <option value="summary" selected="true">Summary</option>
        <option value="binned">Binned</option>
    </param>
    <when value="summary"/>
    <when value="binned">
        <param name="bin_size" type="integer" value="128" min="1" max="4096"
               label="Bin size (bp)"/>
    </when>
</conditional>
```

In the command block:

```
--output-mode '$output_mode.mode'
#if str($output_mode.mode) == "binned"
    --bin-size $output_mode.bin_size
#end if
```

### Field Gotchas

Seen on Galaxy 25.1 during wrapper work; not in the standards.

- **Generator expressions do not compile in Cheetah.** `#if any(x for $b in $barcodes)` fails. Use an explicit `#for` loop that sets a flag.
- **`find ... -exec ... \;` gets mangled in a command block.** Use a shell glob or a `for f in dir/*; do ...; done` loop instead (which is also the remote-execution-safe way to enumerate files).
- **`re` is not available until you import it.** `#import re` first, or `re.sub(...)` raises `NameError` (seen in a configfile).

---

## 5. Input Parameter Conventions

**Standard:** [Parameter name, argument and help][std-argument], [Booleans][std-booleans], [Dynamic Options][std-dynamic], [Data parameters][std-data-params], [security item 2][sec-2].

### The `argument=` Attribute

Prefer `argument="--flag"` over bare `name="flag"`. This auto-generates the `name` attribute (stripping leading dashes, replacing `-` with `_`) and displays the flag in the help text. Always use the long form (`--output` not `-o`).

```xml
<!-- GOOD: argument auto-generates name="min_score" and shows --min-score in help -->
<param argument="--min-score" type="float" value="0.5" min="0.0" max="1.0"
       label="Minimum score" help="Filter results below this threshold"/>

<!-- Also fine when there's no direct CLI flag mapping -->
<param name="organism" type="select" label="Organism">
```

`argument=` does not inject the flag. You still write `--min-score $min_score` in the command block.

**IUC standard attribute order:** `name, argument, type, format, min|truevalue, max|falsevalue, value|checked, optional, label, help`.

### Grouping with `<section>`

Group related parameters into logical sections for complex tools:

```xml
<inputs>
    <param argument="--input" type="data" format="vcf" label="Input VCF"/>

    <section name="filtering" title="Filtering options" expanded="false">
        <param argument="--min-score" type="float" value="0.5" .../>
        <param argument="--max-pvalue" type="float" value="0.05" .../>
    </section>

    <section name="output_options" title="Output options" expanded="false">
        <param argument="--format" type="select" .../>
    </section>
</inputs>
```

Params in a collapsed `<section>` still have values and still reach the command. Params in a non-selected `<when>` of a `<conditional>` do not ([Dynamic Options][std-dynamic]).

### Validation

**Standard:** [security item 2][sec-2].

Order of preference: a typed parameter (`integer`, `float`, `boolean`, `select`, `data_column`) over free text; a `select` whenever the valid values are a finite set; free text only when nothing else describes the input, and then with a validator.

Tripwire: Galaxy regex validators match from the start only. **Anchor the end** (`^...$`) or the validator accepts `good;anything`. An unanchored pattern with `*` accepts every string.

Validator vs sanitizer: a validator rejects a bad value before the job starts and tells the user why; a sanitizer silently rewrites characters as the value is rendered. Prefer the validator. The default sanitizer is not shell escaping, and a custom `<sanitizer>` that lets `'` through breaks the `'$value'` quoting boundary. Validation complements quoting; it does not replace it.

Use `min`/`max` attributes on integer/float params. Never use `optional="true"` when a default is appropriate -- just set the default.

```xml
<!-- Multi-select: require at least one selection -->
<param name="items" type="select" multiple="true" label="Items">
    <option value="A" selected="true">A</option>
    <option value="B">B</option>
    <validator type="no_options" message="Select at least one item"/>
</param>

<!-- Text with regex validation: anchored at both ends -->
<param name="terms" type="text" value="" label="Terms (optional)"
    help="Comma-separated terms, e.g. UBERON:0002107">
    <validator type="regex" message="Only alphanumeric, colons, commas, spaces">^[A-Za-z0-9:, ]*$</validator>
</param>

<!-- Integer with range -->
<param name="max_count" type="integer" value="100" min="1" max="10000"
       label="Maximum items" help="Start small to verify results"/>
```

Every free-text param gets a hostile-input test (`references/testing.md`).

### Boolean Parameters

Put the CLI flag in `truevalue`/`falsevalue`. Don't use booleans as conditionals — use `select` instead when other params depend on the choice.

```xml
<!-- GOOD: flag in truevalue -->
<param argument="--gzip" type="boolean" truevalue="--gzip" falsevalue=""
       checked="false" label="Compress output"/>

<!-- BAD: boolean controlling other params — use select + conditional instead -->
<param name="advanced" type="boolean" label="Show advanced options"/>
```

### Compressed Datatype Support

**Standard:** [Data parameters][std-data-params].

If Galaxy has a compressed version of a datatype (`fasta.gz`, `fastqsanger.bz2`), the data param accepts both the compressed and the uncompressed format. If the underlying tool cannot read the compressed form, decompress in `<command>`. Galaxy does not decompress for you. Reviewers ask "are fasta.gz supported?" as a matter of course.

```xml
<param name="input" type="data" format="fasta,fasta.gz" label="Input sequences"/>
<param name="reads" type="data" format="fastqsanger,fastqsanger.gz,fastqsanger.bz2" label="Reads"/>
```

```
#if $input.is_of_type('fasta.gz')
    gunzip -c '$input' > input.fa &&
#else
    ln -s '$input' input.fa &&
#end if
mytool input.fa
```

Include tests with both compressed and uncompressed inputs to verify both paths work.

### Multiple Inputs vs Repeat Parameters

For multiple files of the same type, use `multiple="true"` on the data param:

```xml
<param name="inputs" type="data" format="bam" multiple="true" label="BAM files"/>
```

Use `<repeat>` for variable-length groups of mixed parameters (e.g., a dataset + options per entry):

```xml
<repeat name="samples" title="Sample">
    <param name="input" type="data" format="bam" label="BAM file"/>
    <param name="label" type="text" label="Sample label"/>
</repeat>
```

In command: `#for $s in $samples# --sample '$s.input' --label '$s.label' #end for#`

### Data Collections

For paired-end data and multi-sample workflows:

```xml
<!-- Accept paired collection -->
<param name="paired_input" type="data_collection" collection_type="paired"
       format="fastqsanger,fastqsanger.gz" label="Paired reads"/>

<!-- Accept multiple datasets. Users can pick a list collection for this in the UI,
     but a test <collection> block cannot feed it: tests pass value="a.bam,b.bam". -->
<param name="inputs" type="data" format="bam" multiple="true" label="BAM files"/>
```

Access paired ends: `$paired_input.forward` / `$paired_input.reverse`

### Identifiers as Filenames

**Standard:** [security item 3][sec-3].

`element_identifier`, dataset names, and labels are user-influenced strings. Tripwire, in order of preference:

1. Use a **fixed filename** (`ln -s '$input' input.bam`) unless the tool needs the identifier in the name (for example because it becomes a sample name in the output).
2. If it does, import `re`, convert with `str()`, map to a small alphabet, and **handle collisions**: two identifiers can sanitize to the same name. Add a stable unique suffix or fail.
3. Single-quote the result wherever it reaches the shell.

```
#import re
#for $i, $input in enumerate($inputs)
    #set $safe_id = re.sub(r'[^\w.-]', '_', str($input.element_identifier))
    ln -s '$input' '${safe_id}_${i}.bam' &&
#end for
```

The loop index is the unique suffix here. If the suffix would leak into output the user sees, detect duplicates and `exit 1` with a clear message instead.

### Reference Data and Data Tables

**Standard:** [integration_checklist][std-checklist] ("built-in data", "tool-data" bullets); [pulsar items 3 and 4][pulsar-3].

Reference data (indexes, databases) comes from a data table, never a hardcoded path and never `$__app__.tool_data_tables`:

```xml
<param name="reference" type="select" label="Reference genome">
    <options from_data_table="mytool_indexes">
        <validator type="no_options" message="No indexes are available"/>
    </options>
</param>
```

In the command, use `'$reference.fields.path'` directly. Do not copy it through `str()` into a `#set` first: that bakes a server-side path into the command that remote execution cannot rewrite.

Tripwire: a tool that uses a data table ships **four** files, and missing ones are a routine review finding.

| File | Content |
|---|---|
| `tool-data/mytool_indexes.loc.sample` | Comment-only: column layout and an example row, all `#`-prefixed |
| `tool_data_table_conf.xml.sample` | The `<table>` definition pointing at `tool-data/mytool_indexes.loc` |
| `tool_data_table_conf.xml.test` | Same table, with `<file path="${__HERE__}/test-data/mytool_indexes.loc"/>` |
| `test-data/mytool_indexes.loc` | Real tab-separated rows for the tests; paths use `${__HERE__}/...` |

> **Verify:** the checklist calls `tool_data_table_conf.xml.sample` "comment-only" as well, but existing tools-iuc tools (for example `tools/bowtie2`) ship it with a live `<table>` definition. This skill follows the tools-iuc form. Copy a sibling tool's files rather than writing from scratch.

Test setup is in `references/testing.md`. A tool that *populates* a data table is a data manager: see [data_managers][std-dm].

### Subcommand Strategy

**Standard:** [Wrapping tools with subcommands][std-subcommands]. Separate wrappers when subcommands need different resources; a conditional only when they are closely related and share a resource profile.

### IUC-Specific Rules

Reviewer habits and style points. The first, third, and fifth are skill-side observations from review, not standards text.

- **No `display="checkboxes"`** on multi-select params. Reviewers will flag it. Let Galaxy decide the widget.
- **Mark defaults with `selected="true"`** on the `<option>` element.
- **Unit notation in labels:** Use SI-style lowercase: "kb" not "KB", "Mb" not "MB". The option values can use uppercase (e.g., `value="1MB"`) but display text should read "1 Mb".
- **Help text:** Keep it actionable and short. "Start small to verify results" is good.
- **Citations:** Prefer `type="doi"` over `type="bibtex"` when a DOI is available. Search the upstream repo/paper for the correct DOI.
- **4-space indentation** throughout XML and Cheetah code. Run `planemo format` before submitting to ensure consistent indentation matching Galaxy Language Server style.
- **One line per element:** keep all attributes of an element on one line; only `label` and `help` may wrap ([Coding Style][std-style]).
- **New datatypes** go to the Galaxy repository, not the tool repo ([tool_xml][std-tool-xml], "Datatypes").

---

## 6. Python Wrapper Scripts (When Needed)

**Standard:** [Coding Style][std-style] (Python 3, PEP 8, smarkets import order); [integration_checklist][std-checklist] (flake8 in IUC CI).

Add a wrapper script only when:

- The CLI output format doesn't match what Galaxy expects (needs conversion)
- You need multi-step pipelines that are too complex for shell in the command block
- The tool wraps an API rather than a CLI binary
- You need per-item error handling with partial results

### When You Do Need a Script

Key rules for Galaxy wrapper scripts:

- **All logging to stderr** — stdout is reserved for Galaxy output capture
- **`__version__`** must match `@TOOL_VERSION@` in macros.xml
- **Per-item error handling** — `try/except` inside the loop, continue on failure, exit 1 only if ALL items failed
- **No user-facing `--verbose`** — control via Galaxy's logging, not user params
- **Declare the script** (and anything it imports or sources) in `<required_files>` -- see Section 2
- **Pin a python requirement.** Requirements are referenced by name and version ([package_xml][std-package]); reviewers pin python in practice. Use the version sibling tools use:
  ```xml
  <requirement type="package" version="3.12">python</requirement>
  ```
  A second `<requirement>` makes this a multi-requirement tool, which needs a mulled container (Section 9).
- **`csv.writer` emits CRLF by default.** Pass `lineterminator="\n"`, or the output breaks `$`-anchored test assertions and downstream tabular tools.
- **flake8 clean.** IUC CI lints Python with flake8; run it before hand-back.

```python
#!/usr/bin/env python
"""MyTool converter for Galaxy — transforms native output to tabular."""

import argparse
import logging
import sys

__version__ = "1.2.3"  # Must match @TOOL_VERSION@


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.StreamHandler(sys.stderr)],
    )

    # Process input, write output...
    logging.info("Processed %s -> %s", args.input, args.output)


if __name__ == "__main__":
    main()
```

### API Tool Scripts

**Standard:** [security item 8][sec-8] (credentials), [item 9][sec-9] (downloads).

For tools wrapping external APIs (not the common case), the script also handles:
- **Fixture bypass** — `--test-fixture` loads JSON instead of calling the API
- **Credentials** -- read API key from environment (`os.environ.get("MYTOOL_API_KEY")`); never log it, never accept it as an argument
- **Network access** -- HTTPS with certificate verification on, pinned or versioned endpoints, checksums for downloaded artifacts, and no unconstrained user-supplied URL passed to `curl`, `wget`, or an HTTP client

---

## 7. Resource Management

Skill-only guidance; the standards do not cover it.

### General Principle

Do not expose `--verbose`, `--memory`, `--threads`, or similar operational flags as user-facing inputs. Set sensible defaults or derive them from Galaxy's system variables in the command block. Users control what to compute, admins control how.

### Parallelism (GALAXY_SLOTS)

Use Galaxy's `GALAXY_SLOTS` environment variable for thread/process count. This is admin-controlled, not user-facing.

Available environment variables:

| Variable | Meaning |
|----------|---------|
| `\${GALAXY_SLOTS:-1}` | Allocated CPU cores |
| `\${GALAXY_MEMORY_MB:-4096}` | Total memory in MB |
| `\${GALAXY_MEMORY_MB_PER_SLOT:-4096}` | Memory per core in MB |

In the command block:

```
--threads \${GALAXY_SLOTS:-1}
```

### Memory

If the upstream CLI has no hard memory limit, surface the parameter that proxies it best (chunk size, hash table size, block size) and set conservative defaults. Don't expose raw memory flags.

```xml
<!-- GOOD: proxy for memory usage, meaningful to the user -->
<param argument="--chunk-size" type="integer" value="10000" min="1000" max="1000000"
       label="Chunk size" help="Larger values use more memory but process faster"/>

<!-- BAD: raw memory flag, meaningless to most users -->
<param argument="--memory" type="integer" value="4" label="Memory (GB)"/>
```

Scratch space is the job working directory. No `/tmp` and no other absolute paths ([pulsar item 4][pulsar-4]).

---

## 8. Test Infrastructure

**Standard:** [Tests][std-tests]; [integration_checklist][std-checklist] "Testing Your Tool".

See `references/testing.md` for the full assertion reference, collection testing, compressed output testing, repeat element tests, tool-data test setup, the hostile-input test, gotchas, and failure analysis. Also see the [Galaxy Tool XSD Schema](https://docs.galaxyproject.org/en/latest/dev/schema.html) for the complete assertion specification.

### Test XML Structure

```xml
<tests>
    <test expect_num_outputs="1">
        <param name="input_fastq" value="test_input.fastq.gz"/>
        <param name="min_quality" value="20"/>
        <output name="output_bam" file="expected_output.bam" ftype="bam" compare="sim_size" delta="100"/>
    </test>
    <test expect_num_outputs="1">
        <!-- Second test: different code path -->
        <param name="input_fastq" value="test_input.fastq.gz"/>
        <param name="min_quality" value="0"/>
        <output name="output_bam">
            <assert_contents>
                <has_size min="100"/>
            </assert_contents>
        </output>
    </test>
</tests>
```

### Test Structure Rules

- Give each test a **unique purpose** (e.g., "defaults", "compression on", "filtering active")
- Point expected files to unique golden files — no duplicate outputs for different logic paths
- Always include `expect_num_outputs` — count only outputs actually produced (filters that evaluate to False don't produce output)
- **Test data: aim for under 1 MB in total.** Shrink or downsample inputs first. For large *outputs*, assert on content (`<has_text>`, `<has_size>`, `compare="contains"`) instead of committing a golden file. For a large *input* that cannot be made smaller, deposit it in the [Galaxy Tool Developers Zenodo community](https://zenodo.org/communities/galaxy_tool_dev_sig) and reference it with `location=` (see below). Assertions do not solve large inputs.
- Include tests for output filters to verify filtering actually occurs
- Include tests for error conditions with expected failure
- Include a hostile-input test for every free-text parameter (`references/testing.md`)

Large input from Zenodo:

```xml
<param name="input_file" location="https://zenodo.org/records/RECORD/files/FILE_NAME"/>
```

### Assert Patterns

Use `<assert_contents>` with `<has_text>` for key content verification. Check header lines, known output values, and identifiers from the input data:

```xml
<output name="output_tsv">
    <assert_contents>
        <has_text text="chrom"/>           <!-- column header -->
        <has_text text="mean_signal"/>     <!-- column header -->
        <has_text text="sample_1"/>        <!-- from input -->
    </assert_contents>
</output>
```

For the full assertion catalog (file comparison modes, stream/command assertions, compressed output testing), see [`references/testing.md`](references/testing.md).

Test conditionals and sections with explicit nesting (not pipe syntax):

```xml
<!-- BAD: pipe syntax -->
<param name="filtering|min_score" value="0.5"/>

<!-- GOOD: explicit nesting -->
<section name="filtering">
    <param name="min_score" value="0.5"/>
</section>
```

### Running Tests

Run `planemo test --biocontainers` to execute tests. See **Reference: Useful Planemo Commands** at the end of this document for the full command set. If the run fails before any test starts, check the Gotchas in `references/testing.md` (Galaxy Python version, missing mulled image).

### Generating Expected Output Files

The standard approach for CLI-wrapping tools. Run the tool once via planemo and let it update the expected output files in place:

```bash
planemo test --biocontainers --update_test_data tools/mytool/mytool_align.xml
```

Read the generated files before committing them. `--update_test_data` records whatever the tool produced, including wrong output.

### Fixture-Based Testing (API Tools Only)

For the uncommon case of tools calling external APIs: record real responses as JSON fixtures and replay them in tests. This lets `planemo test` run without API keys. See the `test_fixture_param` macro in the macros.xml Patterns section above.

**Scrub recorded fixtures before committing.** A recorded response or request can carry API keys, bearer tokens, session cookies, signed URLs, account identifiers, or private endpoints. Record the response body only (not request headers), replace anything secret-shaped with an obvious placeholder, and grep the fixture for the key you used ([security item 8][sec-8]).

```bash
# Fixture-based tests (no container needed, no API key)
planemo test tools/mytool/
```

---

## 9. .shed.yml, Dependencies, and Versioning

### .shed.yml

**Standard:** [shed_yml][std-shed]; auto-generated suites: [Advanced Parameters][std-shed-advanced].

```yaml
categories:
- Relevant Category
- Another Category
description: Short one-line description
long_description: |
  Multi-line description of the tool suite.
  What it does, what API/library it wraps.
name: mytool
owner: iuc
homepage_url: https://github.com/original/project
remote_repository_url: https://github.com/galaxyproject/tools-iuc/tree/main/tools/mytool
type: unrestricted
auto_tool_repositories:
  name_template: "{{ tool_id }}"
  description_template: "Wrapper for MyTool application: {{ tool_name }}."
suite:
  name: "suite_mytool"
  description: "A suite of Galaxy tools for MyTool."
  type: repository_suite_definition
```

**Key points:**
- `auto_tool_repositories` creates one Tool Shed repo per tool XML automatically
- `suite` groups all tools under a single installable suite. Let planemo generate it; do not hand-write a suite for one set of related tools
- `owner: iuc` for tools-iuc submissions
- `remote_repository_url` should point to the tools-iuc path (the final location)
- Repository names cannot contain periods; a tool repo is `$name`, a suite `suite_$name`
- `include` / `exclude` for shared or oversized files: see [shed_yml][std-shed]
- Check the file with `planemo shed_lint tools/mytool/`. Do not run any `shed_create` / `shed_upload` / `shed_update` command: IUC CI publishes after merge

### Dependencies

**Standard:** [package_xml][std-package].

- Requirements must be available in **conda-forge** or **bioconda** channels, referenced by name **and version**. Use the latest stable versions.
- A tool with one requirement gets its BioContainer automatically. A tool with **two or more** requirements needs a mulled container for that exact combination, which is only built once the combination is registered ([Multi-requirement containers][std-package-mulled]). For repositories tracked by planemo-monitor (tools-iuc is the usual case) that registration happens automatically once the tool is merged, so `planemo test --biocontainers` can fail locally on a new combination with no image on quay. Build it locally (see Gotchas in `references/testing.md`) and say so in the hand-back.

### Versioning and Lint Hygiene

- The version-bump rule is in Section 11.
- If a `.lint_skip` file exists in the tool directory, try to fix the underlying issues and remove the skip entries rather than adding to them.
- Run `planemo lint` after every change. It catches element ordering, missing attributes, and parameter mismatches.

---

## 10. Pre-PR Passes

Four passes, run in order by Step 8 of Section 12. For every item record **pass**, **fail**, or **N/A** with a one-line reason. "N/A" needs the reason as much as "fail" does (for example "N/A -- tool makes no network calls"). The rationale for each item lives in the standards; the tables below give only the check.

### 10a. Security Pass

**Standard:** [Tool Security Checklist][std-security]. The readiness checklist asks you to trace every user-controlled value through Cheetah, the generated command, helper scripts, files, network clients, and rendered outputs. These ten checks are that trace.

The greps produce candidates to read, not verdicts. They cover `macros.xml` too, since tokens expand into the command. Run them in bash (`$F` is left unquoted on purpose so the glob expands):

```bash
F="tools/mytool/*.xml"; D="tools/mytool"
grep -nE '\$\{?[A-Za-z_]' $F                                               # g1a: every $ reference
grep -nE 'bash -c|sh -c|awk|sed|perl -e|Rscript -e|python -c' $F           # g1b: nested interpreters
grep -nE 'type="text"|<sanitizer' $F                                       # g2: free text, sanitizers
grep -nE 'element_identifier|\.name\b|display_name' $F                     # g3: identifiers
grep -nE '#(set|echo)|\.(split|strip|replace|join|format)\(|sanitize=' $F  # g4: expressions on params
grep -nE '#import (os|subprocess|pathlib|glob|shutil)|open\(|os\.|__app__' $F   # g5: template side effects
grep -nE '<configfile|eval|exec\(' $F                                      # g6: generated code
grep -rniE 'pickle|joblib|rds|rdata|torch|tar |unzip' $D --include='*.xml' --include='*.py' --include='*.R'   # g7
grep -rniE 'api[_-]?key|token|secret|bearer|authorization|password' $D     # g8: secrets, fixtures included
grep -rnE  'curl|wget|https?://|requests\.|urllib' $D --include='*.xml' --include='*.py' --include='*.R'      # g9
grep -nE  'format="(html|svg)"' $F                                         # g10: active content
```

| # | Item | Check |
|---|---|---|
| 1 | [Single-quote dynamic values in shell commands][sec-1] | g1a: each reference in `<command>` and tokens is inside single quotes, or is an `integer`/`float`, a boolean `truevalue`, or a fixed-option `select`. No `"$param"`. g1b: no parameter is interpolated into the inner program. |
| 2 | [Constrain free text to the smallest useful language][sec-2] | g2: for each text param, could it be a typed param or a `select`? If not, it has a regex validator anchored at both ends, or a written reason. No sanitizer admits `'`. A hostile-input test exists. |
| 3 | [Treat identifiers and generated filenames as untrusted][sec-3] | g3: each use is not a filename, or was replaced by a fixed filename, or goes through `#import re` + `str()` + `re.sub` + collision handling. Every one is single-quoted at the shell. |
| 4 | [Preserve sanitization through Cheetah expressions][sec-4] | g4: every method call on a text param is on `str($param)`, and each result that reaches the shell is single-quoted. No `sanitize="false"` without a validator and a reason. |
| 5 | [Keep Cheetah templating side-effect free][sec-5] | g5 returns nothing in `<command>`, tokens, or `<configfile>` Cheetah directives, or only hits you can explain (plain script text in a configfile runs in the job, not in the template). |
| 6 | [Keep user values as data, never generated code][sec-6] | g6: for each `<configfile>`, is its content code (R, Python, shell, awk)? If so, no `$param` appears in it; values come from `<configfiles><inputs>` JSON or arguments. No `eval`/`exec`, no script-body param, no free-text "extra options" param. |
| 7 | [Treat uploaded serialized objects and archives as active input][sec-7] | N/A unless g7 or the input formats show pickle, RData/RDS, a model checkpoint, or an archive the tool extracts. Otherwise: data-only loading mode or a documented trust decision; extraction rejects absolute and `..` paths. |
| 8 | [Keep credentials out of commands, arguments, and logs][sec-8] | N/A unless the tool uses a secret. Otherwise: secrets come only from `<credentials>` into the environment; none in XML, command line, `<configfile>`, `echo`, logs, or outputs. g8 shows placeholders only. |
| 9 | [Verify downloads and preserve transport security][sec-9] | N/A unless g9 shows the tool or its script fetching anything. Otherwise: HTTPS with verification on, versioned URLs, checksum verified, job fails on mismatch, no unconstrained user URL. |
| 10 | [Treat browser-rendered output as active content][sec-10] | N/A unless g10 or the script shows an HTML, SVG, or JavaScript output. Otherwise: user-controlled values are escaped for their context and the trust decision is written in the hand-back. |

### 10b. ToolShed Readiness Pass

**Standard:** [ToolShed Readiness Checklist][std-checklist].

| # | Item | Check |
|---|---|---|
| 1 | Licensing | The wrapped software's license is OSI-approved and you can name it. No dependency, model, or reference dataset the tool needs carries a use restriction, or the restriction is handled per [licensing][std-licensing]. The PR template's "license permits unrestricted use" box can be ticked honestly. |
| 2 | ID, name, version, description | Follow [Tool ids][std-ids], [Tool versions][std-versions]; `<description>` present. Version is built from tokens. |
| 3 | Requirements | Every `<requirement>` has a `version`. Two or more requirements: mulled container noted (Section 9). |
| 4 | Version command | `<version_command>` runs the tool's real version command. |
| 5 | Security review | Pass 10a is complete. |
| 6 | Outputs | Optional outputs use `<filter>`; each filter has a test where it is true and one where it is false. |
| 7 | Tests | At least one test per code path, each with `expect_num_outputs`; at least one expected-failure test. |
| 8 | Help | `<help>` exists, is in CDATA, and is valid in the one format it declares (Section 2). The checklist at the pinned commit says reStructuredText; Markdown rests on [PR #95][std-pr95] and review practice. |
| 9 | Citation | `<citations>` present, DOI preferred. |
| 10 | Built-in data | N/A unless a data table is used. Otherwise all four files from Section 5 exist: `ls tool-data/*.loc.sample tool_data_table_conf.xml.sample tool_data_table_conf.xml.test test-data/*.loc`. |
| 11 | Order and indentation | Elements and param attributes in the [Coding Style][std-style] order; 4 spaces. `planemo format` then `planemo lint` clean. |
| 12 | Test data size | `du -sh tools/mytool/test-data`. Over 1 MB: shrunk, replaced by assertions, or moved to Zenodo with `location=`. State the size. |
| 13 | `planemo lint` | Passes with no errors; warnings are fixed or explained. |
| 14 | `planemo test` | All tests pass. Quote the summary line. |
| 15 | `planemo serve` | Tool form looks right when served locally. If you cannot open a browser, record N/A with that reason and ask the human to check. |
| 16 | `.shed.yml` | Present; `planemo shed_lint` passes. |
| 17 | Python lint | N/A with no Python. Otherwise `flake8` passes on every shipped `.py`. |
| 18 | Tracking issue | A tools-iuc issue for the tool exists or was checked for (`references/tool-placement.md`). |

### 10c. Remote Execution Pass (Pulsar)

**Standard:** [Remote Execution Checklist][std-pulsar]. These fail only on a production server, so no local test catches them.

| # | Item | Check |
|---|---|---|
| 1 | [Declare every tool-directory file with `<required_files>`][pulsar-1] | r1 over `<command>`, `<version_command>`, `<configfiles>`, and tokens: every referenced file has an `<include>` with the same relative path, including the entry script and files reached by `cp`, `source()`, or `import`. No include points at a file that is never referenced. Nothing is read from `$__tool_directory__/test-data/` at runtime. |
| 2 | [No filesystem I/O in Cheetah templates][pulsar-2] | g5 from 10a is clean. File sizes and types come from dataset metadata; directory listings happen in the shell or in script text at runtime. |
| 3 | [Do not reach into `$__app__` at template time][pulsar-3] | g5 from 10a shows no `__app__`. Reference data uses `from_data_table` and `$param.fields.path` directly, not through a `#set`. |
| 4 | [No hardcoded absolute paths][pulsar-4] | r4 over command, tokens, configfiles, and scripts returns nothing unexpected. No `/tmp`. |
| 5 | [Outputs must be real files inside the job output directory][pulsar-5] | r5: no `ln -s` creates an `$output`, an entry in an output's `files_path`, or a file later collected as an output. Outputs are written, moved, or copied. |
| 6 | [Keep discovered outputs inside the job working directory][pulsar-6] | r5: every `from_work_dir` and `discover_datasets directory=` is a relative path under the working directory. |

```bash
grep -n   '__tool_directory__' $F                                              # r1
grep -rnE '/(tmp|data|mnt|opt|usr|home|scratch)/' $D --include='*.xml' --include='*.py' --include='*.R'   # r4
grep -nE  'ln -s|from_work_dir|discover_datasets' $F                           # r5
```

### 10d. Review Flags

Reviewer habits and wrapper conventions not already covered by 10a-10c. Rows that became items above (help and citation present, test data size, text validators, credentials) were removed from these tables.

#### Will Definitely Be Flagged

| Issue | What Reviewers Say | Fix |
|-------|-------------------|-----|
| Cheetah snippet in `<xml>` macro | "This should be a `<token>`, not `<xml>`" | Move to `<token name="@CMD_...@">` with CDATA |
| `display="checkboxes"` on multi-select | "Remove display attribute" | Delete it; let Galaxy pick the widget |
| Missing `detect_errors` | "Add error detection" / "Why not aggressive?" | Add `detect_errors="aggressive"`; `exit_code` only for the two cases in Section 2, with a comment |
| Version not from tokens | "Use @TOOL_VERSION@ token" | Replace hardcoded version strings |
| Profile hardcoded per tool in a suite | "Please use a macro token for the profile" | `@PROFILE@` token in `macros.xml` |
| `version_command` that echoes a constant | "Use the tool's version command" | `<tool> --version` |
| `KB` / `MB` in display text | "Use kb / Mb (SI lowercase)" | Fix option labels to lowercase units |
| Missing `expect_num_outputs` on test | "Add expect_num_outputs" | Add `expect_num_outputs="1"` to `<test>` |
| Test data outside `test-data/` | "Move to test-data directory" | Move files, update paths |
| `optional="true"` with a default | "Just use the default, remove optional" | Remove `optional`, set `value` |
| stdout used for logging | "Use stderr for logging" | `logging.StreamHandler(sys.stderr)` |
| Missing `argument=` on params | "Use argument= instead of bare name=" | `<param argument="--flag" .../>` |
| Bare `name=` duplicating `argument=` | "argument= auto-generates the name" | Remove redundant `name=`, keep `argument=` (flag still needed in command) |
| Param not used in command | "Orphaned parameter" | Remove param or wire it into command |
| Boolean used as conditional | "Use select + conditional" | Replace boolean with select param when other params depend on choice |
| Boolean wrapped in `#if` just to emit a flag | "Use truevalue and falsevalue and simplify the command" | Flag in `truevalue`, emit `$param` directly |
| Missing compressed format support | "Are fasta.gz supported?" | `format="fasta,fasta.gz"`, decompress in `<command>` if the tool cannot |
| Missing bio.tools xref | "Add bio.tools cross-reference" | Add `<xrefs><xref type="bio.tools">id</xref></xrefs>` |
| Unpinned `python` requirement | "we pin python" | Add `version=` |

#### Commonly Requested Improvements

| Improvement | Typical Comment | Resolution |
|-------------|----------------|------------|
| Validator on multi-select | "Add no_options validator" | `<validator type="no_options" message="..."/>` |
| Help text structure | "Break this into sections" | Headings in Markdown help; titles or `-----` rules in RST help |
| Help format | "You can use markdown if this is any easier" | New tool without images: `<help format="markdown">`, profile >= 24.2. Existing RST help: leave it |
| Shared XML duplicated across a suite | "This can go in macros.xml" | Move to a macro; do not add a macro for something used once |
| Consistent output labels | "Use standard label pattern" | `${tool.name} on ${on_string}` |
| Test assert_contents | "Add content assertions, not just expect_num_outputs" | Add `<has_text>` checks |
| Error handling in script | "What happens if one item fails?" | Per-item try/except with continue (wrapper scripts) |
| Fixture-based tests (API tools) | "Tests shouldn't need network" | Add `--test-fixture` hidden param pattern |
| Section grouping | "Group related params" | Use `<section name="..." title="...">` |
| Pipe syntax in tests | "Use explicit nesting" | Replace `section\|param` with nested XML |
| `planemo tool_init` scaffold | "Use tool_init for boilerplate" | `planemo tool_init --id ... --requirement ...` |
| Test output `lines_diff` | "Use lines_diff for non-deterministic outputs" | `<output ... lines_diff="2"/>` |
| Missing stderr/stdout assertions | "Assert on expected warnings" | `<assert_stderr><has_text .../></assert_stderr>` |
| `<creator>` | "You can add yourself as creator" | Optional; add when the human asks (Section 2) |

---

## 11. Updating Existing Tools

**Standard:** [Tool versions][std-versions].

When updating a tool to a new upstream version, follow this workflow.

### When to Bump the Version

| Change | `@TOOL_VERSION@` | `@VERSION_SUFFIX@` |
|---|---|---|
| New upstream version | New version | Reset to `0` |
| Wrapper-only change that will be released (command, params, outputs, requirements, a fix) | Unchanged | **+1** |
| Cosmetic only (help text, labels, whitespace) | Unchanged | Unchanged |

The standards require the bump for every wrapper change released to the Tool Shed, so that two different wrappers never share a version. Versions follow PEP 440. In a suite with one shared `@VERSION_SUFFIX@`, bumping it bumps every tool; follow the suite's existing practice (the standards allow a per-tool `@TOOL_VERSION@+galaxyN` instead).

### Research Upstream Changes

Before touching the XML, check what changed:

1. **Check release notes** — GitHub releases, CHANGELOG.md, migration guides
2. **Compare --help output** between old and new versions (run the container if available)
3. **Search for breaking changes** — renamed/removed flags, changed defaults, output format changes

```bash
docker run quay.io/biocontainers/<package>:<new_version> <command> --help
```

| Change Type | Impact | Action |
|-------------|--------|--------|
| Flag renamed | Tool breaks | Update command section |
| Flag removed | Tool breaks | Remove from XML or make conditional |
| New required flag | Tool breaks | Add to command section |
| New optional flag | None | Consider adding to inputs |
| Default changed | Output changes | Update tests, maybe document |
| Output format changed | Tests break | Update golden files and assertions |

### Version Bump Procedure

1. Update `@TOOL_VERSION@` (in `macros.xml`, or inline in a single-tool wrapper)
2. Reset `@VERSION_SUFFIX@` to `0`
3. Review command section for deprecated/renamed flags
4. Review output filters for correctness
5. Update help section if upstream added significant features. Keep the help format the tool already uses
6. Run `planemo lint` and `planemo test`
7. Fix test failures — update golden files, adjust assertions
8. Re-run the passes in Section 10 for everything the change touched
9. Commit with descriptive message

### Common Update Bugs

**Repeat element access** — the most common bug found during updates:

```cheetah
## WRONG: accesses outer scope, not loop variable
#for item in $filters.search:
    --search '$filters.search_term'
#end for

## CORRECT: access param through loop variable
#for item in $filters.search:
    --search '$item.search_term'
#end for
```

**Wrong filter value** — output filters checking the wrong condition:

```xml
<!-- BUG: 3' UTR output filtering on 5p-utr -->
<data name="threep_utr" format="fasta">
    <filter>"5p-utr" in file_choices['include']</filter>
</data>

<!-- CORRECT -->
<data name="threep_utr" format="fasta">
    <filter>"3p-utr" in file_choices['include']</filter>
</data>
```

---

## 12. Step-by-Step Workflow: Creating a New Tool

### Step 0: Go / No-Go

Two gates before any XML exists.

**Existing wrappers.** Search the [Main Tool Shed](https://toolshed.g2.bx.psu.edu/) and [Test Tool Shed](https://testtoolshed.g2.bx.psu.edu/) web UIs, and search `tools-iuc`, `tools-devteam`, and other IUC-maintained repos on GitHub for existing wrappers.

Also check `bioconda` and `conda-forge` for an existing recipe — if none exists, you may need to create one first.

**Licensing.** **Standard:** [licensing][std-licensing], items [1][lic-osi] to [3][lic-redistribution]. The IUC accepts wrappers for software under an OSI-approved license, and the PR template asks you to confirm the license permits unrestricted use. Before writing anything:

- Find the license of the wrapped software and confirm it is OSI-approved. If it is not, stop and ask the human; wrapper engineering does not fix this.
- Check the other layers: dependencies, and any reference data or models the tool needs. An open-source tool can depend on a restricted model or database.
- Confirm redistribution is permitted (normally: a Bioconda package exists or can exist).
- The `license` attribute on `<tool>` describes the wrapper, not the wrapped software. Do not use it to record the upstream license.

If any layer is restricted, follow the licensing page rather than improvising; it has the declaration and affirmation pattern. Record the outcome in the hand-back either way.

### Step 1: Inspect the Upstream Tool

```bash
docker run quay.io/biocontainers/<package>:<version> <command> --help
```

Identify: required inputs, outputs, key parameters, threading flags, version string. Capture the stderr of a successful run too; it decides `detect_errors` (Section 2).

### Step 2: Set Up Structure

Use `planemo tool_init` for the initial scaffold, then refine:

```bash
planemo tool_init --id mytool_function --name "MyTool Function" \
    --description "Brief description" \
    --requirement mytool@1.2.3 \
    --example_command "mytool function --input input.fa --output output.tsv" \
    --example_input input.fa --example_output output.tsv \
    --doi "10.1234/example" \
    tools/mytool/
```

Or manually:

```bash
mkdir -p tools/mytool/test-data
touch tools/mytool/.shed.yml
```

Look at two or three sibling tools in the target repo now and note their conventions: help format, `detect_errors`, profile, inline tokens vs `macros.xml`. Those win over this skill (see Authority).

### Step 3: Write the Tokens (and macros.xml if Shared)

One tool with nothing shared: put `@TOOL_VERSION@` and `@VERSION_SUFFIX@` tokens inline in the tool's `<macros>` block and skip the file. Two or more tools, or anything reused: create `macros.xml` with version tokens, `@PROFILE@`, requirements (from bioconda/conda-forge), shared params, and citations, using the template in Section 3.

### Step 4: Write the Tool XML

Follow element order from Section 2:
- `<tool>` with token-based version and a profile chosen by the rule in Section 2
- `<macros>` import or inline tokens
- `<xrefs>` bio.tools
- `<expand macro="requirements">` or inline `<requirements>`, every requirement pinned
- `<required_files>` for every wrapper-shipped file
- `<version_command>` -- the tool's real version command; CDATA if it contains `&`
- `<command>` with CDATA, `detect_errors`, and direct CLI call; everything single-quoted
- `<inputs>` -- use `argument=` for CLI flags, group with `<section>`, typed params first, anchored validators on text
- `<outputs>` with `${tool.name} on ${on_string}` labels, `from_work_dir` if needed
- `<tests>` — at least one test per code path, with `expect_num_outputs` and `assert_contents`
- `<help>` in the format chosen by the table in Section 2
- `<expand macro="citations"/>` or inline `<citations>`
- `<creator>` only if asked

### Step 5: Create Test Data and Golden Files

- Create minimal input files that exercise the tool
- Run `planemo test --biocontainers --update_test_data` to generate expected output files in place, then read them
- For a data table, add the four files from Section 5

### Step 6: Write .shed.yml

Use `auto_tool_repositories` and `suite` for multi-tool packages. Check it:

```bash
planemo shed_lint tools/mytool/
```

### Step 7: Format, Lint, and Test

```bash
planemo format tools/mytool/
planemo lint tools/mytool/
planemo test --biocontainers tools/mytool/
flake8 tools/mytool/          # only if the wrapper ships Python
planemo serve tools/mytool/   # manual check of the tool form
```

### Step 8: Run the Pre-PR Passes

Run **10a, 10b, 10c, 10d** in that order. For each item in 10a-10c record pass, fail, or N/A with a one-line reason; for 10d list only the rows that applied and what you did. Fix every fail and re-run the affected item. Do not hand back with an unexplained fail.

Put the result in the hand-back or the PR description, so a reviewer (human or agent) can audit it:

```
10a Security
 1 pass -- 6 params, all single-quoted; no nested interpreters
 2 pass -- 1 text param (sample label), validator ^[A-Za-z0-9_.-]+$, hostile-input test added
 3 pass -- element_identifier used for link names via re.sub + index suffix
 ...
 9 N/A  -- no network access
10b Readiness            (18 items, same form)
10c Remote execution     (6 items, same form)
10d Review flags         rows that applied, and what was done
Deviations (target repo convention): RST help and inline tokens, matching tools/<sibling>
Not verified: planemo serve (no browser in this session)
```

---

## 13. Quick Reference — Galaxy XML Element Ordering

The full order, including rarely used elements (`parallelism`, `code`, `stdio`, `request_param_translation`), is in [Coding Style][std-style]. The elements wrappers normally contain, in order:

```
tool
├── description
├── macros
├── xrefs
├── requirements
├── required_files (when the wrapper ships files)
├── version_command
├── command
├── environment_variables (optional)
├── configfiles (optional)
├── inputs
├── outputs
├── tests
├── help
├── citations
└── creator (optional)
```

`edam_topics` and `edam_operations` sit between `macros` and `xrefs` in the standards' order, but annotate through bio.tools instead (Section 2).

---

## Reference: Useful Planemo Commands

```bash
# Format XML files (consistent indentation, matches Galaxy Language Server)
planemo format tools/mytool/

# Format with preview (show diff without writing)
planemo format --dry-run tools/mytool/

# Lint a tool directory
planemo lint tools/mytool/

# Lint the .shed.yml (lint only; never upload by hand)
planemo shed_lint tools/mytool/

# Test with biocontainers (standard for CLI tools)
planemo test --biocontainers tools/mytool/

# Test a single tool
planemo test --biocontainers tools/mytool/mytool_align.xml

# Test without containers (API tools with fixtures)
planemo test tools/mytool/

# Serve locally for manual testing
planemo serve tools/mytool/
```

<!-- IUC standards link targets. Deep-link anchors are unverified; see the note under Authority. -->
[std]: https://galaxy-iuc-standards.readthedocs.io/en/latest/
[std-pr95]: https://github.com/galaxy-iuc/standards/pull/95
[std-tool-xml]: https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/tool_xml.html
[std-security]: https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/security.html
[std-checklist]: https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/integration_checklist.html
[std-pulsar]: https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/pulsar.html
[std-package]: https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/package_xml.html
[std-package-mulled]: https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/package_xml.html#multi-requirement-containers
[std-licensing]: https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/licensing.html
[std-shed]: https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/shed_yml.html
[std-shed-advanced]: https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/shed_yml.html#advanced-parameters
[std-repos]: https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/repositories.html
[std-dm]: https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/data_managers.html
[std-management]: https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/management.html
[std-versions]: https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/tool_xml.html#tool-versions
[std-ids]: https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/tool_xml.html#tool-ids
[std-profile]: https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/tool_xml.html#tool-profile
[std-biotools]: https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/tool_xml.html#tool-cross-references-bio-tools
[std-required-files]: https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/tool_xml.html#required-files
[std-creator]: https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/tool_xml.html#creator
[std-argument]: https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/tool_xml.html#parameter-name-argument-and-help
[std-tests]: https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/tool_xml.html#tests
[std-data-params]: https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/tool_xml.html#data-parameters
[std-booleans]: https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/tool_xml.html#booleans
[std-dynamic]: https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/tool_xml.html#dynamic-options
[std-command]: https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/tool_xml.html#command-tag
[std-exit]: https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/tool_xml.html#exit-code-detection
[std-help]: https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/tool_xml.html#help-tag
[std-images]: https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/tool_xml.html#including-images
[std-indices]: https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/tool_xml.html#generating-indices
[std-subcommands]: https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/tool_xml.html#wrapping-tools-with-subcommands
[std-style]: https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/tool_xml.html#coding-style
[sec-1]: https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/security.html#security-shell-quoting
[sec-2]: https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/security.html#security-free-text
[sec-3]: https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/security.html#security-identifiers
[sec-4]: https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/security.html#security-sanitizer-semantics
[sec-5]: https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/security.html#security-cheetah-side-effects
[sec-6]: https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/security.html#security-data-not-code
[sec-7]: https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/security.html#security-unsafe-formats
[sec-8]: https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/security.html#security-credentials
[sec-9]: https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/security.html#security-downloads
[sec-10]: https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/security.html#security-active-content
[pulsar-1]: https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/pulsar.html#pulsar-required-files
[pulsar-2]: https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/pulsar.html#pulsar-no-fs-io
[pulsar-3]: https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/pulsar.html#pulsar-data-tables
[pulsar-4]: https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/pulsar.html#pulsar-no-abs-paths
[pulsar-5]: https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/pulsar.html#pulsar-outputs-in-jobdir
[pulsar-6]: https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/pulsar.html#pulsar-discovered-outputs
[lic-osi]: https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/licensing.html#licensing-osi
[lic-redistribution]: https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/licensing.html#licensing-redistribution
