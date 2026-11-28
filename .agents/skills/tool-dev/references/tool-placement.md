# Galaxy Tool Placement Decision Guide

Where to create your Galaxy tool wrapper.

---

## Decision Tree

```
Where should I create this tool?

1. Check if tool already exists
   ├─ tools-iuc? → Use it!
   ├─ Other known repo? → Use it!
   └─ ToolShed? → Evaluate quality, use if good

2. Tool doesn't exist - where to create?
   ├─ Community-useful? → tools-iuc (preferred)
   ├─ Domain-specific? → Appropriate repo (genouest, bgruening, etc.)
   └─ Project-specific? → Custom/local
```

---

## Option 1: tools-iuc (Preferred)

**Repository**: https://github.com/galaxyproject/tools-iuc

**Create here when**:
- The software is under an OSI-approved license and nothing it needs (dependencies, models, reference data) carries a use restriction. This is a gate, not a preference: check it before writing anything. See [Licensing and Use Restrictions](https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/licensing.html)
- Tool is widely used in the community
- Tool is well-maintained upstream
- Tool would benefit multiple users
- Tool is not project-specific

**Advantages**:
- High quality, well-tested
- Actively maintained
- Available on usegalaxy.* servers
- Automatic CI/CD
- Community support

**Process**:
1. Check existing PRs — someone may already be working on it
2. Open an issue to discuss the tool before implementing
3. Fork the repository
4. Create tool following the [IUC standards](https://galaxy-iuc-standards.readthedocs.io/en/latest/) (how-to in `../SKILL.md`)
5. Run the four pre-PR passes (`../SKILL.md` Section 10) and put the result in the PR description
6. Submit PR
7. Address review feedback

**IUC requirements** (summary only; the standards are authoritative, and the target repo's own convention wins where it differs):
- Follow directory structure: `tools/toolname/`
- Version from `@TOOL_VERSION@` / `@VERSION_SUFFIX@` tokens. A shared `macros.xml` only when two or more tools share something; a single tool can keep tokens inline
- Include comprehensive tests with `expect_num_outputs` and `assert_contents`
- Test data: aim for under 1 MB. Unavoidable large inputs go to the Galaxy Tool Developers Zenodo community and are referenced with `location=`
- Write a clear help section. Markdown (`format="markdown"`, profile >= 24.2) is preferred for new tools; use RST when the tool ships help images or the repo's sibling tools use RST. Existing RST help need not be converted
- Add citations (`type="doi"`)
- Use `argument=` on CLI-mapped params
- Bioconda/conda-forge dependencies, pinned by version

---

## Option 2: Domain-Specific Repositories

| Repository | Focus | When to Use |
|------------|-------|-------------|
| **tools-iuc** | General bioinformatics | Default choice for community tools |
| **genouest/galaxy-tools** | Genomics, annotation | Genomics-specific tools (BRAKER3, etc.) |
| **bgruening/galaxytools** | Cheminformatics, misc | Chemistry, specialized tools |
| **ARTbio/tools-artbio** | RNA-seq, small RNA | RNA analysis tools |

---

## Option 3: Custom/Local Tools

**Create custom when**:
- Tool is project-specific
- Tool is not useful to broader community
- Quick prototyping needed
- Workflow-embedded tool is sufficient

### Standalone Custom Tool

```
my-project/
└── tools/
    └── mytool/
        ├── mytool.xml
        ├── macros.xml
        └── test-data/
```

Advantages: full control, fast iteration, no review process.
Disadvantages: not available to community, manual maintenance, no automatic CI/CD.

### Workflow-Embedded Tool

Tool definition embedded in workflow `.ga` file. Only for very simple tools (< 50 lines) used in a single workflow. Limited to expression tools.

---

## Evaluation Criteria

### Choose tools-iuc if:

- Tool is used by multiple research groups
- Tool is actively maintained upstream
- Tool has stable API/interface
- You can commit to maintaining it
- Tool follows bioinformatics best practices

### Choose custom if:

- Tool is project-specific
- Tool is experimental/prototype
- Tool is simple wrapper script
- No community benefit expected

---

## Quality Standards

Regardless of where you create the tool:

**Required**: valid XML structure, version information, at least one test case, help section, proper input/output definitions.

**Recommended**: multiple test cases, comprehensive help, citations, version tokens, proper error handling.

**Best practice**: test data under 1 MB where possible, flexible test assertions, clear parameter descriptions, sensible defaults.

The security and remote-execution passes in `../SKILL.md` Section 10 apply wherever the tool will run on a shared server, not only to tools-iuc submissions.

---

## Migration Path

It's OK to start with a custom tool for development, then submit to tools-iuc once stable. This approach allows rapid iteration, proves tool usefulness, and results in better quality submissions.

---

## Related

- **Full tool development guide**: `../SKILL.md`
- **Testing guide**: `testing.md`
