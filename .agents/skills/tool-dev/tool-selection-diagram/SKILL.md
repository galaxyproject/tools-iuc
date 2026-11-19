---
name: tool-selection-diagram
description: Generate a tool selection flowchart diagram (PNG) for a Galaxy tool suite from a JSON definition. Use when creating visual guides that help users pick the right tool in a multi-tool suite.
---

# Tool Selection Diagram Generator

Generate "which tool should I use?" flowchart diagrams for multi-tool Galaxy suites (e.g., MMseqs2, Samtools). Produces PNG images using Galaxy's color palette and Atkinson Hyperlegible font.

## When to Use

- Creating a new visual guide for a multi-tool suite
- Updating an existing tool selection diagram after adding/removing tools
- A tool suite has 3+ tools and users frequently ask "which one should I use?"

## Quick Reference

| Item | Value |
|------|-------|
| **Input** | JSON definition file |
| **Output** | PNG flowchart image |
| **Script** | `tool-selection-diagram/scripts/render_tool_diagram.py` |
| **Dependencies** | `Pillow` (Python) |
| **Font** | Atkinson Hyperlegible (falls back to DejaVu Sans) |
| **Colors** | Galaxy "Paired" colormap from `gxy-colors.svg` |

---

## JSON Schema

The JSON definition describes a tree: start question -> goals -> (optional criteria) -> tools.

```json
{
  "title": "string (required) — diagram caption, shown at bottom",
  "start_question": "string (required) — root node text, e.g. 'What do you need to do?'",
  "goals": [
    {
      "label": "string (required) — analysis goal name",
      "criteria": [
        {
          "label": "string (required) — decision criterion",
          "tools": [
            {
              "name": "string (required) — Galaxy tool name as shown in UI",
              "description": "string (optional) — short description below node"
            }
          ]
        }
      ],
      "tools": [
        {
          "name": "string (required)",
          "description": "string (optional)"
        }
      ]
    }
  ]
}
```

### Field Details

| Field | Type | Required | Notes |
|-------|------|----------|-------|
| `title` | string | Yes | Diagram caption (e.g., "MMseqs2 Tool Selection Guide") |
| `start_question` | string | Yes | Root node — the question users start with |
| `goals[].label` | string | Yes | Analysis goal (tier 2 node) |
| `goals[].criteria` | array | No | Decision criteria (tier 3 nodes). Mutually exclusive with `goals[].tools` at the same level. |
| `goals[].criteria[].label` | string | Yes | Criterion text |
| `goals[].criteria[].tools` | array | Yes | Tools under this criterion |
| `goals[].tools` | array | No | Tools directly under a goal (3-tier path, no criteria). Mutually exclusive with `goals[].criteria`. |
| `tools[].name` | string | Yes | Tool name as it appears in Galaxy |
| `tools[].description` | string | No | Short description shown below the tool node |

Each goal uses **either** `criteria` (4-tier path) **or** `tools` (3-tier path), never both. A single diagram can mix both patterns across different goals.

---

## Workflow

### Step 1: Identify tools in the suite

List all tools in the suite. Group them by analysis goal (what the user is trying to accomplish). Each group becomes a goal node.

### Step 2: Decide tier structure per goal

For each goal, decide whether users need a decision criterion to choose between tools:

| Situation | Tier Structure | Example |
|-----------|---------------|---------|
| Only 1-2 tools, choice is obvious | 3-tier (goal -> tools) | Sort/Index -> samtools sort, samtools index |
| Users need guidance to pick | 4-tier (goal -> criterion -> tools) | Search -> "Large DB" -> linsearch, "Sensitive" -> search |

Present the proposed structure to the user for approval before drafting JSON.

### Step 3: Draft JSON definition

Use the examples in `examples/` as templates:
- `simple_example.json` — minimal 3-tier (no criteria)
- `mmseqs2_tool_selection.json` — mixed 3-tier and 4-tier
- `samtools_tool_selection.json` — large suite, mixed tiers, multi-tool criteria

### Step 4: Run the renderer

```bash
python3 tool-selection-diagram/scripts/render_tool_diagram.py \
    --input definition.json \
    --output tools/<suite>/static/images/<suite>_tool_selection.png \
    --dpi 150
```

The script requires `Pillow`. Install if needed: `pip install Pillow`.

### Step 5: Verify the PNG

Open the PNG and check:
- All tool names match their Galaxy tool names exactly
- No text is clipped or overlapping
- The tree structure matches the intended decision flow
- Descriptions are concise and readable

### Step 6: Embed in tool help

Reference the diagram in each tool's `<help>` section so users see it in Galaxy. See **Embedding in Tool Help** below.

---

## Tier Structure Guide

The renderer supports mixed-tier diagrams. Each goal independently chooses its path:

| Pattern | When to Use | Structure |
|---------|-------------|-----------|
| **3-tier** | Tools under a goal are self-explanatory or few | Start -> Goal -> Tool(s) |
| **4-tier** | Users need a decision criterion to choose | Start -> Goal -> Criterion -> Tool(s) |
| **Mixed** | Some goals need criteria, others don't | Mix of both in one diagram |

A criterion node can have multiple tools beneath it (e.g., "Quick summary" -> flagstat + idxstats).

---

## Color Palette

Colors come from Galaxy's "Paired" colormap (`gxy-colors.svg`):

| Tier | Color | Hex | Usage |
|------|-------|-----|-------|
| Start | Dark navy | `#2c3143` | Root question node |
| Goal | Blue | `#2077b3` | Analysis goal nodes |
| Criterion | Light orange | `#fdbf6f` | Decision criterion nodes |
| Tool | Green | `#74c376` | Galaxy tool nodes |

Text colors: white on dark backgrounds (Start, Goal), dark navy on light backgrounds (Criterion, Tool).

---

## Common Mistakes

| Mistake | Why It's Wrong | Fix |
|---------|----------------|-----|
| Tool names don't match Galaxy UI | Users can't find the tool | Use exact Galaxy tool names (e.g., "samtools view" not "SAMtools View") |
| Too many tools under one criterion | Diagram becomes unreadable | Split into multiple criteria or group differently |
| Missing `start_question` | Script will error | Always include a root question |
| Labels too long (>40 chars) | Text gets clipped in nodes | Shorten labels; use `description` for details |
| Using both `criteria` and `tools` on one goal | Invalid JSON structure | Each goal uses one or the other, not both |
| Missing `title` | No caption on diagram | Always set a descriptive title |
| Descriptions on non-tool nodes | Only tool nodes show descriptions | Move extra text to criterion labels |
| Not verifying PNG after generation | Overlapping text or layout issues | Always open and visually check the output |

---

## Examples

### Simple 3-tier (no criteria)

`examples/simple_example.json` — 2 goals, 3 tools, direct goal-to-tool paths:

```json
{
  "title": "My Toolkit Selection Guide",
  "start_question": "What do you want to do?",
  "goals": [
    {
      "label": "Analyze data",
      "tools": [
        {"name": "tool-analyze", "description": "Run full analysis"},
        {"name": "tool-summarize", "description": "Quick summary stats"}
      ]
    },
    {
      "label": "Visualize results",
      "tools": [
        {"name": "tool-plot", "description": "Generate charts"}
      ]
    }
  ]
}
```

### Mixed 3-tier and 4-tier

`examples/mmseqs2_tool_selection.json` — 4 goals, some with criteria (search, cluster, taxonomy) and one without (orthologs):

```json
{
  "goals": [
    {
      "label": "Search sequences",
      "criteria": [
        {"label": "Large DB / high ID", "tools": [{"name": "easy-linsearch"}]},
        {"label": "General / remote", "tools": [{"name": "easy-search"}]}
      ]
    },
    {
      "label": "Find orthologs",
      "tools": [{"name": "easy-rbh", "description": "1 : 1 reciprocal best hit"}]
    }
  ]
}
```

### Large suite with multi-tool criteria

`examples/samtools_tool_selection.json` — 5 goals, 12 tools. Some criteria point to multiple tools (e.g., "Quick summary" -> flagstat + idxstats).

---

## Embedding in Tool Help

Add the diagram to each tool's `<help>` section using a relative image path:

```xml
<help format="markdown"><![CDATA[
**Which Samtools tool should I use?**

![Samtools Tool Selection Guide](static/images/samtools_tool_selection.png)

-----

**samtools view** filters and converts SAM/BAM/CRAM files...
]]></help>
```

Place the PNG at:
```
tools/<suite>/
├── static/
│   └── images/
│       └── <suite>_tool_selection.png
├── macros.xml
├── <suite>_tool_a.xml
└── ...
```

Galaxy serves files from `static/` relative to the tool directory, so the `![](static/images/...)` path works in the rendered help.

---

## Script Reference

```
render_tool_diagram.py --input <JSON> --output <PNG> [--dpi N]
```

| Flag | Required | Default | Description |
|------|----------|---------|-------------|
| `--input` | Yes | — | Path to JSON definition file |
| `--output` | Yes | — | Output PNG path |
| `--dpi` | No | 150 | Output resolution |

The script auto-detects tier structure (3-tier, 4-tier, or mixed) from the JSON. No configuration needed beyond the JSON definition.
