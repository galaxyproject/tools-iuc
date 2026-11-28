# Tool Help Diagrams (SVG)

Hand-coded SVG schematics that show input → operation → output for Galaxy tools. Embedded in tool XML `<help>` sections via reStructuredText `.. image::` directives.

**Embedding a diagram means RST help.** Galaxy rewrites `.. image::` paths to the tool's `static/images/` directory only for reStructuredText help (still true in 26.2). In `<help format="markdown">` the image does not resolve. So a tool that ships a help image keeps RST help for the whole section, even though Markdown is otherwise preferred for new tools (`../SKILL.md` Section 2). Do not convert an RST help section that contains images to Markdown.

**Standard:** [tool_xml: Including Images](https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/tool_xml.html#including-images).

## Where Diagrams Live

**Built-in tools** (tools bundled with the Galaxy codebase):

```
galaxy/static/images/tools/collection_ops/   # collection operation diagrams
```

**External tools** (the vast majority — tools-iuc, tools-devteam, etc.):

```
tools/<tool_name>/static/images/             # per-tool static directory
```

For example:

```
tools/phyloseq/static/images/standard_plot.png
tools/bandage/static/images/bandage_graph.png
tool_collections/bamtools/bamtools_filter/static/images/complex-filters.png
```

## How to Embed in Tool XML

Inside an RST `<help>` section (no `format` attribute, or `format="restructuredtext"`; not `format="markdown"`):

**External tools** -- reference the bare filename, as the IUC standards do. Galaxy resolves it to the tool's `static/images/` directory:

```rst
.. image:: my_diagram.svg
  :alt: Description of the diagram
  :width: 500
```

The older form `.. image:: $PATH_TO_IMAGES/my_diagram.svg` still works: Galaxy strips the `$PATH_TO_IMAGES` prefix before resolving the path, and its own comment calls the prefix no longer required. Most existing tools-iuc help sections still use it. When editing an existing tool, keep whichever form it has; for a new tool use the bare filename.

> **Verify:** the path rewriting was read from a Galaxy 26.1 development checkout (`set_image_paths`), where it runs for tools installed from a Tool Shed. A local `planemo serve` may therefore show a broken image even when the directive is right. Check a sibling tool with images the same way before concluding the path is wrong.

**Built-in Galaxy tools** — use `${static_path}` (resolves to Galaxy's central static dir):

```rst
.. image:: ${static_path}/images/tools/collection_ops/flatten.svg
  :alt: Flatten a nested collection into a simple list
  :width: 500
```

General rules:
- `:width:` — typically 500 for simple diagrams, 620–800 for complex ones
- `:alt:` — always include for accessibility
- Place after the text description, usually at the end of the help section
- Separate from preceding text with a `--------` RST horizontal rule

## SVG Structure

Every diagram follows this template:

```xml
<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 W H" width="W" height="H">
  <defs>
    <!-- Drop shadow for container boxes -->
    <filter id="shadow" x="-4%" y="-4%" width="110%" height="110%">
      <feGaussianBlur in="SourceAlpha" stdDeviation="2" result="blur"/>
      <feComponentTransfer in="blur">
        <feFuncA type="linear" slope="0.3"/>
      </feComponentTransfer>
      <feMerge>
        <feMergeNode/>
        <feMergeNode in="SourceGraphic"/>
      </feMerge>
    </filter>
    <!-- Arrow marker for connector lines -->
    <marker id="arrowhead" markerWidth="10" markerHeight="7"
            refX="9" refY="3.5" orient="auto">
      <polygon points="0 0, 10 3.5, 0 7" fill="#555"/>
    </marker>
  </defs>

  <!-- LEFT: input -->
  ...
  <!-- CENTER: arrow + operation label -->
  ...
  <!-- RIGHT: output -->
  ...
</svg>
```

### Typical viewBox Sizes

| Complexity | viewBox | When |
|---|---|---|
| Simple (1→1) | `620 210` | filter, sort, relabel |
| Medium | `620 260` or `680 270` | flatten, zip, merge |
| Complex / multi-case | `800 640` | build_list (3 cases) |

## Color Palette

### Container Boxes (Collections)

| Element | Fill | Stroke | Notes |
|---|---|---|---|
| Collection (outer) | `#d9ead3` (light green) | `#888` | `filter="url(#shadow)"`, `rx="8"` |
| Nested container | `#e8f0e3` (lighter green) | `#aaa` | `rx="6"` |
| Forward reads | `#fef3e2` (light orange) | `#d4760a` | Paired: forward/R1 |
| Reverse reads | `#f0e8f5` (light purple) | `#7b4fa0` | Paired: reverse/R2 |
| Warning/special state | `#fff3cd` (light yellow) | `#e6a817` | Empty, error, null |

### Dataset Boxes (Inside Collections)

| Element | Fill | Stroke | Notes |
|---|---|---|---|
| Normal dataset | `#fff` (white) | `#aaa`, `stroke-width="0.7"` | `rx="4"` |
| Forward dataset | `#fef3e2` | `#d4760a`, `stroke-width="0.5"` | Inside paired container |
| Reverse dataset | `#f0e8f5` | `#7b4fa0`, `stroke-width="0.5"` | Inside paired container |
| Highlighted (empty/error) | `#fff3cd` | `#e6a817`, `stroke-width="1.2"` | Draws attention |

### Text

| Role | Font | Size | Weight | Fill |
|---|---|---|---|---|
| Collection title | sans-serif | 12–13 | bold | `#333` |
| Dataset label | sans-serif | 10–12 | normal | `#333` or color-matched |
| Nested pair label | sans-serif | 10 | bold | `#555` |
| Operation label (on arrow) | sans-serif | 12 | bold | `#555` |
| Arrow parameter text | sans-serif | 10 | normal | `#888` |
| Annotation (bottom) | sans-serif | 10 | italic | `#666` |
| Size/metadata | sans-serif | 9 | normal | `#999` or `#b45309` |

### Arrows and Lines

| Element | Stroke | Width | Notes |
|---|---|---|---|
| Main arrow | `#555` | 2 | `marker-end="url(#arrowhead)"` |
| Dashed separator | `#ccc` | 1 | `stroke-dasharray="6,4"` — between cases |

## Layout Conventions

### Standard Layout (left → right)

```
[Input collection]  ——arrow——>  [Output collection]
     x≈10–50          x≈240–400       x≈380–450
```

- Input on left, output on right
- Arrow centered vertically between input and output
- Operation name above the arrow, optional parameter text below
- Bottom annotation in italic explains the transformation rule

### Spacing

- Outer collection padding: 12–15px
- Between dataset boxes: 6–8px vertical gap
- Dataset box height: 26–28px (normal), 14–20px (nested/compact)
- Collection border radius: `rx="8"` (outer), `rx="6"` (nested), `rx="4"` (dataset)

### Multi-Case Diagrams

When a tool has multiple usage patterns (e.g., build_list):

1. Stack cases vertically, labeled A, B, C
2. Separate with dashed lines: `<line ... stroke="#ccc" stroke-dasharray="6,4"/>`
3. Each case is a self-contained input→arrow→output row

## Step-by-Step: Creating a New Diagram

1. **Decide what to show** — pick the simplest example that demonstrates the tool's core transformation. Use 2–3 elements, not more.

2. **Sketch the layout** — input structure on left, output on right. Determine viewBox size from the complexity table above.

3. **Start from the template** — copy the `<defs>` block (shadow + arrowhead) from above. These are identical in every diagram.

4. **Draw input collection** — outer `<rect>` with `fill="#d9ead3"`, then dataset `<rect>`s inside with `fill="#fff"`. Add `<text>` labels.

5. **Draw the arrow** — horizontal `<line>` with `marker-end="url(#arrowhead)"`. Add operation name above.

6. **Draw output collection** — same style as input, showing the result.

7. **Add annotation** — italic text at the bottom explaining the transformation rule (e.g., "hierarchy removed, identifiers joined").

8. **Save as SVG** in `<tool_dir>/static/images/` (external tools) or `static/images/tools/<category>/` (built-in tools).

9. **Add `.. image::` directive** to the tool's `<help>` section, which must be RST (see the top of this file).

## Complete Example: Filter Empty

A simple 1-input → 1-output diagram highlighting the removed element:

```xml
<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 620 210" width="620" height="210">
  <defs>
    <filter id="shadow" x="-4%" y="-4%" width="110%" height="110%">
      <feGaussianBlur in="SourceAlpha" stdDeviation="2" result="blur"/>
      <feComponentTransfer in="blur">
        <feFuncA type="linear" slope="0.3"/>
      </feComponentTransfer>
      <feMerge>
        <feMergeNode/>
        <feMergeNode in="SourceGraphic"/>
      </feMerge>
    </filter>
    <marker id="arrowhead" markerWidth="10" markerHeight="7"
            refX="9" refY="3.5" orient="auto">
      <polygon points="0 0, 10 3.5, 0 7" fill="#555"/>
    </marker>
  </defs>

  <!-- Input collection -->
  <rect x="10" y="10" width="210" height="185" rx="8" ry="8"
        fill="#d9ead3" stroke="#888" stroke-width="1" filter="url(#shadow)"/>
  <text x="115" y="30" text-anchor="middle" font-family="sans-serif"
        font-size="13" font-weight="bold" fill="#333">List</text>

  <rect x="25" y="42" width="180" height="28" rx="4" ry="4"
        fill="#fff" stroke="#aaa" stroke-width="0.7"/>
  <text x="115" y="61" text-anchor="middle" font-family="sans-serif"
        font-size="12" fill="#333">sample1</text>
  <text x="195" y="61" text-anchor="end" font-family="sans-serif"
        font-size="9" fill="#999">1.2 MB</text>

  <!-- Highlighted empty element -->
  <rect x="25" y="78" width="180" height="28" rx="4" ry="4"
        fill="#fff3cd" stroke="#e6a817" stroke-width="1.2"/>
  <text x="115" y="97" text-anchor="middle" font-family="sans-serif"
        font-size="12" fill="#b45309">sample2</text>
  <text x="195" y="97" text-anchor="end" font-family="sans-serif"
        font-size="9" fill="#b45309">0 bytes</text>

  <rect x="25" y="114" width="180" height="28" rx="4" ry="4"
        fill="#fff" stroke="#aaa" stroke-width="0.7"/>
  <text x="115" y="133" text-anchor="middle" font-family="sans-serif"
        font-size="12" fill="#333">sample3</text>

  <text x="115" y="166" text-anchor="middle" font-family="sans-serif"
        font-size="10" font-style="italic" fill="#b45309">sample2 is empty</text>

  <!-- Arrow -->
  <line x1="240" y1="92" x2="360" y2="92" stroke="#555" stroke-width="2"
        marker-end="url(#arrowhead)"/>
  <text x="300" y="82" text-anchor="middle" font-family="sans-serif"
        font-size="12" font-weight="bold" fill="#555">Filter</text>

  <!-- Output collection -->
  <rect x="380" y="10" width="220" height="155" rx="8" ry="8"
        fill="#d9ead3" stroke="#888" stroke-width="1" filter="url(#shadow)"/>
  <text x="490" y="30" text-anchor="middle" font-family="sans-serif"
        font-size="13" font-weight="bold" fill="#333">List (filtered)</text>

  <rect x="395" y="42" width="190" height="28" rx="4" ry="4"
        fill="#fff" stroke="#aaa" stroke-width="0.7"/>
  <text x="490" y="61" text-anchor="middle" font-family="sans-serif"
        font-size="12" fill="#333">sample1</text>

  <rect x="395" y="78" width="190" height="28" rx="4" ry="4"
        fill="#fff" stroke="#aaa" stroke-width="0.7"/>
  <text x="490" y="97" text-anchor="middle" font-family="sans-serif"
        font-size="12" fill="#333">sample3</text>

  <text x="490" y="140" text-anchor="middle" font-family="sans-serif"
        font-size="10" font-style="italic" fill="#666">empty element removed</text>
</svg>
```

## Checklist

- [ ] SVG uses the standard `<defs>` block (shadow filter + arrowhead marker)
- [ ] Colors match the palette above
- [ ] Layout is left-to-right: input → arrow → output
- [ ] Uses 2–3 example elements (not more)
- [ ] Bottom annotation explains the transformation
- [ ] Saved to `static/images/` (external tools) or `static/images/tools/<category>/` (built-in)
- [ ] Tool `<help>` is RST, not `format="markdown"`
- [ ] `.. image::` directive added to tool `<help>` using the bare filename (external) or `${static_path}` (built-in), with `:alt:` and `:width:`
- [ ] Renders correctly in a browser before committing
