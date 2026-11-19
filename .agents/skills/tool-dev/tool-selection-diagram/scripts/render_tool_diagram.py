#!/usr/bin/env python3
"""
Render a tool selection flowchart diagram from a JSON definition.

Produces a PNG image using Galaxy's color palette (from gxy-colors.svg):
  Start question -> Analysis goals -> Decision criteria -> Galaxy tools

Color palette:  gxy-colors.svg "Paired" colormap + Galaxy logo colors
Font:           Atkinson Hyperlegible ($font-family-base)

Usage:
    python3 render_tool_diagram.py --input definition.json --output diagram.png [--dpi 150]
"""

import argparse
import collections
import json
import math
import os

from PIL import Image, ImageDraw, ImageFont

# ── Galaxy palette (from gxy-colors.svg) ─────────────────────────────────
# "Paired" colormap
GXY_ORANGE_LIGHT = "#fdbf6f"   # paired light orange
GXY_BLUE_LIGHT = "#a3cbdf"     # paired light blue
GXY_GREEN_LIGHT = "#9fd99a"    # paired light green
GXY_RED_LIGHT = "#fb9998"      # paired light red
GXY_ORANGE = "#fe7f02"         # paired dark orange
GXY_BLUE = "#2077b3"           # paired dark blue
GXY_GREEN = "#74c376"          # paired dark green
GXY_RED = "#e31a1e"            # paired dark red
# "True Galaxy logo colors"
GXY_YELLOW = "#ffd602"         # logo yellow
GXY_DARK = "#2c3143"           # logo dark navy
GXY_GREY = "#cccccc"           # logo grey
GXY_LIGHT_GREY = "#f2f2f2"     # logo light grey

BG_COLOR = "#FFFFFF"

# ── Fonts ─────────────────────────────────────────────────────────────────
# Galaxy $font-family-base: "Atkinson Hyperlegible", -apple-system, …
_ATKINSON_DIR = os.path.expanduser("~/.local/share/fonts")
_FONT_SEARCH_DIRS = [
    _ATKINSON_DIR,
    "/usr/share/fonts/truetype/dejavu",           # Debian/Ubuntu
    "/usr/share/fonts/dejavu-sans-fonts",          # Fedora/RHEL
    "/Library/Fonts",                              # macOS system
    "/System/Library/Fonts/Supplemental",          # macOS supplemental
    os.path.expanduser("~/Library/Fonts"),         # macOS user
]


def _pick_font(bold=False):
    # Preferred: Atkinson Hyperlegible
    suffix = "Bold" if bold else "Regular"
    atkinson = f"{_ATKINSON_DIR}/AtkinsonHyperlegible-{suffix}.ttf"
    if os.path.isfile(atkinson):
        return atkinson

    # Fallback: search for common sans-serif fonts across platforms
    bold_suffix = "-Bold" if bold else ""
    candidates = [
        f"DejaVuSans{bold_suffix}.ttf",
        f"Arial{' Bold' if bold else ''}.ttf",
        f"Helvetica{' Bold' if bold else ''}.ttc",
    ]
    for font_dir in _FONT_SEARCH_DIRS:
        for name in candidates:
            path = os.path.join(font_dir, name)
            if os.path.isfile(path):
                return path
    return None


FONT_BOLD = _pick_font(bold=True)
FONT_REGULAR = _pick_font(bold=False)
FONT_ITALIC = (
    f"{_ATKINSON_DIR}/AtkinsonHyperlegible-Italic.ttf"
    if os.path.isfile(f"{_ATKINSON_DIR}/AtkinsonHyperlegible-Italic.ttf")
    else FONT_REGULAR
)

# ── Tier styles ───────────────────────────────────────────────────────────
TIER_STYLES = {
    "start":     {"bg": GXY_DARK,         "fg": "#f8f9fa",  "font": FONT_BOLD, "size": 22},
    "goal":      {"bg": GXY_BLUE,         "fg": "#f8f9fa",  "font": FONT_BOLD, "size": 20},
    "criterion": {"bg": GXY_ORANGE_LIGHT, "fg": GXY_DARK,   "font": FONT_BOLD, "size": 18},
    "tool":      {"bg": GXY_GREEN,        "fg": GXY_DARK,   "font": FONT_BOLD, "size": 20},
}

# ── Layout constants ──────────────────────────────────────────────────────
DESC_FONT_SIZE = 14
LEGEND_FONT_SIZE = 16
TITLE_FONT_SIZE = 18

NODE_HEIGHT = 60
NODE_CORNER_RADIUS = 11
NODE_H_PAD = 30
NODE_MIN_WIDTH = 180
NODE_MAX_WIDTH = 400

ROW_SPACING = 190
MARGIN_X = 80
MARGIN_TOP = 60

TOOL_H_GAP = 40
CLUSTER_GAP = 80

DESC_OFFSET_Y = 12
ARROW_COLOR = GXY_DARK
ARROW_WIDTH = 2
ARROWHEAD_LENGTH = 12
ARROWHEAD_HALF_WIDTH = 5

LEGEND_SWATCH_SIZE = 28
LEGEND_SWATCH_RADIUS = 5
LEGEND_LABEL_GAP = 10
LEGEND_ITEM_GAP = 50
LEGEND_MARGIN_TOP = 50

TITLE_MARGIN_TOP = 20
BOTTOM_MARGIN = 50


# ── Helpers ───────────────────────────────────────────────────────────────

def hex_to_rgb(hex_color):
    h = hex_color.lstrip("#")
    return tuple(int(h[i:i + 2], 16) for i in (0, 2, 4))


def load_font(path, size):
    if path is None:
        return ImageFont.load_default()
    try:
        return ImageFont.truetype(path, size)
    except (OSError, TypeError, AttributeError):
        return ImageFont.load_default()


def text_width(draw, text, font):
    bbox = draw.textbbox((0, 0), text, font=font)
    return bbox[2] - bbox[0]


def node_width_for_text(draw, text, font, max_width=NODE_MAX_WIDTH):
    tw = text_width(draw, text, font)
    w = tw + NODE_H_PAD * 2
    return max(NODE_MIN_WIDTH, min(max_width, w))


def draw_arrow(draw, x1, y1, x2, y2):
    """Line from (x1,y1) to (x2,y2) with a proper angled arrowhead at the tip."""
    color = hex_to_rgb(ARROW_COLOR)
    draw.line([(x1, y1), (x2, y2)], fill=color, width=ARROW_WIDTH)

    # Angle of the line
    angle = math.atan2(y2 - y1, x2 - x1)

    # Arrowhead: a triangle whose tip is at (x2, y2), pointing along the line
    # Two base points are behind the tip, spread perpendicular to the line
    back_x = x2 - ARROWHEAD_LENGTH * math.cos(angle)
    back_y = y2 - ARROWHEAD_LENGTH * math.sin(angle)
    perp_angle = angle + math.pi / 2
    dx = ARROWHEAD_HALF_WIDTH * math.cos(perp_angle)
    dy = ARROWHEAD_HALF_WIDTH * math.sin(perp_angle)

    draw.polygon(
        [(x2, y2), (back_x + dx, back_y + dy), (back_x - dx, back_y - dy)],
        fill=color,
    )


# ── Node ──────────────────────────────────────────────────────────────────

class Node:
    __slots__ = ("label", "tier", "description", "x", "y", "w", "h",
                 "children", "parent", "span")

    def __init__(self, label, tier, description=None):
        self.label = label
        self.tier = tier
        self.description = description
        self.x = 0.0
        self.y = 0.0
        self.w = 0
        self.h = NODE_HEIGHT
        self.children = []
        self.parent = None
        self.span = 0.0

    def top_center(self):
        return self.x, self.y - self.h / 2

    def bottom_center(self):
        return self.x, self.y + self.h / 2

    def bbox(self):
        hw = self.w / 2
        hh = self.h / 2
        return (self.x - hw, self.y - hh, self.x + hw, self.y + hh)

    def leaf_count(self):
        if self.tier == "tool":
            return 1
        if not self.children:
            return 1
        return sum(c.leaf_count() for c in self.children)


# ── Tree construction ─────────────────────────────────────────────────────

def build_tree(definition):
    root = Node(definition["start_question"], "start")
    for goal_def in definition["goals"]:
        goal = Node(goal_def["label"], "goal")
        goal.parent = root
        root.children.append(goal)
        if "criteria" in goal_def:
            for crit_def in goal_def["criteria"]:
                crit = Node(crit_def["label"], "criterion")
                crit.parent = goal
                goal.children.append(crit)
                for tool_def in crit_def["tools"]:
                    tool = Node(tool_def["name"], "tool", tool_def.get("description"))
                    tool.parent = crit
                    crit.children.append(tool)
        elif "tools" in goal_def:
            for tool_def in goal_def["tools"]:
                tool = Node(tool_def["name"], "tool", tool_def.get("description"))
                tool.parent = goal
                goal.children.append(tool)
    return root


# ── Layout ────────────────────────────────────────────────────────────────

def compute_widths(draw, root):
    stack = [root]
    while stack:
        node = stack.pop()
        style = TIER_STYLES[node.tier]
        font = load_font(style["font"], style["size"])
        max_w = 9999 if node.tier == "start" else NODE_MAX_WIDTH
        node.w = node_width_for_text(draw, node.label, font, max_width=max_w)
        stack.extend(node.children)


def compute_spans(node):
    if node.tier == "tool":
        node.span = node.w
        return
    for child in node.children:
        compute_spans(child)
    if not node.children:
        node.span = node.w
        return
    children_span = sum(c.span for c in node.children)
    gaps = TOOL_H_GAP * (len(node.children) - 1) if node.tier == "criterion" else \
           CLUSTER_GAP * (len(node.children) - 1)
    node.span = max(node.w, children_span + gaps)


def position_subtree(node, center_x):
    node.x = center_x
    if not node.children:
        return
    gap = TOOL_H_GAP if node.tier == "criterion" else CLUSTER_GAP
    total_child_span = sum(c.span for c in node.children)
    total_gaps = gap * (len(node.children) - 1)
    used = total_child_span + total_gaps
    left = center_x - used / 2
    for child in node.children:
        child_center = left + child.span / 2
        position_subtree(child, child_center)
        left += child.span + gap


def has_criteria(root):
    stack = [root]
    while stack:
        node = stack.pop()
        if node.tier == "criterion":
            return True
        stack.extend(node.children)
    return False


def assign_y(root):
    top = MARGIN_TOP + NODE_HEIGHT / 2
    use_criteria_row = has_criteria(root)
    if use_criteria_row:
        y_map = {
            "start":     top,
            "goal":      top + ROW_SPACING,
            "criterion": top + ROW_SPACING * 2,
            "tool":      top + ROW_SPACING * 3,
        }
    else:
        y_map = {
            "start":     top,
            "goal":      top + ROW_SPACING,
            "criterion": top + ROW_SPACING * 2,
            "tool":      top + ROW_SPACING * 2,
        }
    stack = [root]
    while stack:
        node = stack.pop()
        node.y = y_map[node.tier]
        stack.extend(node.children)
    return y_map


def layout(draw, root):
    compute_widths(draw, root)
    compute_spans(root)
    canvas_width = root.span + MARGIN_X * 2
    canvas_width = max(canvas_width, 600)
    position_subtree(root, canvas_width / 2)
    y_map = assign_y(root)
    tool_bottom = y_map["tool"] + NODE_HEIGHT / 2
    desc_space = DESC_OFFSET_Y + 20
    legend_space = LEGEND_MARGIN_TOP + LEGEND_SWATCH_SIZE
    title_space = TITLE_MARGIN_TOP + 25
    canvas_height = int(tool_bottom + desc_space + legend_space + title_space + BOTTOM_MARGIN)
    return int(canvas_width), canvas_height, y_map


# ── Rendering ─────────────────────────────────────────────────────────────

def collect_nodes(root):
    nodes = []
    queue = collections.deque([root])
    while queue:
        node = queue.popleft()
        nodes.append(node)
        queue.extend(node.children)
    return nodes


def render(definition, output_path, dpi=150):
    root = build_tree(definition)

    tmp = Image.new("RGB", (1, 1))
    tmp_draw = ImageDraw.Draw(tmp)
    canvas_w, canvas_h, y_map = layout(tmp_draw, root)

    img = Image.new("RGB", (canvas_w, canvas_h), BG_COLOR)
    draw = ImageDraw.Draw(img)
    all_nodes = collect_nodes(root)

    # Draw arrows (behind nodes)
    for node in all_nodes:
        for child in node.children:
            bx, by = node.bottom_center()
            tx, ty = child.top_center()
            draw_arrow(draw, bx, by, tx, ty)

    # Draw nodes
    desc_font = load_font(FONT_ITALIC, DESC_FONT_SIZE)

    for node in all_nodes:
        style = TIER_STYLES[node.tier]
        font = load_font(style["font"], style["size"])
        x0, y0, x1, y1 = node.bbox()

        draw.rounded_rectangle(
            (x0, y0, x1, y1),
            radius=NODE_CORNER_RADIUS,
            fill=hex_to_rgb(style["bg"]),
        )

        # Text centered with anchor="mm" (middle-middle) for proper vertical centering
        draw.text(
            (node.x, node.y), node.label,
            fill=hex_to_rgb(style["fg"]), font=font, anchor="mm",
        )

        # Description below tool nodes
        if node.tier == "tool" and node.description:
            dy = y1 + DESC_OFFSET_Y
            draw.text(
                (node.x, dy), node.description,
                fill=hex_to_rgb(GXY_DARK), font=desc_font, anchor="mt",
            )

    # Draw legend
    legend_font = load_font(FONT_REGULAR, LEGEND_FONT_SIZE)
    legend_items = [
        ("Start", TIER_STYLES["start"]["bg"]),
        ("Analysis goal", TIER_STYLES["goal"]["bg"]),
        ("Decision criterion", TIER_STYLES["criterion"]["bg"]),
        ("Galaxy tool", TIER_STYLES["tool"]["bg"]),
    ]

    legend_total_w = 0
    for label, _ in legend_items:
        legend_total_w += LEGEND_SWATCH_SIZE + LEGEND_LABEL_GAP + text_width(draw, label, legend_font)
    legend_total_w += LEGEND_ITEM_GAP * (len(legend_items) - 1)

    tool_bottom = y_map["tool"] + NODE_HEIGHT / 2
    legend_y = tool_bottom + DESC_OFFSET_Y + 20 + LEGEND_MARGIN_TOP
    lx = (canvas_w - legend_total_w) / 2

    for label, color in legend_items:
        sy = legend_y
        draw.rounded_rectangle(
            (lx, sy, lx + LEGEND_SWATCH_SIZE, sy + LEGEND_SWATCH_SIZE),
            radius=LEGEND_SWATCH_RADIUS,
            fill=hex_to_rgb(color),
        )
        lw = text_width(draw, label, legend_font)
        draw.text(
            (lx + LEGEND_SWATCH_SIZE + LEGEND_LABEL_GAP,
             sy + LEGEND_SWATCH_SIZE / 2),
            label,
            fill=hex_to_rgb(GXY_DARK), font=legend_font, anchor="lm",
        )
        lx += LEGEND_SWATCH_SIZE + LEGEND_LABEL_GAP + lw + LEGEND_ITEM_GAP

    # Draw title
    title = definition.get("title", "")
    if title:
        title_font = load_font(FONT_REGULAR, TITLE_FONT_SIZE)
        title_y = legend_y + LEGEND_SWATCH_SIZE + TITLE_MARGIN_TOP
        draw.text(
            (canvas_w / 2, title_y), title,
            fill=hex_to_rgb(GXY_DARK), font=title_font, anchor="mt",
        )

    img.save(output_path, dpi=(dpi, dpi))
    print(f"Saved {output_path}  ({canvas_w}x{canvas_h} px, {dpi} DPI)")


# ── CLI ───────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Render a tool selection flowchart diagram from a JSON definition."
    )
    parser.add_argument("--input", required=True, help="JSON definition file")
    parser.add_argument("--output", required=True, help="Output PNG path")
    parser.add_argument("--dpi", type=int, default=150, help="Output DPI (default: 150)")
    args = parser.parse_args()

    with open(args.input) as f:
        definition = json.load(f)

    render(definition, args.output, args.dpi)


if __name__ == "__main__":
    main()
