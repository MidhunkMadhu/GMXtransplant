#!/usr/bin/env python3
"""Render the GMXtransplant Markdown documentation into one PDF.

Written against the PDF format directly so that building the documentation
needs nothing beyond the standard library: the documentation has to ship with
the package, and a build-time dependency on a converter would make that
unreliable on the machines this package is installed on.

Usage:
    python3 tools/build_docs_pdf.py [--output docs/GMXtransplant.pdf]
"""
from __future__ import annotations

import argparse
import re
import zlib
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# Order matters: this is the reading order of the finished document.
SOURCES = [
    ("Overview", "README.md"),
    ("Desktop application", "GUI_USAGE.md"),
    ("CHARMM-GUI transplant mode", "CHARMPROT.md"),
    ("Comparison views", "VISUALIZATION.md"),
    ("Portable and offline use", "PORTABLE_USAGE.md"),
]

PAGE_WIDTH, PAGE_HEIGHT = 595, 842
LEFT, RIGHT = 58, PAGE_WIDTH - 58
TOP, BOTTOM = PAGE_HEIGHT - 64, 68
BODY_WIDTH = RIGHT - LEFT

BODY_SIZE, BODY_LEAD = 9.8, 13.4
CODE_SIZE, CODE_LEAD = 8.2, 10.8

REGULAR, BOLD, MONO = b"/F1", b"/F2", b"/F3"

# Helvetica and Helvetica-Bold advance widths (1/1000 em) for ASCII 32-126.
_HELVETICA = (
    "278 278 355 556 556 889 667 191 333 333 389 584 278 333 278 278 "
    "556 556 556 556 556 556 556 556 556 556 278 278 584 584 584 556 "
    "1015 667 667 722 722 667 611 778 722 278 500 667 556 833 722 778 "
    "667 778 722 667 611 722 667 944 667 667 611 278 278 278 469 556 "
    "333 556 556 500 556 556 278 556 556 222 222 500 222 833 556 556 "
    "556 556 333 500 278 556 500 722 500 500 500 334 260 334 584"
)
_HELVETICA_BOLD = (
    "278 333 474 556 556 889 722 238 333 333 389 584 278 333 278 278 "
    "556 556 556 556 556 556 556 556 556 556 333 333 584 584 584 611 "
    "975 722 722 722 722 667 611 778 722 278 556 722 611 833 722 778 "
    "667 778 722 667 611 722 667 944 667 667 611 333 278 333 584 556 "
    "333 556 611 556 611 556 333 611 611 278 278 556 278 889 611 611 "
    "611 611 389 556 333 611 556 778 556 556 500 389 280 389 584"
)
WIDTHS = {
    REGULAR: [int(w) for w in _HELVETICA.split()],
    BOLD: [int(w) for w in _HELVETICA_BOLD.split()],
}

# Characters the docs use that Windows-1252 cannot represent.
TRANSLITERATE = {
    "≤": "<=", "≥": ">=", "→": "->", "←": "<-", "≈": "~", "∞": "inf",
    "‑": "-", "─": "-", "│": "|", "├": "+", "└": "+", "•": "•",
}


def encode(text: str) -> bytes:
    for source, target in TRANSLITERATE.items():
        text = text.replace(source, target)
    return text.encode("cp1252", "replace")


def escape(text: str) -> bytes:
    data = encode(text)
    for old, new in ((b"\\", b"\\\\"), (b"(", b"\\("), (b")", b"\\)")):
        data = data.replace(old, new)
    return data


def text_width(text: str, font: bytes, size: float) -> float:
    if font == MONO:
        return len(text) * 0.6 * size
    table = WIDTHS[font]
    total = 0
    for character in encode(text).decode("cp1252"):
        index = ord(character) - 32
        total += table[index] if 0 <= index < len(table) else 556
    return total * size / 1000.0


def wrap(text: str, font: bytes, size: float, width: float):
    """Greedy word wrap; a single over-long word is left to overhang."""
    lines, current = [], ""
    for word in text.split():
        candidate = f"{current} {word}".strip()
        if current and text_width(candidate, font, size) > width:
            lines.append(current)
            current = word
        else:
            current = candidate
    lines.append(current)
    return lines or [""]


# ---------------------------------------------------------------------------
# Markdown
# ---------------------------------------------------------------------------

INLINE = [
    (re.compile(r"!\[([^\]]*)\]\([^)]*\)"), r"\1"),
    (re.compile(r"\[([^\]]+)\]\(([^)]+)\)"), r"\1 (\2)"),
    (re.compile(r"\*\*\*(.+?)\*\*\*"), r"\1"),
    (re.compile(r"\*\*(.+?)\*\*"), r"\1"),
    (re.compile(r"(?<!\w)_(.+?)_(?!\w)"), r"\1"),
    (re.compile(r"`([^`]+)`"), r"\1"),
    (re.compile(r"<[^>]+>"), ""),
]


def inline(text: str) -> str:
    for pattern, replacement in INLINE:
        text = pattern.sub(replacement, text)
    return text.replace("&amp;", "&").replace("&lt;", "<").replace("&gt;", ">").strip()


def parse(markdown: str):
    """Reduce Markdown to the block kinds this document actually uses."""
    blocks, lines, index = [], markdown.splitlines(), 0
    paragraph: list[str] = []

    def flush():
        if paragraph:
            blocks.append(("paragraph", inline(" ".join(paragraph))))
            paragraph.clear()

    while index < len(lines):
        line = lines[index]
        stripped = line.strip()
        if stripped.startswith("```"):
            flush()
            index += 1
            code = []
            while index < len(lines) and not lines[index].strip().startswith("```"):
                code.append(lines[index].rstrip())
                index += 1
            blocks.append(("code", code))
            index += 1
            continue
        if not stripped:
            flush()
            index += 1
            continue
        heading = re.match(r"(#{1,6})\s+(.*)", stripped)
        if heading:
            flush()
            blocks.append(("heading", len(heading[1]), inline(heading[2])))
            index += 1
            continue
        if re.fullmatch(r"([-*_])\1{2,}", stripped):
            flush()
            blocks.append(("rule",))
            index += 1
            continue
        if stripped.startswith("|") and stripped.endswith("|"):
            flush()
            rows = []
            while index < len(lines) and lines[index].strip().startswith("|"):
                cells = [inline(c) for c in lines[index].strip().strip("|").split("|")]
                if not all(re.fullmatch(r":?-{2,}:?", c.strip() or "-") for c in cells):
                    rows.append(cells)
                index += 1
            blocks.append(("table", rows))
            continue
        bullet = re.match(r"(\s*)([-*+]|\d+[.)])\s+(.*)", line)
        if bullet:
            flush()
            marker = "•" if bullet[2] in "-*+" else bullet[2]
            blocks.append(("item", len(bullet[1]) // 2, marker, inline(bullet[3])))
            index += 1
            continue
        if stripped.startswith(">"):
            flush()
            blocks.append(("quote", inline(stripped.lstrip("> "))))
            index += 1
            continue
        paragraph.append(stripped)
        index += 1
    flush()
    return blocks


# ---------------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------------

class Document:
    def __init__(self):
        self.pages: list[list[bytes]] = []
        self.operations: list[bytes] = []
        self.y = 0.0
        self.entries: list[tuple[int, str, int, float]] = []
        self.break_before_heading = False

    def new_page(self):
        self.operations = []
        self.pages.append(self.operations)
        self.y = TOP

    def space(self, amount: float):
        if self.pages and self.y - amount > BOTTOM:
            self.y -= amount

    def room(self, height: float):
        if not self.pages or self.y - height < BOTTOM:
            self.new_page()

    def text(self, content: str, x: float, font: bytes, size: float, gray: float = 0.0):
        self.operations.append(
            b"BT %.3f g %b %.2f Tf 1 0 0 1 %.2f %.2f Tm (%b) Tj ET\n"
            % (gray, font, size, x, self.y, escape(content))
        )

    def rectangle(self, x, y, width, height, gray):
        self.operations.append(b"%.3f g %.2f %.2f %.2f %.2f re f\n" % (gray, x, y, width, height))

    def line(self, x1, y1, x2, y2, gray=0.75, thickness=0.6):
        self.operations.append(
            b"%.3f G %.2f w %.2f %.2f m %.2f %.2f l S\n" % (gray, thickness, x1, y1, x2, y2)
        )

    # -- block renderers ---------------------------------------------------

    def heading(self, level: int, title: str):
        if level == 1:
            self.new_page()
            size, lead = 19.0, 24.0
        else:
            size = {2: 14.0, 3: 11.2}.get(level, 10.2)
            lead = size + 4
            self.space(16 if level == 2 else 11)
            self.room(lead * 2 + BODY_LEAD)
        if level <= 3:
            self.entries.append((level, title, len(self.pages) - 1, self.y + lead))
        for line in wrap(title, BOLD, size, BODY_WIDTH):
            self.room(lead)
            self.y -= lead
            self.text(line, LEFT, BOLD, size, 0.09 if level > 1 else 0.05)
        if level == 1:
            self.y -= 8
            self.line(LEFT, self.y, RIGHT, self.y, 0.55, 1.0)
        self.y -= 8 if level == 1 else 5

    def paragraph(self, text: str, indent: float = 0.0, gray: float = 0.15):
        for line in wrap(text, REGULAR, BODY_SIZE, BODY_WIDTH - indent):
            self.room(BODY_LEAD)
            self.y -= BODY_LEAD
            self.text(line, LEFT + indent, REGULAR, BODY_SIZE, gray)
        self.y -= 3

    def item(self, depth: int, marker: str, text: str):
        indent = 14 + depth * 16
        hanging = indent + text_width(marker + " ", REGULAR, BODY_SIZE) + 3
        lines = wrap(text, REGULAR, BODY_SIZE, BODY_WIDTH - hanging)
        self.room(BODY_LEAD)
        self.y -= BODY_LEAD
        self.text(marker, LEFT + indent, REGULAR, BODY_SIZE, 0.35)
        self.text(lines[0], LEFT + hanging, REGULAR, BODY_SIZE, 0.15)
        for line in lines[1:]:
            self.room(BODY_LEAD)
            self.y -= BODY_LEAD
            self.text(line, LEFT + hanging, REGULAR, BODY_SIZE, 0.15)
        self.y -= 2

    def quote(self, text: str):
        top = self.y
        self.paragraph(text, indent=16, gray=0.35)
        self.line(LEFT + 5, self.y + 3, LEFT + 5, top - 2, 0.7, 2.0)

    def code(self, lines: list[str]):
        size = CODE_SIZE
        longest = max((text_width(line, MONO, size) for line in lines), default=0)
        while longest > BODY_WIDTH - 20 and size > 5.6:
            size -= 0.3
            longest = max((text_width(line, MONO, size) for line in lines), default=0)
        lead = size + 2.6
        self.space(6)
        for index, line in enumerate(lines):
            self.room(lead + 8)
            if index == 0 or self.y == TOP:
                self.rectangle(LEFT, self.y - lead * _remaining(lines, index, self.y, lead) - 5,
                               BODY_WIDTH, lead * _remaining(lines, index, self.y, lead) + 9, 0.945)
            self.y -= lead
            self.text(line.replace("\t", "    "), LEFT + 8, MONO, size, 0.12)
        self.y -= 9

    def table(self, rows: list[list[str]]):
        if not rows:
            return
        columns = max(len(row) for row in rows)
        rows = [row + [""] * (columns - len(row)) for row in rows]
        size = CODE_SIZE
        widths = [max(len(row[c]) for row in rows) for c in range(columns)]
        while sum(widths) * 0.6 * size + 3 * size * (columns - 1) > BODY_WIDTH and size > 5.4:
            size -= 0.3
        lead = size + 3.4
        self.space(6)
        for index, row in enumerate(rows):
            self.room(lead)
            self.y -= lead
            x = LEFT + 3
            font = BOLD if index == 0 else MONO
            for column, cell in enumerate(row):
                self.text(cell, x, font if index == 0 else MONO, size, 0.12)
                x += widths[column] * 0.6 * size + 3 * size
            if index == 0:
                self.line(LEFT, self.y - 3, RIGHT, self.y - 3, 0.7)
                self.y -= 3
        self.y -= 8

    def rule(self):
        self.space(8)
        self.room(10)
        self.line(LEFT, self.y, RIGHT, self.y, 0.85)
        self.y -= 8


def _remaining(lines, index, y, lead):
    """How many of the remaining code lines fit on this page (for the backdrop)."""
    fits = int((y - BOTTOM) // lead)
    return max(1, min(len(lines) - index, fits))


def render(document: Document, blocks):
    for block in blocks:
        kind = block[0]
        if kind == "heading":
            document.heading(block[1], block[2])
        elif kind == "paragraph":
            document.paragraph(block[1])
        elif kind == "item":
            document.item(block[1], block[2], block[3])
        elif kind == "code":
            document.code(block[1])
        elif kind == "table":
            document.table(block[1])
        elif kind == "quote":
            document.quote(block[1])
        elif kind == "rule":
            document.rule()


# ---------------------------------------------------------------------------
# PDF assembly
# ---------------------------------------------------------------------------

class PDF:
    def __init__(self):
        self.objects: list[bytes | None] = [None]

    def reserve(self) -> int:
        self.objects.append(None)
        return len(self.objects) - 1

    def add(self, data: bytes) -> int:
        number = self.reserve()
        self.objects[number] = data
        return number

    def stream(self, data: bytes) -> int:
        compressed = zlib.compress(data, 9)
        return self.add(b"<< /Length %d /Filter /FlateDecode >>\nstream\n%b\nendstream"
                        % (len(compressed), compressed))

    def build(self, root: int) -> bytes:
        out = bytearray(b"%PDF-1.4\n%\xe2\xe3\xcf\xd3\n")
        offsets = [0] * len(self.objects)
        for number in range(1, len(self.objects)):
            offsets[number] = len(out)
            out += b"%d 0 obj\n%b\nendobj\n" % (number, self.objects[number])
        start = len(out)
        out += b"xref\n0 %d\n0000000000 65535 f \n" % len(self.objects)
        for number in range(1, len(self.objects)):
            out += b"%010d 00000 n \n" % offsets[number]
        out += (b"trailer\n<< /Size %d /Root %d 0 R >>\nstartxref\n%d\n%%%%EOF\n"
                % (len(self.objects), root, start))
        return bytes(out)


def title_page(document: Document, title: str, subtitle: str, version: str):
    document.new_page()
    document.rectangle(0, PAGE_HEIGHT - 300, PAGE_WIDTH, 300, 0.93)
    document.y = PAGE_HEIGHT - 170
    document.text(title, LEFT, BOLD, 30, 0.05)
    document.y -= 34
    document.text(subtitle, LEFT, REGULAR, 15, 0.3)
    document.y -= 26
    document.text(f"Version {version}", LEFT, REGULAR, 11, 0.4)
    document.y = PAGE_HEIGHT - 360
    for line in [
        "Complete documentation for the command-line interface and the desktop",
        "application, assembled from the guides shipped with this release.",
        "",
        f"Generated {date.today().isoformat()}.",
    ]:
        document.y -= 16
        if line:
            document.text(line, LEFT, REGULAR, 10.5, 0.25)


def contents_page(document: Document, entries, offset: int):
    document.new_page()
    document.y -= 26
    document.text("Contents", LEFT, BOLD, 19, 0.05)
    document.y -= 12
    document.line(LEFT, document.y, RIGHT, document.y, 0.55, 1.0)
    document.y -= 14
    for level, title, page_index, _ in entries:
        if level > 2:
            continue
        size = 10.6 if level == 1 else 9.6
        indent = 0 if level == 1 else 18
        document.room(16)
        document.y -= 15 if level == 1 else 13
        font = BOLD if level == 1 else REGULAR
        number = str(page_index + offset + 1)
        label = title
        while (text_width(label, font, size) + text_width(number, REGULAR, size) + indent
               > BODY_WIDTH - 20 and len(label) > 8):
            label = label[:-2]
        document.text(label, LEFT + indent, font, size, 0.1 if level == 1 else 0.25)
        document.text(number, RIGHT - text_width(number, REGULAR, size), REGULAR, size, 0.3)
        dots_left = LEFT + indent + text_width(label, font, size) + 6
        dots_right = RIGHT - text_width(number, REGULAR, size) - 6
        if dots_right > dots_left:
            document.line(dots_left, document.y + 3, dots_right, document.y + 3, 0.82, 0.5)


def build(output: Path, version: str) -> Path:
    body = Document()
    for title, filename in SOURCES:
        source = ROOT / filename
        if not source.is_file():
            continue
        blocks = parse(source.read_text(encoding="utf-8"))
        # One top-level heading per source file, using our own section name.
        blocks = [b for b in blocks if not (b[0] == "heading" and b[1] == 1)]
        blocks = [("heading", 1, title)] + [
            (b[0], max(2, b[1]), *b[2:]) if b[0] == "heading" else b for b in blocks
        ]
        render(body, blocks)

    front = Document()
    title_page(front, "GMXTRANSPLANT", "MD System Builder", version)
    contents_page(front, body.entries, 0)
    # The contents page count is now known, so renumber against the real offset.
    front = Document()
    title_page(front, "GMXTRANSPLANT", "MD System Builder", version)
    contents_page(front, body.entries, 1)
    offset = len(front.pages)
    front = Document()
    title_page(front, "GMXTRANSPLANT", "MD System Builder", version)
    contents_page(front, body.entries, offset)

    pages = front.pages + body.pages
    pdf = PDF()
    page_numbers = [pdf.reserve() for _ in pages]
    pages_number = pdf.reserve()
    fonts = {
        b"/F1": pdf.add(b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica /Encoding /WinAnsiEncoding >>"),
        b"/F2": pdf.add(b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica-Bold /Encoding /WinAnsiEncoding >>"),
        b"/F3": pdf.add(b"<< /Type /Font /Subtype /Type1 /BaseFont /Courier /Encoding /WinAnsiEncoding >>"),
    }
    resources = pdf.add(b"<< /Font << %b >> >>" % b" ".join(
        b"%b %d 0 R" % (name, number) for name, number in fonts.items()))

    for index, operations in enumerate(pages):
        footer = []
        if index >= len(front.pages):
            label = f"{index + 1}"
            footer.append(b"%.3f G 0.5 w %.2f %.2f m %.2f %.2f l S\n"
                          % (0.85, LEFT, BOTTOM - 14, RIGHT, BOTTOM - 14))
            footer.append(b"BT 0.45 g /F1 8.4 Tf 1 0 0 1 %.2f %.2f Tm (%b) Tj ET\n"
                          % (LEFT, BOTTOM - 26, escape("GMXTRANSPLANT MD System Builder")))
            footer.append(b"BT 0.45 g /F1 8.4 Tf 1 0 0 1 %.2f %.2f Tm (%b) Tj ET\n"
                          % (RIGHT - text_width(label, REGULAR, 8.4), BOTTOM - 26, escape(label)))
        content = pdf.stream(b"".join(operations + footer))
        pdf.objects[page_numbers[index]] = (
            b"<< /Type /Page /Parent %d 0 R /MediaBox [0 0 %d %d] /Resources %d 0 R "
            b"/Contents %d 0 R >>" % (pages_number, PAGE_WIDTH, PAGE_HEIGHT, resources, content))

    pdf.objects[pages_number] = (
        b"<< /Type /Pages /Count %d /Kids [%b] >>"
        % (len(pages), b" ".join(b"%d 0 R" % number for number in page_numbers)))

    # Bookmarks, so every viewer offers the same navigation as the contents page.
    outline_root = pdf.reserve()
    items, previous_top = [], None
    tops = []
    for level, title, page_index, y in body.entries:
        if level > 2:
            continue
        number = pdf.reserve()
        items.append((number, level, title, page_index + len(front.pages), y))
    for position, (number, level, title, page_index, y) in enumerate(items):
        if level == 1:
            tops.append(position)
    for position, (number, level, title, page_index, y) in enumerate(items):
        parent = outline_root
        if level == 2:
            owners = [p for p in tops if p < position]
            parent = items[owners[-1]][0] if owners else outline_root
        siblings = [p for p, item in enumerate(items)
                    if item[1] == level and _same_parent(items, tops, p, position, level)]
        order = siblings.index(position)
        previous = b"/Prev %d 0 R " % items[siblings[order - 1]][0] if order else b""
        following = b"/Next %d 0 R " % items[siblings[order + 1]][0] if order + 1 < len(siblings) else b""
        children = [p for p, item in enumerate(items)
                    if item[1] == 2 and _same_parent(items, tops, p, position, 2) and level == 1]
        child = b""
        if level == 1 and children:
            child = (b"/First %d 0 R /Last %d 0 R /Count %d "
                     % (items[children[0]][0], items[children[-1]][0], len(children)))
        pdf.objects[number] = (
            b"<< /Title (%b) /Parent %d 0 R %b%b%b/Dest [%d 0 R /XYZ %.2f %.2f 0] >>"
            % (escape(title), parent, previous, following, child,
               page_numbers[page_index], float(LEFT - 10), y + 12))
    top_items = [items[p][0] for p in tops]
    pdf.objects[outline_root] = (
        b"<< /Type /Outlines /First %d 0 R /Last %d 0 R /Count %d >>"
        % (top_items[0], top_items[-1], len(top_items)) if top_items
        else b"<< /Type /Outlines /Count 0 >>")

    info = pdf.add(b"<< /Title (GMXTRANSPLANT MD System Builder) /Producer (GMXtransplant) >>")
    root = pdf.add(b"<< /Type /Catalog /Pages %d 0 R /Outlines %d 0 R /PageMode /UseOutlines >>"
                   % (pages_number, outline_root))
    data = pdf.build(root)
    # Info is referenced from the trailer, which build() writes without it.
    data = data.replace(b"/Root %d 0 R >>" % root, b"/Root %d 0 R /Info %d 0 R >>" % (root, info))
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_bytes(data)
    return output


def _same_parent(items, tops, candidate, reference, level):
    if level == 1:
        return True
    owner = [p for p in tops if p < reference]
    owner = owner[-1] if owner else -1
    other = [p for p in tops if p < candidate]
    other = other[-1] if other else -1
    return owner == other


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default=str(ROOT / "docs" / "GMXtransplant.pdf"))
    arguments = parser.parse_args()
    version = re.search(r'__version__\s*=\s*"([^"]+)"',
                        (ROOT / "gmxtransplant" / "__init__.py").read_text())[1]
    path = build(Path(arguments.output), version)
    print(f"Wrote {path} ({path.stat().st_size / 1024:.0f} kB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
