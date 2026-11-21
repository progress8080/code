"""Build the bachelor's thesis as PDF and Word (DOCX).

Run from the project root:
    python thesis/build_thesis.py               -> thesis/output/thesis.pdf and thesis.docx
    python thesis/build_thesis.py --no-docx     -> PDF only
    python thesis/build_thesis.py --keep-build  -> keep the intermediate files in thesis/_build
    python thesis/build_thesis.py --allow-missing  -> skip chapters that do not exist yet

Requirements: pandoc ≥ 3 (with built-in citeproc), Microsoft Edge, Google Chrome or
Chromium (headless printing), and the Python packages pyyaml, matplotlib and pypdfium2.

How the PDF is produced
  1. The main matter (chapters, references, appendices) is converted to HTML and
     printed to PDF with Arabic page numbers starting at 1.
  2. The page of every chapter, section, figure and table is read back from that
     PDF to fill in the table of contents and the lists of figures and tables.
  3. The front matter (title page, declaration, abstract, acknowledgements,
     contents, lists, nomenclature) is printed with Roman page numbers, twice,
     so that its own entries in the contents also get correct page numbers.
  4. Both PDFs are merged into thesis.pdf.
"""

from __future__ import annotations

import html
import os
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

import yaml

THESIS = Path(__file__).resolve().parent
ROOT = THESIS.parent
OUT = THESIS / "output"
BUILD = THESIS / "_build"

CHAPTERS = [
    "chapters/01_introduction.md",
    "chapters/02_background.md",
    "chapters/03_methodology.md",
    "chapters/04_implementation.md",
    "chapters/05_verification.md",
    "chapters/06_results.md",
    "chapters/07_conclusions.md",
]
FRONT = ["front/declaration.md", "front/abstract.md", "front/acknowledgements.md"]
NOMENCLATURE = "front/nomenclature.md"

# Source files listed in Appendix B, in reading order.
LISTED = [
    ("Project configuration", ["pyproject.toml"]),
    ("Calculation package", [
        "src/shipsupport/__init__.py", "src/shipsupport/__main__.py", "src/shipsupport/units.py",
        "src/shipsupport/model.py", "src/shipsupport/loader.py", "src/shipsupport/weights.py",
        "src/shipsupport/deformation.py",
        "src/shipsupport/solvers/__init__.py", "src/shipsupport/solvers/contact.py",
        "src/shipsupport/solvers/result.py", "src/shipsupport/solvers/rigid.py",
        "src/shipsupport/solvers/fem.py", "src/shipsupport/strength.py",
        "src/shipsupport/spacing.py",
        "src/shipsupport/analysis.py", "src/shipsupport/report.py", "src/shipsupport/cli.py"]),
    ("Graphical user interface", [
        "src/shipsupport_gui/__init__.py", "src/shipsupport_gui/__main__.py",
        "src/shipsupport_gui/app.py", "src/shipsupport_gui/style.py",
        "src/shipsupport_gui/widgets.py",
        "src/shipsupport_gui/plan_panel.py", "src/shipsupport_gui/input_panels.py",
        "src/shipsupport_gui/results_panels.py",
        "src/shipsupport_gui/worker.py", "src/shipsupport_gui/main_window.py"]),
    ("Automated tests", [
        "tests/helpers.py", "tests/test_weights.py", "tests/test_loader.py", "tests/test_rigid.py",
        "tests/test_fem.py", "tests/test_strength.py", "tests/test_deformation.py",
        "tests/test_spacing.py", "tests/test_example_case.py", "tests/test_gui.py"]),
]

CSS = r"""
@page { size: A4; margin: 25mm 20mm 25mm 30mm; }
@page main  { @bottom-center { content: counter(page); font: 10pt "Times New Roman", serif; } }
@page front { @bottom-center { content: counter(page, lower-roman);
                               font: 10pt "Times New Roman", serif; } }
@page title { @bottom-center { content: none; } }
body { font-family: "Times New Roman", Times, serif; font-size: 12pt; line-height: 1.5;
       text-align: justify; hyphens: auto; color: #000; margin: 0; }
body.main { page: main; }
body.front { page: front; }
.titlepage { page: title; text-align: center; break-after: page; line-height: 1.35; }
.titlepage .university { font-size: 16pt; font-weight: bold; margin-top: 10mm; }
.titlepage .department { font-size: 13pt; margin-top: 2mm; }
.titlepage .title { font-size: 21pt; font-weight: bold; margin-top: 38mm; line-height: 1.25; }
.titlepage .subtitle { font-size: 13.5pt; font-style: italic; margin-top: 6mm; }
.titlepage .kind { font-size: 12pt; margin-top: 26mm; }
.titlepage .degree { font-size: 13pt; font-weight: bold; margin-top: 3mm; }
.titlepage .by { margin-top: 18mm; font-size: 12pt; }
.titlepage .author { font-size: 15pt; font-weight: bold; }
.titlepage table { margin: 16mm auto 0; border: none; font-size: 12pt; }
.titlepage td { border: none; padding: 1mm 4mm; text-align: left; }
.titlepage .placedate { margin-top: 18mm; font-size: 12pt; }
h1 { font-size: 18pt; break-before: page; margin: 0 0 8mm; padding-top: 6mm; line-height: 1.25; }
h2 { font-size: 14pt; margin: 7mm 0 2mm; break-after: avoid; }
h3 { font-size: 12.5pt; margin: 5mm 0 1.5mm; break-after: avoid; }
h4 { font-size: 12pt; font-style: italic; break-after: avoid; }
p { margin: 0 0 3mm; orphans: 3; widows: 3; }
figure { margin: 5mm 0; text-align: center; break-inside: avoid; }
figure img { max-width: 100%; max-height: 190mm; }
figcaption { font-size: 10.5pt; margin-top: 1.5mm; text-align: center; }
table { border-collapse: collapse; margin: 4mm auto; font-size: 10.5pt; line-height: 1.3;
        break-inside: auto; }
caption { caption-side: top; font-size: 10.5pt; margin-bottom: 1.5mm; text-align: center; }
th { border-top: 1.5px solid #000; border-bottom: 1px solid #000; padding: 1mm 2.5mm; }
td { padding: .8mm 2.5mm; border-bottom: .5px solid #bbb; }
tbody tr:last-child td { border-bottom: 1.5px solid #000; }
tr { break-inside: avoid; }
pre { font-family: Consolas, "Courier New", monospace; font-size: 8pt; line-height: 1.3;
      white-space: pre-wrap; word-break: break-all; background: #f6f6f6; padding: 2mm 3mm;
      border-left: 2px solid #999; text-align: left; }
code { font-family: Consolas, "Courier New", monospace; font-size: .9em; }
.equation { display: flex; align-items: center; break-inside: avoid; margin: 2mm 0; }
.equation > p { flex: 1 1 auto; text-align: center; margin: 0; }
.equation > p:last-child { flex: 0 0 16mm; text-align: right; }
math[display="block"] { margin: 1.5mm 0; }
blockquote { margin: 3mm 8mm; font-style: italic; }
a { color: inherit; text-decoration: none; }
.toc { list-style: none; padding: 0; margin: 0; }
.toc li { display: flex; align-items: baseline; margin: 1mm 0; }
.toc li .t { flex: 0 1 auto; }
.toc li .dots { flex: 1; border-bottom: 1px dotted #555; margin: 0 2mm; transform: translateY(-1mm); }
.toc li .pg { flex: 0 0 auto; }
.toc li.l1 { font-weight: bold; margin-top: 3mm; }
.toc li.l2 { padding-left: 8mm; }
.toc li.lf { padding-left: 0; }
#refs .csl-entry { padding-left: 8mm; text-indent: -8mm; margin-bottom: 2mm; text-align: left; }
.codeinfo { font-size: 10pt; color: #333; }
"""


@dataclass
class Entry:
    level: int          # 1, 2 for headings; 0 for figure/table list entries
    text: str
    needle: str         # normalised text searched for in the PDF
    page: str = ""


# --------------------------------------------------------------------------- helpers

def _norm(text: str) -> str:
    text = re.sub(r"\$[^$]*\$", "", text)                 # drop inline math
    text = re.sub(r"[*_`\\]", "", text)
    text = text.replace("–", "-").replace("—", "-").replace("−", "-")
    return re.sub(r"\s+", "", text)


def _strip_attrs(heading: str) -> str:
    return re.sub(r"\s*\{[^}]*\}\s*$", "", heading).strip()


def _md_inline_to_html(text: str) -> str:
    """Tiny inline converter for list entries (math dropped to plain text)."""
    text = re.sub(r"\$([^$]*)\$", lambda m: m.group(1).replace("\\", ""), text)
    text = text.replace("`", "")
    text = html.escape(text)
    return re.sub(r"\*\*(.+?)\*\*|\*(.+?)\*", lambda m: m.group(1) or m.group(2), text)


def _iter_outside_code(lines):
    fence = None
    for line in lines:
        stripped = line.lstrip()
        match = re.match(r"(`{3,}|~{3,})", stripped)
        if match:
            if fence is None:
                fence = match.group(1)
            elif stripped.startswith(fence):
                fence = None
            yield line, True
            continue
        yield line, fence is not None


def equations(md: str, target: str) -> str:
    """Turn  $$ ... \\tag{3.1} $$  into a numbered equation for the given target.

    Only text outside code blocks is converted, so code listings stay verbatim.
    """
    def repl(match):
        body = match.group(1)
        tag = re.search(r"\\tag\{([^}]*)\}", body)
        if not tag:
            return match.group(0)
        body = body.replace(tag.group(0), "").strip()
        number = tag.group(1)
        if target == "docx":
            return f"$$ {body} \\qquad\\qquad ({number}) $$"
        return f'::: {{.equation}}\n$$ {body} $$\n\n[({number})]{{.eqno}}\n:::'

    out, chunk = [], []

    def flush():
        if chunk:
            out.append(re.sub(r"\$\$(.*?)\$\$", repl, "\n".join(chunk), flags=re.S))
            chunk.clear()

    for line, in_code in _iter_outside_code(md.split("\n")):
        if in_code:
            flush()
            out.append(line)
        else:
            chunk.append(line)
    flush()
    return "\n".join(out)


# --------------------------------------------------------------------------- appendices

# Folder tree in Appendix B: the software only — the thesis folder itself and generated,
# cached and output files are left out.
SKIP_DIRS = {"__pycache__", ".pytest_cache", ".git", ".vscode", ".idea", "node_modules",
             "build", "dist", "results", "_build", "output", "thesis"}
SKIP_DIR_SUFFIXES = ("_results", ".egg-info")
SKIP_FILE_SUFFIXES = (".pyc",)
LANGUAGES = {".py": "python", ".toml": "toml", ".yaml": "yaml", ".yml": "yaml"}
DESCRIPTIONS = {
    "pyproject.toml": "package metadata, dependencies, entry points",
    ".gitignore": "files not tracked by version control",
    "src": "source code",
    "shipsupport": "calculation package",
    "solvers": "rigid-hull and FEM solvers, contact iteration",
    "shipsupport_gui": "PyQt5 desktop application",
    "examples": "example cases (placeholder data)",
    "tests": "automated test suite",
    "README.md": "installation and usage",
}


def _skip(path: Path) -> bool:
    if path.is_dir():
        return (path.name in SKIP_DIRS or path.name.endswith(SKIP_DIR_SUFFIXES)
                or path.name.startswith("pytest-cache-files"))
    return path.name.endswith(SKIP_FILE_SUFFIXES)


def folder_tree(root: Path) -> str:
    lines = [f"{root.name}/"]

    def walk(folder: Path, prefix: str):
        entries = sorted((p for p in folder.iterdir() if not _skip(p)),
                         key=lambda p: (p.is_file(), p.name.lower()))
        for i, entry in enumerate(entries):
            last = i == len(entries) - 1
            label = f"{prefix}{'└── ' if last else '├── '}{entry.name}{'/' if entry.is_dir() else ''}"
            note = DESCRIPTIONS.get(entry.name)
            lines.append(f"{label:<58}# {note}" if note else label)
            if entry.is_dir():
                walk(entry, prefix + ("    " if last else "│   "))

    walk(root, "")
    return "\n".join(lines)


def _fence(text: str) -> str:
    """A code fence longer than any backtick run inside the listed file."""
    longest = max((len(m) for m in re.findall(r"`{3,}", text)), default=0)
    return "`" * max(3, longest + 1)


def find_browser() -> str | None:
    for name in ("msedge", "chrome", "google-chrome", "chromium", "chromium-browser"):
        if shutil.which(name):
            return shutil.which(name)
    candidates = [
        r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
        r"C:\Program Files\Microsoft\Edge\Application\msedge.exe",
        r"C:\Program Files\Google\Chrome\Application\chrome.exe",
        "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
        "/Applications/Microsoft Edge.app/Contents/MacOS/Microsoft Edge",
    ]
    return next((c for c in candidates if Path(c).exists()), None)


def appendix_a() -> str:
    base = (ROOT / "examples" / "example_case.yaml").read_text(encoding="utf-8").rstrip()
    welding = (ROOT / "examples" / "welding_case.yaml").read_text(encoding="utf-8")
    marker = "# Welding deformation of the assembled hull"
    start = welding.find(marker)
    extra = welding[welding.rfind("\n# ----", 0, start) + 1:].rstrip() if start >= 0 else welding
    return (
        "# Appendix A – Example Input Files\n\n"
        "All values in these files are illustrative placeholders chosen to demonstrate the "
        "methods; they are not data of a real ship. Lengths are in metres from the aft "
        "perpendicular, weights in tonnes, moments in tonne-metres, stresses in MPa and gaps "
        "in millimetres, "
        "as indicated by the key suffixes.\n\n"
        "## A.1 Example vessel (`examples/example_case.yaml`)\n\n"
        f"```yaml\n{base}\n```\n\n"
        "## A.2 Welded-hull case (`examples/welding_case.yaml`)\n\n"
        "The welded-hull case uses the same ship, blocks and settings as A.1 and adds the "
        "following sections.\n\n"
        f"```yaml\n{extra}\n```\n\n"
        "## A.3 Spacing design for a block plan (`examples/spacing_case.yaml`)\n\n"
        "The block plan and the spacing-design settings of the main design example; the ship, "
        "weights, sections, structure and limits are those of A.1.\n\n"
        f"```yaml\n{_section(ROOT / 'examples' / 'spacing_case.yaml', 'blocks:')}\n```\n\n"
        "## A.4 Weight-distribution curve (`examples/weight_curve_case.yaml`)\n\n"
        "The lightship weight given as a piecewise-linear curve and deadweight items given with "
        "their longitudinal moments.\n\n"
        f"```yaml\n{_section(ROOT / 'examples' / 'weight_curve_case.yaml', 'weights:', 'weight_curve:')}\n```\n")


def _section(path: Path, *keys: str) -> str:
    """The top-level YAML sections starting with the given keys (with their comments)."""
    lines = path.read_text(encoding="utf-8").splitlines()
    out, keep = [], False
    for line in lines:
        if line and not line[0].isspace() and not line.startswith("#"):
            keep = line.startswith(keys) or line.startswith("spacing_design:") and "blocks:" in keys
        if keep:
            out.append(line)
    return "\n".join(out).rstrip()


def appendix_b() -> str:
    parts = [
        "# Appendix B – Source Code\n",
        "This appendix lists the folder structure of the software developed in this thesis and "
        "the complete source code of the calculation package, the graphical user interface and "
        "the automated tests. The listing is generated from the files on disk when the thesis is "
        "built, so it always corresponds to the version of the software described in Chapters 4 "
        "and 5. The thesis sources themselves and generated files are omitted.\n",
        "## B.1 Folder Structure\n",
        "```text\n" + folder_tree(ROOT) + "\n```\n",
    ]
    total = 0
    for s, (section, files) in enumerate(LISTED, start=2):
        parts.append(f"## B.{s} {section}\n")
        for i, rel in enumerate(files, start=1):
            text = (ROOT / rel).read_text(encoding="utf-8").rstrip("\n")
            n = text.count("\n") + 1 if text else 0
            total += n
            fence = _fence(text)
            language = LANGUAGES.get(Path(rel).suffix, "text")
            parts.append(f"### B.{s}.{i} {rel}\n")
            parts.append(f"[{n} lines]{{.codeinfo}}\n")
            parts.append(f"{fence}{language}\n{text}\n{fence}\n")
    parts.insert(2, f"In total the listed files contain {total} lines.\n")
    return "\n".join(parts)


# --------------------------------------------------------------------------- assembly

ALLOW_MISSING = "--allow-missing" in sys.argv


def main_markdown(target: str) -> str:
    parts = []
    for chapter in CHAPTERS:
        path = THESIS / chapter
        if not path.exists():
            if not ALLOW_MISSING:
                raise SystemExit(f"Missing chapter: {chapter}")
            print(f"warning: skipping missing chapter {chapter}")
            continue
        parts.append(path.read_text(encoding="utf-8").strip())
    parts.append("# References {.unnumbered}\n\n::: {#refs}\n:::")
    parts.append(appendix_a())
    parts.append(appendix_b())
    md = "\n\n".join(parts) + "\n"
    return equations(md, target)


def collect_entries(md: str) -> tuple[list[Entry], list[Entry], list[Entry]]:
    toc, figures, tables = [], [], []
    for line, in_code in _iter_outside_code(md.splitlines()):
        if in_code:
            continue
        heading = re.match(r"^(#{1,2})\s+(.+)$", line)
        if heading:
            text = _strip_attrs(heading.group(2))
            toc.append(Entry(len(heading.group(1)), text, _norm(text)[:40]))
            continue
        fig = re.match(r"^!\[(Figure\s+[A-Z0-9]+\.\d+):\s*(.+?)\]\(", line)
        if fig:
            figures.append(Entry(0, f"{fig.group(1)}: {fig.group(2)}", _norm(fig.group(1) + ":")))
            continue
        tab = re.match(r"^Table:\s*(Table\s+[A-Z0-9]+\.\d+)\s*[–—-]\s*(.+)$", line)
        if tab:
            tables.append(Entry(0, f"{tab.group(1)} – {_strip_attrs(tab.group(2))}",
                                _norm(tab.group(1) + "–")))
    return toc, figures, tables


def find_pages(pdf_path: Path, entries: list[Entry], first_page: int = 0, roman: bool = False):
    """Fill entry.page with the printed page number where each entry first appears."""
    import pypdfium2 as pdfium

    pdf = pdfium.PdfDocument(str(pdf_path))
    pages = [_norm(pdf[i].get_textpage().get_text_range()) for i in range(len(pdf))]
    pdf.close()
    cursor = first_page
    missing = []
    for entry in entries:
        for i in range(cursor, len(pages)):
            if entry.needle and entry.needle in pages[i]:
                entry.page = _roman(i + 1) if roman else str(i + 1)
                cursor = i
                break
        else:
            missing.append(entry.text)
    if missing:
        print(f"warning: page not found for {len(missing)} entries, e.g. {missing[:3]}")
    return len(pages)


def _roman(n: int) -> str:
    out = ""
    for value, numeral in ((1000, "m"), (900, "cm"), (500, "d"), (400, "cd"), (100, "c"),
                           (90, "xc"), (50, "l"), (40, "xl"), (10, "x"), (9, "ix"), (5, "v"),
                           (4, "iv"), (1, "i")):
        while n >= value:
            out += numeral
            n -= value
    return out


def _list_html(entries: list[Entry], cls: str) -> str:
    rows = []
    for e in entries:
        level = {1: "l1", 2: "l2"}.get(e.level, "lf")
        rows.append(f'<li class="{level}"><span class="t">{_md_inline_to_html(e.text)}</span>'
                    f'<span class="dots"></span><span class="pg">{e.page}</span></li>')
    return f'<ul class="toc {cls}">\n' + "\n".join(rows) + "\n</ul>"


def title_page_html(meta: dict) -> str:
    rows = [("Student ID:", meta.get("student_id")), ("Supervisor:", meta.get("supervisor")),
            ("Co-supervisor:", meta.get("co_supervisor"))]
    rows = "".join(f"<tr><td>{html.escape(k)}</td><td>{html.escape(str(v))}</td></tr>"
                   for k, v in rows if v)
    e = lambda key: html.escape(str(meta.get(key, "")))  # noqa: E731
    return f"""<section class="titlepage">
<div class="university">{e('university')}</div>
<div class="department">{e('department')}</div>
<div class="title">{e('title')}</div>
<div class="subtitle">{e('subtitle')}</div>
<div class="kind">A thesis submitted in partial fulfilment of the requirements for the degree of</div>
<div class="degree">{e('degree')}</div>
<div class="by">by</div>
<div class="author">{e('author')}</div>
<table>{rows}</table>
<div class="placedate">{e('place')}, {e('date')}</div>
</section>"""


def front_markdown(toc_html: str, lof_html: str, lot_html: str) -> str:
    parts = [(THESIS / f).read_text(encoding="utf-8").strip() for f in FRONT]
    parts.append("# Table of Contents {.unnumbered}\n\n```{=html}\n" + toc_html + "\n```")
    parts.append("# List of Figures {.unnumbered}\n\n```{=html}\n" + lof_html + "\n```")
    parts.append("# List of Tables {.unnumbered}\n\n```{=html}\n" + lot_html + "\n```")
    parts.append((THESIS / NOMENCLATURE).read_text(encoding="utf-8").strip())
    return "\n\n".join(parts) + "\n"


# --------------------------------------------------------------------------- tools

def pandoc(md: str, output: Path, fmt: str, extra: list[str], body_class: str | None = None):
    BUILD.mkdir(exist_ok=True)
    source = BUILD / (output.stem + ".md")
    source.write_text(md, encoding="utf-8")
    args = ["pandoc", source.name, "-o", str(output), "--from", "markdown+tex_math_dollars",
            "--citeproc", "--bibliography", str(THESIS / "references.bib"),
            "--resource-path", f"{THESIS}{';' if sys.platform == 'win32' else ':'}{BUILD}", *extra]
    if fmt == "html":
        css = BUILD / "thesis.css"
        css.write_text(CSS, encoding="utf-8")
        args += ["--to", "html5", "--standalone", "--mathml", "--embed-resources",
                 "--css", str(css), "--metadata", "pagetitle=Thesis",
                 "--metadata", "link-citations=true"]
    subprocess.run(args, cwd=BUILD, check=True)
    if fmt == "html" and body_class:
        text = output.read_text(encoding="utf-8")
        output.write_text(text.replace("<body>", f'<body class="{body_class}">', 1), encoding="utf-8")


def print_pdf(html_path: Path, pdf_path: Path):
    browser = find_browser()
    if browser is None:
        raise SystemExit("Edge/Chrome/Chromium not found for PDF printing")
    pdf_path.unlink(missing_ok=True)           # never reuse a PDF from an earlier build
    with tempfile.TemporaryDirectory() as profile:
        subprocess.run([browser, "--headless=new", "--disable-gpu", "--no-first-run",
                        f"--user-data-dir={profile}", "--no-pdf-header-footer",
                        "--run-all-compositor-stages-before-draw", "--virtual-time-budget=60000",
                        f"--print-to-pdf={pdf_path}", html_path.as_uri()],
                       check=True, timeout=900, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if not pdf_path.exists():
        raise SystemExit(f"PDF printing failed: {pdf_path}")


def merge_pdfs(parts: list[Path], output: Path, meta: dict):
    import pypdfium2 as pdfium

    merged = pdfium.PdfDocument.new()
    for part in parts:
        doc = pdfium.PdfDocument(str(part))
        merged.import_pages(doc)
        doc.close()
    merged.save(str(output))
    merged.close()


# --------------------------------------------------------------------------- docx

def build_docx(meta: dict, figures: list[Entry], tables: list[Entry]):
    page_break = '```{=openxml}\n<w:p><w:r><w:br w:type="page"/></w:r></w:p>\n```'
    toc_field = ('```{=openxml}\n<w:p><w:r><w:fldChar w:fldCharType="begin" w:dirty="true"/></w:r>'
                 '<w:r><w:instrText xml:space="preserve">TOC \\o "1-2" \\h \\z \\u</w:instrText></w:r>'
                 '<w:r><w:fldChar w:fldCharType="separate"/></w:r><w:r><w:t>Right-click and choose '
                 '"Update Field" to build the table of contents.</w:t></w:r>'
                 '<w:r><w:fldChar w:fldCharType="end"/></w:r></w:p>\n```')
    t = lambda k: str(meta.get(k, ""))  # noqa: E731
    title = "\n\n".join([
        f'::: {{custom-style="Title"}}\n{t("title")}\n:::',
        f'::: {{custom-style="Subtitle"}}\n{t("subtitle")}\n:::',
        f"**{t('university')}**  \n{t('department')}",
        f"A thesis submitted in partial fulfilment of the requirements for the degree of  \n"
        f"**{t('degree')}**",
        f"by **{t('author')}** ({t('student_id')})",
        f"Supervisor: {t('supervisor')}  \nCo-supervisor: {t('co_supervisor')}",
        f"{t('place')}, {t('date')}",
    ])
    parts = [title]
    for f in FRONT:
        parts += [page_break, (THESIS / f).read_text(encoding="utf-8").strip()]
    parts += [page_break, "# Table of Contents {.unnumbered}", toc_field]
    parts += [page_break, "# List of Figures {.unnumbered}",
              "\n".join(f"- {e.text}" for e in figures)]
    parts += [page_break, "# List of Tables {.unnumbered}",
              "\n".join(f"- {e.text}" for e in tables)]
    parts += [page_break, (THESIS / NOMENCLATURE).read_text(encoding="utf-8").strip()]
    body_lines = []
    for line, in_code in _iter_outside_code(main_markdown("docx").splitlines()):
        if not in_code and line.startswith("# "):     # chapter heading, not a code comment
            body_lines += [page_break, ""]
        body_lines.append(line)
    body = "\n".join(body_lines)
    md = "\n\n".join(parts) + "\n\n" + body
    pandoc(md, OUT / "thesis.docx", "docx", ["--metadata", "link-citations=true"])


# --------------------------------------------------------------------------- main

def main(argv: list[str]) -> int:
    if shutil.which("pandoc") is None:
        print("pandoc is required: https://pandoc.org", file=sys.stderr)
        return 1
    try:
        import pypdfium2  # noqa: F401  (page numbers and merging of the PDF parts)
    except ImportError:
        print("pypdfium2 is required: pip install pypdfium2", file=sys.stderr)
        return 1
    meta = yaml.safe_load((THESIS / "metadata.yaml").read_text(encoding="utf-8"))
    OUT.mkdir(exist_ok=True)
    BUILD.mkdir(exist_ok=True)

    env = {**os.environ, "PYTHONDONTWRITEBYTECODE": "1"}   # keep the project free of caches
    figure_args = [sys.executable, str(THESIS / "make_thesis_figures.py")]
    subprocess.run(figure_args + (["--no-gui"] if "--no-gui" in argv else []), check=True, env=env)

    # 1. main matter
    md = main_markdown("html")
    toc, figures, tables = collect_entries(md)
    main_html, main_pdf = BUILD / "main.html", BUILD / "main.pdf"
    pandoc(md, main_html, "html", [], body_class="main")
    print_pdf(main_html, main_pdf)
    find_pages(main_pdf, toc)
    find_pages(main_pdf, figures)
    find_pages(main_pdf, tables)

    # 2. front matter, printed twice so its own contents entries get page numbers
    front_entries = [Entry(1, "Declaration of Authorship", _norm("Declaration of Authorship")),
                     Entry(1, "Abstract", _norm("Abstract")),
                     Entry(1, "Acknowledgements", _norm("Acknowledgements")),
                     Entry(1, "List of Figures", _norm("List of Figures")),
                     Entry(1, "List of Tables", _norm("List of Tables")),
                     Entry(1, "Nomenclature", _norm("Nomenclature"))]
    front_html, front_pdf = BUILD / "front.html", BUILD / "front.pdf"
    for _ in range(2):
        toc_html = _list_html(front_entries + toc, "contents")
        md_front = front_markdown(toc_html, _list_html(figures, "lof"), _list_html(tables, "lot"))
        pandoc(md_front, front_html, "html", ["--include-before-body",
                                             str(_write(BUILD / "title.html", title_page_html(meta)))],
               body_class="front")
        print_pdf(front_html, front_pdf)
        _locate_front(front_pdf, front_entries)

    merge_pdfs([front_pdf, main_pdf], OUT / "thesis.pdf", meta)
    print("wrote", OUT / "thesis.pdf")

    if "--no-docx" not in argv:
        build_docx(meta, figures, tables)
        print("wrote", OUT / "thesis.docx")
    if "--keep-build" not in argv:              # intermediate files are not needed afterwards
        shutil.rmtree(BUILD, ignore_errors=True)
    return 0


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def _locate_front(pdf_path: Path, entries: list[Entry]):
    """Front-matter headings sit at the top of their page; skip the contents pages."""
    import pypdfium2 as pdfium

    pdf = pdfium.PdfDocument(str(pdf_path))
    texts = [_norm(pdf[i].get_textpage().get_text_range()) for i in range(len(pdf))]
    pdf.close()
    for entry in entries:
        for i, text in enumerate(texts):
            # extracted text starts with the page number (e.g. "iii"), then the heading
            if 0 <= text[:12 + len(entry.needle)].find(entry.needle) <= len(_roman(i + 1)):
                entry.page = _roman(i + 1)
                break


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
