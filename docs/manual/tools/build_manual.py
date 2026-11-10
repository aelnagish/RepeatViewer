"""Build the RepeatViewer user manual from the Markdown chapters.

    python docs/manual/tools/build_manual.py            web pages + in-app copy
    python docs/manual/tools/build_manual.py --pdf      also the PDF

Sources: docs/manual/NN-name.md (one chapter each) and docs/manual/img/.
Outputs:
    docs/manual/*.html                          the web manual (GitHub Pages)
    repeatviewer/assets/manuals/user_manual/    the same pages for Help in the app
    docs/manual/RepeatViewer_User_Manual.pdf    with --pdf (copied into the app too)

Needs the "markdown" package (pip install markdown). The PDF needs Playwright
with Chromium (pip install playwright; playwright install chromium).
"""

from __future__ import annotations

import argparse
import html
import json
import re
import shutil
import struct
import sys
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "docs" / "manual"
IMG = SRC / "img"
APP_OUT = ROOT / "repeatviewer" / "assets" / "manuals" / "user_manual"
PDF_NAME = "RepeatViewer_User_Manual.pdf"
TITLE = "RepeatViewer User Manual"
AUTHORS = "Abdullah Yusuf, Nusrat Sultana, Jiangzhao Qian, Tony Heitkam"
AFFILIATION = "Molecular Botany, RWTH Aachen University, Aachen, Germany"

try:
    import markdown
except ImportError:  # pragma: no cover - build-time tool
    sys.exit("The manual build needs the markdown package: pip install markdown")


def version() -> str:
    text = (ROOT / "repeatviewer" / "__init__.py").read_text(encoding="utf-8")
    m = re.search(r'__version__\s*=\s*"([^"]+)"', text)
    return m.group(1) if m else ""


def png_size(path: Path) -> tuple[int, int] | None:
    try:
        with path.open("rb") as fh:
            head = fh.read(24)
        if head[:8] != b"\x89PNG\r\n\x1a\n":
            return None
        return struct.unpack(">II", head[16:24])
    except OSError:
        return None


def slug(text: str, sep: str = "-") -> str:
    text = re.sub(r"<[^>]+>", "", text)
    text = html.unescape(text).lower()
    text = re.sub(r"[^\w\s-]", "", text)
    return re.sub(r"[\s_-]+", sep, text).strip(sep)


class Chapter:
    def __init__(self, path: Path):
        self.path = path
        self.num = int(path.name[:2])
        self.stem = path.stem
        self.text = path.read_text(encoding="utf-8")
        m = re.search(r"^# (.+)$", self.text, re.M)
        self.title = m.group(1).strip() if m else self.stem
        self.html = ""
        self.sections: list[tuple[str, str, str]] = []  # (number, id, title)

    @property
    def page(self) -> str:
        return f"{self.stem}.html"


def render(ch: Chapter, *, print_mode: bool) -> str:
    md = markdown.Markdown(extensions=["tables", "fenced_code", "toc", "attr_list", "sane_lists"],
                           extension_configs={"toc": {"slugify": lambda v, s: slug(v)}})
    body = md.convert(ch.text)
    prefix = f"ch{ch.num:02d}-" if print_mode else ""

    # Chapter title.
    body = re.sub(r"<h1 id=\"([^\"]+)\">(.*?)</h1>",
                  lambda m: (f'<h1 id="{prefix}{m.group(1)}" class="chapter-title">'
                             f'<span class="chapter-label">Chapter {ch.num}</span>'
                             f'{m.group(2)}</h1>'), body, count=1)

    # Number the sections (h2) and keep them for the tables of contents.
    counter = {"h2": 0}
    ch.sections = []

    def number_h2(m):
        counter["h2"] += 1
        num = f"{ch.num}.{counter['h2']}"
        ch.sections.append((num, m.group(1), re.sub(r"<[^>]+>", "", m.group(2))))
        return (f'<h2 id="{prefix}{m.group(1)}"><span class="secnum">{num}</span>'
                f'{m.group(2)}</h2>')

    body = re.sub(r"<h2 id=\"([^\"]+)\">(.*?)</h2>", number_h2, body)
    if prefix:
        body = re.sub(r"<h([3-6]) id=\"([^\"]+)\">", rf'<h\1 id="{prefix}\2">', body)

    # Figures: a paragraph holding only an image becomes a numbered figure.
    fig = {"n": 0}

    def figure(m):
        fig["n"] += 1
        alt, src = html.unescape(m.group(1)), m.group(2)
        size = png_size(SRC / src)
        width = f' width="{size[0] // 2}"' if size and size[0] > 400 else (
            f' width="{size[0]}"' if size else "")
        return (f'<figure id="{prefix}fig-{ch.num}-{fig["n"]}">'
                f'<img src="{src}" alt="{html.escape(alt)}"{width} loading="lazy">'
                f'<figcaption><b>Figure {ch.num}.{fig["n"]}.</b> {html.escape(alt)}</figcaption>'
                f'</figure>')

    body = re.sub(r'<p><img alt="([^"]*)" src="([^"]+)"\s*/?></p>', figure, body)

    # Links between chapters.
    def link(m):
        target, anchor = m.group(1), m.group(2) or ""
        if print_mode:
            ref = f"ch{int(target[:2]):02d}-" + (anchor[1:] if anchor else
                                                   slug(CHAPTER_TITLES.get(target, target)))
            return f'href="#{ref}"'
        return f'href="{target}.html{anchor}"'

    body = re.sub(r'href="(\d\d-[\w-]+)\.md(#[\w-]+)?"', link, body)
    # Other repository files (PLUGINS.md) point to GitHub.
    body = re.sub(r'href="\.\./([\w.-]+\.md)"',
                  r'href="https://github.com/qjiangzhao/RepeatViewer/blob/main/docs/\1"', body)
    # Blockquotes that start with a bold word are notes.
    body = body.replace("<blockquote>", '<blockquote class="note">')
    ch.html = body
    return body


CHAPTER_TITLES: dict[str, str] = {}


def sidebar(chapters: list[Chapter], current: Chapter | None) -> str:
    items = []
    for ch in chapters:
        cls = ' class="current"' if ch is current else ""
        sub = ""
        if ch is current and ch.sections:
            sub = "<ul>" + "".join(
                f'<li><a href="#{sid}">{num} {html.escape(t)}</a></li>'
                for num, sid, t in ch.sections) + "</ul>"
        items.append(f'<li{cls}><a href="{ch.page}"><span class="num">{ch.num}</span>'
                     f'{html.escape(ch.title)}</a>{sub}</li>')
    return "<ol class=\"toc\">" + "".join(items) + "</ol>"


def page(title: str, body: str, nav: str, *, prev=None, nxt=None, ver: str) -> str:
    pn = '<nav class="pager">'
    pn += (f'<a class="prev" href="{prev.page}"><small>Previous</small>{html.escape(prev.title)}</a>'
           if prev else "<span></span>")
    pn += (f'<a class="next" href="{nxt.page}"><small>Next</small>{html.escape(nxt.title)}</a>'
           if nxt else "<span></span>")
    pn += "</nav>"
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(title)} · {TITLE}</title>
<link rel="stylesheet" href="manual.css">
</head><body>
<header class="topbar"><a class="brand" href="index.html">RepeatViewer <span>User Manual {ver}</span></a>
<input id="search" type="search" placeholder="Search the manual" autocomplete="off">
<a class="pdf" href="{PDF_NAME}">PDF</a></header>
<div class="layout"><aside class="sidebar">{nav}</aside>
<main id="content"><div id="results" hidden></div>{body}{pn}</main></div>
<script src="search_index.js"></script><script src="manual.js"></script>
</body></html>"""


CSS = r"""
:root{--ink:#1d1d1f;--muted:#5f6368;--line:#dcdde1;--bg:#fff;--side:#f6f7f9;--accent:#1f5fbf;
--code:#f3f4f6;--note:#eef4fc;--font:-apple-system,BlinkMacSystemFont,"Inter","Segoe UI","Helvetica Neue",Arial,sans-serif;
--mono:"SF Mono",Menlo,Consolas,"DejaVu Sans Mono",monospace}
@media (prefers-color-scheme:dark){:root:not([data-theme=light]){--ink:#e8e8ea;--muted:#a0a3a8;--line:#3a3b3f;
--bg:#1b1c1e;--side:#232427;--accent:#7fb0ff;--code:#2a2b2f;--note:#22303f}}
*{box-sizing:border-box}html{scroll-padding-top:64px}
body{margin:0;background:var(--bg);color:var(--ink);font:15px/1.6 var(--font)}
a{color:var(--accent);text-decoration:none}a:hover{text-decoration:underline}
.topbar{position:sticky;top:0;z-index:5;display:flex;gap:16px;align-items:center;height:52px;padding:0 20px;
background:var(--bg);border-bottom:1px solid var(--line)}
.brand{font-weight:650;color:var(--ink)}.brand span{font-weight:400;color:var(--muted);margin-left:6px}
#search{margin-left:auto;width:min(320px,40vw);padding:6px 10px;border:1px solid var(--line);border-radius:8px;
background:var(--bg);color:var(--ink);font:inherit}
.pdf{border:1px solid var(--line);border-radius:8px;padding:4px 12px;color:var(--ink)}
.layout{display:flex;max-width:1280px;margin:0 auto}
.sidebar{position:sticky;top:52px;align-self:flex-start;width:280px;height:calc(100vh - 52px);overflow:auto;
padding:20px 12px 40px 20px;background:var(--side);border-right:1px solid var(--line);font-size:14px}
.toc{list-style:none;margin:0;padding:0}.toc>li{margin:2px 0}
.toc a{display:block;padding:4px 8px;border-radius:6px;color:var(--ink)}
.toc .num{display:inline-block;width:24px;color:var(--muted)}
.toc li.current>a{background:var(--bg);font-weight:600;box-shadow:0 0 0 1px var(--line)}
.toc ul{list-style:none;margin:4px 0 8px 24px;padding:0}.toc ul a{color:var(--muted);padding:2px 8px}
main{flex:1;min-width:0;padding:32px 48px 64px;max-width:900px}
h1{font-size:30px;line-height:1.2;margin:0 0 20px}.chapter-label{display:block;font-size:13px;font-weight:600;
color:var(--accent);text-transform:uppercase;letter-spacing:.06em;margin-bottom:6px}
h2{font-size:22px;margin:40px 0 12px;padding-top:12px;border-top:1px solid var(--line)}
.secnum{color:var(--muted);font-weight:500;margin-right:10px}
h3{font-size:17px;margin:28px 0 8px}h4{font-size:15px;margin:20px 0 6px}
p,ul,ol{margin:0 0 12px}li{margin:3px 0}
table{border-collapse:collapse;width:100%;margin:12px 0 18px;font-size:14px}
th,td{border:1px solid var(--line);padding:6px 9px;text-align:left;vertical-align:top}
th{background:var(--side);font-weight:600}
code{font:13px var(--mono);background:var(--code);padding:1px 4px;border-radius:4px}
pre{background:var(--code);padding:12px 14px;border-radius:8px;overflow:auto}pre code{padding:0;background:none}
figure{margin:18px 0 24px;text-align:center}
figure img{max-width:100%;height:auto;border:1px solid var(--line);border-radius:8px;
box-shadow:0 2px 10px rgba(0,0,0,.08)}
figcaption{font-size:13px;color:var(--muted);margin-top:8px;text-align:left}
blockquote.note{margin:16px 0;padding:10px 14px;background:var(--note);border-left:3px solid var(--accent);
border-radius:6px}blockquote.note p{margin:0}
.pager{display:flex;justify-content:space-between;gap:16px;margin-top:48px;padding-top:20px;
border-top:1px solid var(--line)}.pager a{display:block;max-width:48%}.pager small{display:block;color:var(--muted)}
.pager .next{text-align:right;margin-left:auto}
#results{margin-bottom:28px;padding:12px 16px;border:1px solid var(--line);border-radius:8px;background:var(--side)}
#results h4{margin:0 0 6px}#results ol{margin:0;padding-left:20px}#results p{margin:2px 0 8px;color:var(--muted);font-size:13px}
.cover{padding:40px 0 10px}.cover h1{font-size:40px}.cover p.lead{font-size:18px;color:var(--muted)}.cover p.authors{font-size:16px;font-weight:500;margin:-8px 0 12px}.cover .affil{font-size:14px;font-weight:400;color:var(--muted)}
.chapters{display:grid;grid-template-columns:repeat(auto-fill,minmax(250px,1fr));gap:12px;margin-top:24px}
.chapters a{display:block;padding:14px 16px;border:1px solid var(--line);border-radius:10px;color:var(--ink)}
.chapters a:hover{border-color:var(--accent);text-decoration:none}.chapters small{display:block;color:var(--muted)}
@media (max-width:860px){.sidebar{display:none}main{padding:24px 16px}#search{width:45vw}}
"""

JS = r"""
(function(){
  var box=document.getElementById('search'),out=document.getElementById('results');
  if(!box||!window.MANUAL_INDEX)return;
  function esc(s){return s.replace(/[&<>]/g,function(c){return{'&':'&amp;','<':'&lt;','>':'&gt;'}[c]})}
  box.addEventListener('input',function(){
    var q=box.value.trim().toLowerCase();
    if(q.length<2){out.hidden=true;out.innerHTML='';return}
    var words=q.split(/\s+/),hits=[];
    MANUAL_INDEX.forEach(function(s){
      var t=(s.t+' '+s.x).toLowerCase(),ok=words.every(function(w){return t.indexOf(w)>=0});
      if(ok){var i=s.x.toLowerCase().indexOf(words[0]);
        hits.push({s:s,snip:i<0?s.x.slice(0,140):s.x.slice(Math.max(0,i-60),i+100)})}
    });
    out.hidden=false;
    out.innerHTML='<h4>'+hits.length+' section'+(hits.length==1?'':'s')+' found</h4><ol>'+hits.slice(0,30).map(function(h){
      return '<li><a href="'+h.s.u+'">'+esc(h.s.t)+'</a><p>…'+esc(h.snip)+'…</p></li>'}).join('')+'</ol>';
  });
  document.addEventListener('keydown',function(e){if(e.key=='/'&&document.activeElement!==box){e.preventDefault();box.focus()}});
})();
"""

PRINT_CSS = r"""
@page{size:A4;margin:20mm 18mm 20mm 18mm}
body{font:10pt/1.5 -apple-system,"Inter","Helvetica Neue",Arial,sans-serif;color:#1d1d1f;margin:0}
a{color:#1f5fbf;text-decoration:none}
.coverpage{height:250mm;display:flex;flex-direction:column;justify-content:center;page-break-after:always}
.coverpage h1{font-size:34pt;margin:0 0 8pt}.coverpage .sub{font-size:14pt;color:#5f6368}
.coverpage .authors{font-size:12pt;font-weight:500;margin-top:14pt}
.coverpage .affil{font-size:10pt;color:#5f6368;margin-top:3pt}
.coverpage img{width:100%;border:1px solid #dcdde1;border-radius:6pt;margin:24pt 0}
.coverpage .meta{font-size:10pt;color:#5f6368}
.tocpage{page-break-after:always}.tocpage h1{font-size:22pt}
.tocpage ol{list-style:none;padding:0}.tocpage>ol>li{margin:8pt 0 2pt;font-weight:600}
.tocpage ol ol{margin:2pt 0 0 18pt;font-weight:400}.tocpage ol ol li{margin:1pt 0}
.tocpage .n{display:inline-block;min-width:28pt;color:#5f6368}
section.chapter{page-break-before:always}
h1{font-size:22pt;margin:0 0 14pt;line-height:1.2}.chapter-label{display:block;font-size:9pt;color:#1f5fbf;
text-transform:uppercase;letter-spacing:.08em;margin-bottom:4pt}
h2{font-size:14pt;margin:20pt 0 6pt;page-break-after:avoid}.secnum{color:#5f6368;margin-right:8pt}
h3{font-size:11.5pt;margin:14pt 0 4pt;page-break-after:avoid}h4{font-size:10pt;page-break-after:avoid}
p,ul,ol{margin:0 0 7pt}li{margin:1.5pt 0}
table{border-collapse:collapse;width:100%;margin:6pt 0 10pt;font-size:8.6pt;page-break-inside:auto}
tr{page-break-inside:avoid}th,td{border:0.5pt solid #c9cbd0;padding:3pt 5pt;text-align:left;vertical-align:top}
th{background:#f3f4f6}
code{font:8.6pt "SF Mono",Menlo,"DejaVu Sans Mono",monospace;background:#f3f4f6;padding:0 2pt;border-radius:2pt}
pre{background:#f3f4f6;padding:6pt 8pt;border-radius:4pt;white-space:pre-wrap;page-break-inside:avoid}
pre code{background:none}
figure{margin:8pt 0 12pt;text-align:center;page-break-inside:avoid}
figure img{max-width:100%;max-height:150mm;height:auto;border:0.5pt solid #c9cbd0;border-radius:4pt}
figcaption{font-size:8.6pt;color:#5f6368;text-align:left;margin-top:4pt}
blockquote.note{margin:8pt 0;padding:6pt 10pt;background:#eef4fc;border-left:2pt solid #1f5fbf}
blockquote.note p{margin:0}
"""


def build(pdf: bool) -> None:
    ver = version()
    chapters = [Chapter(p) for p in sorted(SRC.glob("[0-9][0-9]-*.md"))]
    if not chapters:
        sys.exit(f"No chapters in {SRC}")
    for ch in chapters:
        CHAPTER_TITLES[ch.stem] = ch.title

    # Web pages.
    index = []
    for i, ch in enumerate(chapters):
        body = render(ch, print_mode=False)
        prev = chapters[i - 1] if i else None
        nxt = chapters[i + 1] if i + 1 < len(chapters) else None
        (SRC / ch.page).write_text(page(ch.title, body, sidebar(chapters, ch), prev=prev,
                                        nxt=nxt, ver=ver), encoding="utf-8")
        # Search index: one entry per section.
        parts = re.split(r'(<h2 id="[^"]+">.*?</h2>)', body)
        current_id, current_t = "", ch.title
        for part in parts:
            m = re.match(r'<h2 id="([^"]+)">(.*?)</h2>', part)
            if m:
                current_id = m.group(1)
                current_t = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", m.group(2))).strip()
                continue
            text = re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", " ", part))).strip()
            if text:
                index.append({"t": f"{ch.title}: {current_t}" if current_id else ch.title,
                              "u": f"{ch.page}#{current_id}" if current_id else ch.page,
                              "x": text[:4000]})

    cards = "".join(
        f'<a href="{ch.page}"><small>Chapter {ch.num}</small>{html.escape(ch.title)}</a>'
        for ch in chapters)
    cover = (f'<div class="cover"><h1>RepeatViewer User Manual</h1>'
             f'<p class="authors">{AUTHORS}<br><span class="affil">{AFFILIATION}</span></p>'
             f'<p class="lead">Version {ver}. A free desktop application for repetitive DNA: '
             f'transposable elements, satellites and the work around them.</p>'
             f'<figure><img src="img/main_window.png" alt="RepeatViewer" width="1000"></figure>'
             f'<p>Read it online, or download the <a href="{PDF_NAME}">PDF</a>. '
             f'Press <code>/</code> to search.</p>'
             f'<div class="chapters">{cards}</div></div>')
    (SRC / "index.html").write_text(page("Contents", cover, sidebar(chapters, None),
                                         nxt=chapters[0], ver=ver), encoding="utf-8")
    (SRC / "manual.css").write_text(CSS, encoding="utf-8")
    (SRC / "manual.js").write_text(JS, encoding="utf-8")
    (SRC / "search_index.js").write_text(
        "window.MANUAL_INDEX=" + json.dumps(index, ensure_ascii=False) + ";", encoding="utf-8")
    print(f"Web manual: {len(chapters)} chapters, {len(index)} search entries")

    # README for browsing the chapters on GitHub.
    lines = [f"# {TITLE}", "",
             "Read it as web pages (`index.html` in this folder, also under Help in the app) "
             f"or as a [PDF]({PDF_NAME}). The chapters below are the sources.", ""]
    lines += [f"{ch.num}. [{ch.title}]({ch.stem}.md)" for ch in chapters]
    lines += ["", "To change the manual, edit the chapter files, then rebuild:", "",
              "```", "python docs/manual/tools/build_manual.py --pdf", "```", "",
              "After a change to a window, take its screenshot again:", "",
              "```", "python docs/manual/tools/take_screenshots.py --only NAME", "```", ""]
    (SRC / "README.md").write_text("\n".join(lines), encoding="utf-8")

    if pdf:
        build_pdf(chapters, ver)

    # In-app copy (Help → RepeatViewer Manual).
    if APP_OUT.exists():
        shutil.rmtree(APP_OUT)
    (APP_OUT / "img").mkdir(parents=True)
    used = set()
    for f in SRC.glob("*.html"):
        used.update(re.findall(r'src="img/([^"]+)"', f.read_text(encoding="utf-8")))
    for name in ["index.html", "manual.css", "manual.js", "search_index.js",
                 *[ch.page for ch in chapters]]:
        shutil.copy2(SRC / name, APP_OUT / name)
    for name in sorted(used):
        if (IMG / name).exists():
            shutil.copy2(IMG / name, APP_OUT / "img" / name)
        else:
            print(f"  missing image: img/{name}")
    if (SRC / PDF_NAME).exists():
        shutil.copy2(SRC / PDF_NAME, APP_OUT / PDF_NAME)
    unused = sorted({p.name for p in IMG.glob("*.png")} - used)
    if unused:
        print("  images not used in the text:", ", ".join(unused))
    print(f"In-app copy: {APP_OUT.relative_to(ROOT)}")


def build_pdf(chapters: list[Chapter], ver: str) -> None:
    try:
        from playwright.sync_api import sync_playwright
    except ImportError:
        print("PDF skipped: pip install playwright && playwright install chromium")
        return
    toc = "<ol>" + "".join(
        f'<li><a href="#ch{ch.num:02d}-{slug(ch.title)}"><span class="n">{ch.num}</span>'
        f'{html.escape(ch.title)}</a><ol>' + "".join(
            f'<li><a href="#ch{ch.num:02d}-{sid}"><span class="n">{num}</span>'
            f'{html.escape(t)}</a></li>' for num, sid, t in ch.sections) + "</ol></li>"
        for ch in chapters) + "</ol>"
    body = "".join(f'<section class="chapter">{render(ch, print_mode=True)}</section>'
                   for ch in chapters)
    doc = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<title>{TITLE}</title><style>{PRINT_CSS}</style></head><body>
<div class="coverpage"><h1>RepeatViewer</h1><div class="sub">User Manual · Version {ver}</div>
<div class="authors">{AUTHORS}</div>
<div class="affil">{AFFILIATION}</div>
<img src="img/main_window.png" alt="">
<div class="meta">{date.today():%B %Y} · github.com/qjiangzhao/RepeatViewer</div></div>
<div class="tocpage"><h1>Contents</h1>{toc}</div>{body}</body></html>"""
    tmp = SRC / "_print.html"
    tmp.write_text(doc, encoding="utf-8")
    out = SRC / PDF_NAME
    footer = ('<div style="width:100%;font-size:7.5pt;color:#888;padding:0 18mm;display:flex;'
              'justify-content:space-between;font-family:Helvetica,Arial,sans-serif">'
              f'<span>RepeatViewer User Manual {ver}</span>'
              '<span><span class="pageNumber"></span> / <span class="totalPages"></span></span></div>')
    with sync_playwright() as p:
        browser = p.chromium.launch()
        pg = browser.new_page()
        pg.goto(tmp.as_uri(), wait_until="networkidle")
        kwargs = dict(path=str(out), format="A4", print_background=True,
                      display_header_footer=True, header_template="<span></span>",
                      footer_template=footer,
                      margin={"top": "18mm", "bottom": "18mm", "left": "18mm", "right": "18mm"})
        try:
            pg.pdf(outline=True, tagged=True, **kwargs)
        except TypeError:  # older Playwright
            pg.pdf(**kwargs)
        browser.close()
    tmp.unlink()
    print(f"PDF: {out.relative_to(ROOT)} ({out.stat().st_size // 1024} KB)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Build the RepeatViewer user manual.")
    ap.add_argument("--pdf", action="store_true", help="also build the PDF")
    build(ap.parse_args().pdf)
