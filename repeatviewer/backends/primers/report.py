"""HTML report for a primer result (same look as the Mapping report)."""

from __future__ import annotations

from html import escape as e
from pathlib import Path

from repeatviewer.backends.primers.gel import Lane, gel_svg

CITES = (
    "Untergasser A. et al. (2012) Primer3--new capabilities and interfaces. "
    "Nucleic Acids Res 40:e115.",
    "Koressaar T., Remm M. (2007) Enhancements and modifications of primer design program "
    "Primer3. Bioinformatics 23:1289–1291.",
    "SantaLucia J. (1998) A unified view of polymer, dumbbell, and oligonucleotide DNA "
    "nearest-neighbor thermodynamics. PNAS 95:1460–1465.",
)

CSS = """
:root{--fg:#1d1d1b;--muted:#6b6b67;--line:#e4e4e0;--soft:#f6f6f4}
body{font-family:-apple-system,"Segoe UI",Helvetica,Arial,sans-serif;margin:0;color:var(--fg);background:#fff;font-size:13px}
main{max-width:940px;margin:0 auto;padding:28px 24px 40px}
h1{font-size:19px;font-weight:600;margin:0 0 4px} .sub{color:var(--muted);margin:0 0 20px}
h2{font-size:14px;font-weight:600;margin:28px 0 8px;padding-bottom:6px;border-bottom:1px solid var(--line)}
.tiles{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:12px}
.tile{background:var(--soft);border-radius:8px;padding:12px 14px}
.tile .k{color:var(--muted);font-size:12px} .tile .v{font-size:21px;font-weight:600;margin:2px 0}
.tile .s{color:var(--muted);font-size:12px}
table{border-collapse:collapse;width:100%} th,td{padding:5px 8px;text-align:left;vertical-align:top}
.grid th{color:var(--muted);font-weight:500;border-bottom:1px solid var(--line)}
.grid td{border-bottom:1px solid #f0f0ec} td.n,th.n{text-align:right;font-variant-numeric:tabular-nums;white-space:nowrap}
.scroll{overflow-x:auto}
td.c{font-family:Menlo,Consolas,monospace;font-size:12px;white-space:nowrap} td.m{color:var(--muted)}
td.c i{font-style:normal;color:var(--muted)}
.details th{width:150px;color:var(--muted);font-weight:500;padding-left:0} .details td{padding-left:0}
.warn{background:#fff6e5;border:1px solid #f0c36d;border-radius:8px;padding:10px 14px;margin-bottom:16px}
.methods{line-height:1.55;background:var(--soft);border-radius:8px;padding:12px 14px}
.note{color:var(--muted);font-size:12px;margin:6px 0 0}
pre{font-size:11px;background:var(--soft);border-radius:8px;padding:10px 12px;white-space:pre;overflow-x:auto;margin:0 0 8px}
.two{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:24px}
.gel{display:flex;gap:20px;align-items:flex-start}
.tag{display:inline-block;border-radius:4px;padding:0 5px;font-size:11px;background:#eef3ee;color:#2e7d32}
.tag.bad{background:#fdeeee;color:#b3261e}
ul{margin:4px 0 0 18px;padding:0} li{margin:2px 0}
@media (max-width:700px){.tiles{grid-template-columns:repeat(2,minmax(0,1fr))}.two{grid-template-columns:1fr}}
"""


def _tile(k: str, v: str, s: str = "") -> str:
    return f'<div class="tile"><div class="k">{e(k)}</div><div class="v">{e(v)}</div>' \
           f'<div class="s">{e(s)}</div></div>'


def _seq_html(seq: str) -> str:
    """Lower-case tail grey, binding part plain."""
    i = 0
    while i < len(seq) and seq[i].islower():
        i += 1
    if 0 < i < len(seq):
        return f"<i>{e(seq[:i])}</i>{e(seq[i:])}"
    return e(seq)


def _num(v, fmt="{:.1f}") -> str:
    try:
        return fmt.format(float(v))
    except (TypeError, ValueError):
        return "–"


def _design_sections(d: dict) -> tuple[list[str], str]:
    pairs = d.get("pairs", [])
    tiles = []
    best = pairs[0] if pairs else None
    tiles.append(_tile("Pairs found" if any(p.get("forward") for p in pairs) or not pairs
                       else "Probes found", f"{len(pairs)}", d.get("mode_label", "")))
    if best:
        f, r, pr = best.get("forward"), best.get("reverse"), best.get("probe")
        if f and r:
            tiles.append(_tile("Best pair Tm", f"{f['tm']:.1f} / {r['tm']:.1f}",
                               f"Δ {best.get('tm_diff', 0):.1f} °C"))
            tiles.append(_tile("Best product", f"{best.get('product_size', 0):,} bp",
                               best.get("name", "")))
        elif pr:
            tiles.append(_tile("Best probe Tm", f"{pr['tm']:.1f} °C", f"{pr['length']} nt"))
            cov = best.get("coverage") or {}
            if cov:
                tiles.append(_tile("Copies covered", f"{cov.get('target_amplified', 0)}/"
                                   f"{cov.get('target_total', 0)}", "with the mismatches allowed"))
    genome = d.get("genome") or {}
    if genome and best:
        sp = best.get("specificity") or {}
        tiles.append(_tile("Genome products", f"{sp.get('products', 0):,}",
                           Path(genome.get("genome", "")).name))
    elif best and best.get("coverage") and best.get("forward"):
        cov = best["coverage"]
        tiles.append(_tile("Sequences amplified", f"{cov.get('target_amplified', 0)}/"
                           f"{cov.get('target_total', 0)}", "best pair"))
    while len(tiles) < 4:
        tiles.append(_tile("Templates", f"{len(d.get('templates', []))}", ""))
        break
    rows = []
    for p in pairs:
        notes = "; ".join(p.get("notes", []))
        for key, label in (("forward", "F"), ("reverse", "R"), ("probe", "P")):
            o = p.get(key)
            if not o:
                continue
            first = key == next(k for k in ("forward", "reverse", "probe") if p.get(k))
            rows.append(
                "<tr>"
                + (f"<td rowspan={sum(1 for k in ('forward', 'reverse', 'probe') if p.get(k))}>"
                   f"<b>{e(p['name'])}</b>"
                   + (f" <span class=tag>{e(p['group'])}</span>" if p.get("group") else "")
                   + f"<div class=note>{e(notes)}</div></td>" if first else "")
                + f"<td>{label}</td><td class=c>{_seq_html(o['sequence'])}</td>"
                f"<td class=n>{o.get('start', 0):,}–{o.get('end', 0):,}</td>"
                f"<td class=n>{_num(o['tm'])}</td><td class=n>{_num(o['gc'], '{:.0f}')}</td>"
                f"<td class=n>{_num(o.get('hairpin_dg'))}</td>"
                f"<td class=n>{_num(o.get('self_dimer_dg'))}</td>"
                + (f"<td class=n rowspan={sum(1 for k in ('forward', 'reverse', 'probe') if p.get(k))}>"
                   f"{p.get('product_size', 0) or '–'}</td>" if first else "")
                + "</tr>")
        te = (p.get("extra") or {}).get("te_primer")
        if te:
            rows.append(f"<tr><td></td><td>TE</td><td class=c>{_seq_html(te['sequence'])}</td>"
                        f"<td class=n>{te.get('start', 0):,}–{te.get('end', 0):,}</td>"
                        f"<td class=n>{_num(te['tm'])}</td><td class=n>{_num(te['gc'], '{:.0f}')}</td>"
                        f"<td class=n>{_num(te.get('hairpin_dg'))}</td>"
                        f"<td class=n>{_num(te.get('self_dimer_dg'))}</td><td></td></tr>")
    body = []
    if rows:
        body.append(
            "<h2>Primers</h2><div class=scroll><table class=grid><tr><th>Pair</th><th></th><th>Sequence 5′→3′</th>"
            "<th class=n>Position</th><th class=n>Tm °C</th><th class=n>GC %</th>"
            "<th class=n>Hairpin ΔG</th><th class=n>Self-dimer ΔG</th><th class=n>Product</th>"
            f"</tr>{''.join(rows)}</table></div>"
            "<p class=note>ΔG in kcal/mol at the reaction temperature; grey lower-case letters "
            "are 5′ tails that do not bind the template.</p>")
    dimers = []
    for p in pairs[:5]:
        lines = p.get("dimer_lines") or []
        if lines:
            dimers.append(f"<b>{e(p['name'])}</b> F × R, ΔG {_num(p.get('dimer_dg'))}"
                          f"<pre>{e(chr(10).join(lines))}</pre>")
    if dimers:
        body.append("<h2>Strongest primer-dimers</h2>" + "".join(dimers))
    return tiles, "".join(body)


def _mm_tag(n: int) -> str:
    return '<span class="tag">exact</span>' if n == 0 else \
        f'<span class="tag bad">{n} mm</span>'


def _tested_tag(*values) -> str:
    vals = [v for v in values if v]
    if "failed" in vals:
        return ' <span class="tag bad">failed before</span>'
    if "worked" in vals:
        return ' <span class="tag">worked before</span>'
    return ""


def _stock_section(st: dict) -> str:
    """Map on saved primers: what the lab already has for these sequences."""
    ready, single = st.get("ready", []), st.get("single", [])
    html = "<h2>From your Primer Database</h2>"
    if not ready and not single:
        return html + ("<p>None of the {n} saved primers binds these sequences with the "
                       "mismatches allowed. Design new primers (Tools › Primers and Probes › "
                       "Design), or allow more mismatches.</p>").format(n=st.get("primers", 0))
    if ready:
        rows = "".join(
            f"<tr><td><b>{e(r['forward'])}</b> + <b>{e(r['reverse'])}</b>"
            + (' <span class="tag">saved pair</span>' if r.get("saved_pair") else "")
            + _tested_tag(r.get("forward_tested"), r.get("reverse_tested"))
            + f"</td><td class=n>{', '.join(f'{x:,}' for x in r['sizes'])}"
            + (" …" if r["products"] > len(r["sizes"]) else "")
            + f"</td><td>{_mm_tag(r['mm_forward'])} {_mm_tag(r['mm_reverse'])}</td>"
            f"<td class=n>{r['records']} / {st.get('sequences', 0)}</td>"
            f"<td class=m>{e(r.get('forward_where') or '?')}"
            + (f"<br>{e(r['reverse_where'])}" if r.get("reverse_where")
               and r.get("reverse_where") != r.get("forward_where") else "")
            + "</td></tr>" for r in ready[:60])
        html += ("<p>Pairs of saved primers that give a product. Saved pairs and exact "
                 "matches first; mismatches (mm) in a primer lower the chance that it works, "
                 "most of all near its 3′ end.</p>"
                 "<div class=scroll><table class=grid><tr><th>Primers</th><th class=n>"
                 "Product bp</th><th>Forward / reverse</th><th class=n>Sequences</th>"
                 f"<th>Where</th></tr>{rows}</table></div>")
        if st.get("ready_total", 0) > 60:
            html += f"<p class=note>First 60 of {st['ready_total']} shown.</p>"
    if single:
        rows = "".join(
            f"<tr><td><b>{e(x['name'])}</b></td><td>{_mm_tag(x['mismatches'])}</td>"
            f"<td>{e(x['record'])} {e(x['strand'])} {x['start']:,}</td>"
            f"<td class=n>{x['sites']}</td><td class=m>{e(x.get('where') or '?')}</td></tr>"
            for x in single[:40])
        html += ("<h2>Single primers that bind</h2><p>These bind, but no saved partner gives a "
                 "product here. One new partner primer next to the site is enough to use "
                 "them.</p><div class=scroll><table "
                 "class=grid><tr><th>Primer</th><th>Match</th><th>Best site</th><th class=n>"
                 f"Sites</th><th>Where</th></tr>{rows}</table></div>")
    return html


def _map_sections(d: dict) -> tuple[list[str], str]:
    c = d.get("counts", {})
    stock = d.get("stock")
    if stock is not None:
        tiles = [
            _tile("Saved primers tested", f"{stock.get('primers', 0)}", ""),
            _tile("Ready pairs", f"{stock.get('ready_total', 0)}",
                  f"{sum(1 for r in stock.get('ready', []) if r.get('saved_pair'))} saved pairs"),
            _tile("Single primers", f"{stock.get('single_total', 0)}", "bind, no partner"),
            _tile("Sequences", f"{c.get('sequences', 0)}",
                  f"{c.get('sequences_with_product', 0)} with a product"),
        ]
        body = _stock_section(stock) + _map_sections_body(d)
        return tiles, body
    tiles = [
        _tile("Primers tested", f"{c.get('primers', 0)}", ""),
        _tile("Binding sites", f"{c.get('sites', 0):,}", f"in {c.get('sequences', 0)} sequences"),
        _tile("Products", f"{c.get('products', 0):,}",
              f"in {c.get('sequences_with_product', 0)} sequences"),
    ]
    genome = d.get("genome") or {}
    if genome:
        tiles.append(_tile("Genome products", f"{genome.get('products', 0):,}",
                           Path(genome.get("genome", "")).name))
    else:
        tiles.append(_tile("Sequences", f"{c.get('sequences', 0)}", ""))
    return tiles, _map_sections_body(d)


def _map_sections_body(d: dict) -> str:
    stock = d.get("stock")
    if stock is not None:
        # Only the primers that bind are worth listing.
        bound = {s["primer"] for s in d.get("sites", [])}
        d = dict(d, primers=[p for p in d.get("primers", []) if p["name"] in bound])
    prows = "".join(
        f"<tr><td><b>{e(p['name'])}</b></td><td>{e(p.get('kind', ''))}</td>"
        f"<td class=c>{_seq_html(p['sequence'])}</td><td class=n>{p.get('length', 0)}</td>"
        f"<td class=n>{_num(p.get('tm'))}</td><td class=n>{_num(p.get('gc'), '{:.0f}')}</td>"
        f"<td class=m>{e('; '.join(p.get('warnings', [])))}</td></tr>"
        for p in d.get("primers", []))
    srows = "".join(
        f"<tr><td>{e(s['primer'])}</td><td>{e(s['record'])}</td><td>{e(s['strand'])}</td>"
        f"<td class=n>{s['start']:,}–{s['end']:,}</td><td class=n>{s['mismatches']}</td>"
        f"<td class=n>{_num(s.get('tm')) if s.get('tm') else '–'}</td></tr>"
        for s in d.get("sites", [])[:60])
    rrows = "".join(
        f"<tr><td>{e(p['record'])}</td><td>{e(p['forward'])} + {e(p['reverse'])}</td>"
        f"<td class=n>{p['start']:,}–{p['end']:,}</td><td class=n>{p['size']:,}</td>"
        f"<td class=n>{p['mismatches']}</td></tr>"
        for p in d.get("products", [])[:60])
    body = f"<h2>Primers</h2><div class=scroll><table class=grid><tr><th>Name</th><th>Kind</th><th>Sequence 5′→3′</th>" \
           f"<th class=n>Length</th><th class=n>Tm °C</th><th class=n>GC %</th><th>Check</th>" \
           f"</tr>{prows}</table></div>" if prows else ""
    if rrows:
        body += "<h2>Products</h2><table class=grid><tr><th>Sequence</th><th>Primers</th>" \
                "<th class=n>Position</th><th class=n>Size bp</th><th class=n>Mismatches</th>" \
                f"</tr>{rrows}</table>"
        if len(d.get("products", [])) > 60:
            body += "<p class=note>First 60 shown; all are in the products table.</p>"
    if srows:
        body += "<h2>Binding sites</h2><table class=grid><tr><th>Primer</th><th>Sequence</th>" \
                "<th>Strand</th><th class=n>Position</th><th class=n>Mismatches</th>" \
                f"<th class=n>Binding Tm</th></tr>{srows}</table>"
        if len(d.get("sites", [])) > 60:
            body += "<p class=note>First 60 shown; all are in the binding-sites table.</p>"
    return body


def methods_text(d: dict) -> str:
    s = d.get("settings", {})
    if d.get("kind") == "map":
        return (
            "Primers were tested against the sequences with RepeatViewer's primer search "
            f"(up to {s.get('max_mismatches', 0)} mismatches, the last "
            f"{s.get('three_prime_exact', 0)} bases at the 3′ end matching exactly, both "
            "strands). Products were predicted between facing binding sites"
            + (f" of {s.get('product_min')}–{s.get('product_max')} bp" if s.get("use_product_range") else "")
            + f". Tm was calculated with Primer3 ({d.get('thermo', '')})."
        )
    text = (f"Primers were designed with {d.get('engine', 'Primer3')} in RepeatViewer "
            f"({d.get('mode_label', '')}): {s.get('size_min')}–{s.get('size_max')} nt, Tm "
            f"{s.get('tm_min')}–{s.get('tm_max')} °C (optimum {s.get('tm_opt')} °C), GC "
            f"{s.get('gc_min')}–{s.get('gc_max')} %, products {s.get('product_min')}–"
            f"{s.get('product_max')} bp; {d.get('thermo', '')}.")
    return text


def write_report(path: Path, d: dict) -> Path:
    title = d.get("title", "Primers")
    if d.get("kind") == "map":
        tiles, body = _map_sections(d)
    else:
        tiles, body = _design_sections(d)
    warn = "".join(f"<li>{e(w)}</li>" for w in d.get("warnings", []))
    details = [("Created", d.get("created", "")), ("Inputs", ", ".join(Path(p).name for p in d.get("inputs", []))),
               ("Engine", d.get("engine", "")), ("Conditions", d.get("thermo", ""))]
    if d.get("mode_label"):
        details.insert(1, ("Mode", d["mode_label"]))
    detail_html = "".join(f"<tr><th>{e(k)}</th><td>{e(str(v))}</td></tr>" for k, v in details if v)
    gel = d.get("gel") or {}
    gel_html = ""
    if gel.get("lanes"):
        svg = gel_svg([Lane.from_dict(x) for x in gel["lanes"]], ladder=gel.get("ladder", ""))
        gel_html = (f"<h2>Virtual gel</h2><div class=gel>{svg}<div class=note>"
                    f"Ladder: {e(gel.get('ladder', ''))}. Band positions follow log(size); "
                    "green = the expected product, brighter = more copies of that size."
                    "</div></div>")
    explain = "".join(f"<li>{e(x)}</li>" for x in d.get("explain", []))
    files = "".join(f"<li>{e(name)}</li>" for name in (d.get("files") or {}))
    cites = "".join(f"<li>{e(c)}</li>" for c in CITES)
    html = f"""<!doctype html><html><head><meta charset="utf-8"><title>{e(title)}</title>
<style>{CSS}</style></head><body><main>
<h1>{e(title)}</h1><p class="sub">{e(d.get('created', ''))}</p>
{f'<div class="warn"><b>Check before use</b><ul>{warn}</ul></div>' if warn else ''}
<div class="tiles">{''.join(tiles)}</div>
<h2>Run</h2><table class="details">{detail_html}</table>
{body}
{gel_html}
{f'<h2>Why</h2><p class=note>What Primer3 rejected, per template.</p><ul>{explain}</ul>' if explain else ''}
<h2>Methods text</h2><div class="methods">{e(methods_text(d))}</div>
<div class="two"><div><h2>Cite</h2><ul>{cites}</ul></div><div><h2>Files written</h2><ul>{files}</ul></div></div>
</main></body></html>"""
    Path(path).write_text(html, encoding="utf-8")
    return Path(path)
