"""Tool paper citations shown on Options tabs (Zotero RepeatViewer → Tools)."""

from __future__ import annotations


def _entry(tool: str, authors: str, year: str, title: str, journal: str, doi: str) -> str:
    href = f"https://doi.org/{doi}"
    return (
        f"<p style='margin:0 0 10px 0'>"
        f"<b>{tool}</b><br>"
        f"{authors} ({year}) {title}. "
        f"<i>{journal}</i>. "
        f'<a href="{href}">doi:{doi}</a>'
        f"</p>"
    )


# BLAST (local BLAST+ / Online BLAST)
CITATION_BLAST = _entry(
    "BLAST+",
    "Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL",
    "2009",
    "BLAST+: architecture and applications",
    "BMC Bioinformatics 10:421",
    "10.1186/1471-2105-10-421",
)

# Align → MAFFT / MUSCLE / FastGA
CITATION_ALIGN = (
    _entry(
        "MAFFT",
        "Katoh K, Standley DM",
        "2013",
        "MAFFT multiple sequence alignment software version 7: improvements in "
        "performance and usability",
        "Molecular Biology and Evolution 30:772–780",
        "10.1093/molbev/mst010",
    )
    + _entry(
        "MUSCLE",
        "Edgar RC",
        "2004",
        "MUSCLE: multiple sequence alignment with high accuracy and high throughput",
        "Nucleic Acids Research 32:1792–1797",
        "10.1093/nar/gkh340",
    )
    + _entry(
        "FastGA",
        "Myers G, Durbin R, Zhou C",
        "2025",
        "FastGA: fast genome alignment",
        "Bioinformatics Advances 5:vbaf238",
        "10.1093/bioadv/vbaf238",
    )
)

# Dotplot → FlexiDot
CITATION_DOTPLOT = _entry(
    "FlexiDot",
    "Seibt KM, Schmidt T, Heitkam T",
    "2018",
    "FlexiDot: highly customizable, ambiguity-aware dotplots for visual sequence analyses",
    "Bioinformatics 34:3575–3577",
    "10.1093/bioinformatics/bty395",
)

# TE-Aid → beginner curation guide
CITATION_TEAID = _entry(
    "TE curation",
    "Goubert C, Craig RJ, Bilat AF, Peona V, Vogan AA, Protasio AV",
    "2022",
    "A beginner’s guide to manual curation of transposable elements",
    "Mobile DNA 13:7",
    "10.1186/s13100-022-00267-1",
)

# Tree → IQ-TREE / FastME / FastTree
CITATION_TREE = (
    _entry(
        "IQ-TREE 2",
        "Minh BQ, Schmidt HA, Chernomor O, Schrempf D, Woodhams MD, von Haeseler A, "
        "Lanfear R",
        "2020",
        "IQ-TREE 2: new models and efficient methods for phylogenetic inference in "
        "the genomic era",
        "Molecular Biology and Evolution 37:1530–1534",
        "10.1093/molbev/msaa015",
    )
    + _entry(
        "FastME 2.0",
        "Lefort V, Desper R, Gascuel O",
        "2015",
        "FastME 2.0: a comprehensive, accurate, and fast distance-based phylogeny "
        "inference program",
        "Molecular Biology and Evolution 32:2798–2800",
        "10.1093/molbev/msv150",
    )
    + _entry(
        "FastTree",
        "Price MN, Dehal PS, Arkin AP",
        "2009",
        "FastTree: computing large minimum evolution trees with profiles instead of "
        "a distance matrix",
        "Molecular Biology and Evolution 26:1641–1650",
        "10.1093/molbev/msp077",
    )
)

# Tools → samtools (pysam / Win-NGS samtools)
CITATION_SAMTOOLS = _entry(
    "samtools / HTSlib",
    "Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, "
    "Keane T, McCarthy SA, Davies RM, Li H",
    "2021",
    "Twelve years of SAMtools and BCFtools",
    "GigaScience 10:giab008",
    "10.1093/gigascience/giab008",
)

# Mapping → Minimap2 / HISAT2
CITATION_MAPPING = (
    _entry(
        "Minimap2",
        "Li H",
        "2018",
        "Minimap2: pairwise alignment for nucleotide sequences",
        "Bioinformatics 34:3094–3100",
        "10.1093/bioinformatics/bty191",
    )
    + _entry(
        "HISAT",
        "Kim D, Langmead B, Salzberg SL",
        "2015",
        "HISAT: a fast spliced aligner with low memory requirements",
        "Nature Methods 12:357–360",
        "10.1038/nmeth.3317",
    )
)

# Sequence Annotation → ORF / tandem / TEtrimmer
CITATION_SEQUENCE_ANNOTATION = (
    _entry(
        "orfipy (ORF)",
        "Singh U, Wurtele ES",
        "2021",
        "orfipy: a fast and flexible tool for extracting ORFs",
        "Bioinformatics 37:3019–3020",
        "10.1093/bioinformatics/btab090",
    )
    + _entry(
        "Pytrf (tandem repeats)",
        "Du L, Sun D, Chen J, Zhou X, Zhao K, Zeng Q, et al.",
        "2025",
        "Pytrf: a python package for finding tandem repeats from genomic sequences",
        "BMC Bioinformatics 26:151",
        "10.1186/s12859-025-06168-3",
    )
    + _entry(
        "TEtrimmer",
        "Qian J, Xue H, Ou S, Mann L, Storer J, Fürtauer L, Heitkam T, "
        "Wildermuth MC, Kusch S, Panstruga R",
        "2025",
        "TEtrimmer: a tool to automate the manual curation of transposable elements",
        "Nature Communications 16",
        "10.1038/s41467-025-63889-y",
    )
)
