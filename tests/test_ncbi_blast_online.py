"""NCBI online BLAST: XML→tabular conversion and progress phases."""

from __future__ import annotations

from pathlib import Path

import pytest

from repeatviewer.backends.ncbi_blast_online import (
    NcbiBlastOnlineError,
    blast_xml_to_tabular,
    run_ncbi_online_blast,
)


_SAMPLE_XML = """<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov">
<BlastOutput>
  <BlastOutput_program>blastn</BlastOutput_program>
  <BlastOutput_version>BLASTN 2.16.0+</BlastOutput_version>
  <BlastOutput_db>nt</BlastOutput_db>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_query-def>LINE_frag description</Iteration_query-def>
      <Iteration_query-len>200</Iteration_query-len>
      <Iteration_hits>
        <Hit>
          <Hit_id>gi|1|gb|MZ092301.1|</Hit_id>
          <Hit_def>PREDICTED: Beta vulgaris subsp. vulgaris uncharacterized LOC125498962, ncRNA</Hit_def>
          <Hit_accession>MZ092301</Hit_accession>
          <Hit_len>1500</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_bit-score>200.5</Hsp_bit-score>
              <Hsp_score>400</Hsp_score>
              <Hsp_evalue>1e-50</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>100</Hsp_query-to>
              <Hsp_hit-from>50</Hsp_hit-from>
              <Hsp_hit-to>149</Hsp_hit-to>
              <Hsp_identity>95</Hsp_identity>
              <Hsp_positive>95</Hsp_positive>
              <Hsp_gaps>2</Hsp_gaps>
              <Hsp_align-len>100</Hsp_align-len>
              <Hsp_qseq>ACGTACGTAC</Hsp_qseq>
              <Hsp_hseq>ACGTACGTAC</Hsp_hseq>
              <Hsp_midline>||||||||||</Hsp_midline>
            </Hsp>
            <Hsp>
              <Hsp_bit-score>80.0</Hsp_bit-score>
              <Hsp_score>160</Hsp_score>
              <Hsp_evalue>1e-10</Hsp_evalue>
              <Hsp_query-from>10</Hsp_query-from>
              <Hsp_query-to>60</Hsp_query-to>
              <Hsp_hit-from>200</Hsp_hit-from>
              <Hsp_hit-to>150</Hsp_hit-to>
              <Hsp_identity>48</Hsp_identity>
              <Hsp_positive>48</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>50</Hsp_align-len>
              <Hsp_qseq>ACGTAC</Hsp_qseq>
              <Hsp_hseq>ACGTAT</Hsp_hseq>
              <Hsp_midline>||||| </Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>
"""


def test_blast_xml_to_tabular_rows() -> None:
    table, n = blast_xml_to_tabular(_SAMPLE_XML)
    assert n == 2
    lines = [ln for ln in table.splitlines() if ln.strip()]
    assert len(lines) == 2
    cols0 = lines[0].split("\t")
    assert cols0[0] == "LINE_frag"
    assert cols0[1] == "MZ092301"
    assert cols0[2] == "95.000"
    assert cols0[3] == "100"
    assert cols0[10] == "plus"
    assert cols0[11] == "1e-50"
    assert "Beta vulgaris" in cols0[13]
    assert "LOC125498962" in cols0[13]
    # Reverse-strand HSP: coordinates normalised, strand=minus.
    cols1 = lines[1].split("\t")
    assert cols1[8] == "150"
    assert cols1[9] == "200"
    assert cols1[10] == "minus"


def test_blast_xml_to_tabular_empty() -> None:
    empty = """<?xml version="1.0"?>
    <BlastOutput><BlastOutput_iterations>
      <Iteration><Iteration_query-def>q</Iteration_query-def>
      <Iteration_hits/></Iteration>
    </BlastOutput_iterations></BlastOutput>"""
    table, n = blast_xml_to_tabular(empty)
    assert n == 0
    assert table == ""


def test_run_ncbi_online_blast_writes_hits(tmp_path: Path, monkeypatch) -> None:
    query = tmp_path / "q.fa"
    query.write_text(">q\nACGTACGTACGTACGTACGT\n")

    put = "    RID = TEST RID\n    RTOE = 1\n"
    info = "Status=READY\n"
    phases: list[str] = []

    def fake_post(params):
        assert params["CMD"] == "Put"
        assert params["PROGRAM"] == "blastn"
        return put

    def fake_get(params, timeout=60.0):
        if params.get("FORMAT_OBJECT") == "SearchInfo":
            return info
        assert params.get("FORMAT_TYPE") == "XML"
        return _SAMPLE_XML

    monkeypatch.setattr(
        "repeatviewer.backends.ncbi_blast_online._ncbi_post", fake_post
    )
    monkeypatch.setattr(
        "repeatviewer.backends.ncbi_blast_online._ncbi_get", fake_get
    )
    monkeypatch.setattr(
        "repeatviewer.backends.ncbi_blast_online.time.sleep", lambda *_a, **_k: None
    )

    out_paths = run_ncbi_online_blast(
        query,
        tmp_path,
        program="blastn",
        database="nt",
        evalue=0.001,
        hitlist_size=10,
        write_text=True,
        progress=lambda d, t, label: phases.append(label),
    )
    assert len(out_paths) == 2
    html_path, txt_path = out_paths
    assert html_path.name == "q.fa_ob.html"
    assert txt_path.name == "q.fa_ob.txt"
    html = html_path.read_text()
    assert "Descriptions" in html
    assert "Graphic Summary" in html
    assert "Lineage" in html
    assert 'data-tab="tab-msa"' not in html
    assert 'data-tab="tab-taxonomy"' not in html
    assert "MZ092301" in html
    text = txt_path.read_text()
    assert "RID=TEST" in text
    assert "MZ092301" in text
    assert "Beta vulgaris" in text
    assert "gi|" not in text
    assert any("connecting to NCBI" in p for p in phases)
    assert any("submitted" in p for p in phases)
    assert any("downloading" in p for p in phases)
    assert any("2 hit row" in p for p in phases)


def test_online_blast_html_report_from_xml() -> None:
    from repeatviewer.backends.online_blast_report import (
        build_online_blast_html,
        parse_blast_output_xml,
    )

    result = parse_blast_output_xml(
        _SAMPLE_XML,
        rid="TEST",
        query_name="q.fa",
        program="blastn",
        database="nt",
    )
    assert result.primary is not None
    assert result.primary.query_len == 200
    assert len(result.primary.hits) == 1
    assert result.primary.hits[0].accession == "MZ092301"
    assert len(result.primary.hits[0].hsps) == 2
    html = build_online_blast_html(result, fetch_taxonomy=False)
    assert 'id="tab-descriptions"' in html
    assert 'id="tab-graphic"' in html
    assert 'id="tab-lineage"' in html
    assert 'id="tab-msa"' not in html
    assert 'id="tab-alignments"' not in html
    assert 'id="tab-taxonomy"' not in html
    assert "lineage-table" in html
    assert "MZ092301" in html
    assert "200.5" in html or "200" in html
    assert "Beta vulgaris" in html
    assert "<svg" in html
    assert "foreignObject" in html
    assert "Alignment Scores" in html
    assert "-webkit-line-clamp" not in html
    assert "org-header" in html
    assert "sortable" in html
    assert "group-by-organism" in html
    assert "Organism</b> is the species" not in html
    assert "Organism" in html
    assert "Number of Hits" in html
    assert "Blast Name" in html
    assert "data-lineage=" in html
    assert "<b>Description:</b>" not in html
    assert "Description:</b>" not in html


def test_taxonomy_parse_efetch_includes_eukaryota() -> None:
    from repeatviewer.backends.online_blast_report import _taxonomy_parse_efetch

    xml = """<?xml version="1.0"?>
    <TaxaSet><Taxon>
      <TaxId>3555</TaxId>
      <ScientificName>Beta vulgaris subsp. vulgaris</ScientificName>
      <LineageEx>
        <Taxon><TaxId>131567</TaxId><ScientificName>cellular organisms</ScientificName></Taxon>
        <Taxon><TaxId>2759</TaxId><ScientificName>Eukaryota</ScientificName></Taxon>
        <Taxon><TaxId>33090</TaxId><ScientificName>Viridiplantae</ScientificName></Taxon>
        <Taxon><TaxId>161934</TaxId><ScientificName>Beta vulgaris</ScientificName></Taxon>
      </LineageEx>
    </Taxon></TaxaSet>"""
    parsed = _taxonomy_parse_efetch(xml)
    assert "3555" in parsed
    names = parsed["3555"].lineage_names
    assert names[0] == "Eukaryota"
    assert "cellular organisms" not in names
    assert names[-1] == "Beta vulgaris subsp. vulgaris"


def test_lineage_uses_ncbi_taxonomy_path(monkeypatch) -> None:
    from repeatviewer.backends.online_blast_report import (
        build_online_blast_html,
        parse_blast_output_xml,
    )
    import repeatviewer.backends.online_blast_report as report

    result = parse_blast_output_xml(
        _SAMPLE_XML, rid="TEST", query_name="q.fa", program="blastn", database="nt"
    )

    def fake_lookup(organisms, enable_network=True):
        return {
            "Beta vulgaris subsp. vulgaris": [
                "Eukaryota",
                "Viridiplantae",
                "Streptophyta",
                "Beta",
                "Beta vulgaris",
                "Beta vulgaris subsp. vulgaris",
            ]
        }

    monkeypatch.setattr(report, "_lookup_organism_lineages", fake_lookup)
    html = build_online_blast_html(result, fetch_taxonomy=True)
    assert "Eukaryota" in html
    assert "eukaryotes" in html
    assert "Viridiplantae" in html


def test_alignment_score_color_bins() -> None:
    from repeatviewer.backends.online_blast_report import _alignment_score_color

    assert _alignment_score_color(10) == "#000000"
    assert _alignment_score_color(45) == "#0000FF"
    assert _alignment_score_color(60) == "#00FF00"
    assert _alignment_score_color(100) == "#FF00FF"
    assert _alignment_score_color(200) == "#FF0000"


def test_pack_hsps_into_lanes_splits_overlaps() -> None:
    from repeatviewer.backends.online_blast_report import BlastHsp, _pack_hsps_into_lanes

    def hsp(qs: int, qe: int, score: float = 10.0) -> BlastHsp:
        return BlastHsp(
            bit_score=score,
            score=score,
            evalue="0",
            identity=1,
            positive=1,
            gaps=0,
            align_len=qe - qs + 1,
            qstart=qs,
            qend=qe,
            sstart=1,
            send=qe - qs + 1,
            sstrand="plus",
            pident=90.0,
        )

    lanes = _pack_hsps_into_lanes([hsp(1, 100), hsp(50, 150), hsp(200, 250)])
    assert len(lanes) == 2
    assert {(h.qstart, h.qend) for h in lanes[0]} == {(1, 100), (200, 250)}
    assert {(h.qstart, h.qend) for h in lanes[1]} == {(50, 150)}


def test_organism_lineage_parts() -> None:
    from repeatviewer.backends.online_blast_report import _organism_lineage_parts

    assert _organism_lineage_parts("Beta vulgaris subsp. vulgaris") == [
        "Beta",
        "Beta vulgaris",
        "Beta vulgaris subsp. vulgaris",
    ]
    assert _organism_lineage_parts("Spinacia oleracea") == [
        "Spinacia",
        "Spinacia oleracea",
    ]


def test_strip_online_blast_from_annotations(tmp_path: Path) -> None:
    from repeatviewer.backends.ncbi_blast_online import (
        ONLINE_BLAST_SOURCE,
        strip_online_blast_from_annotations,
    )
    from repeatviewer.backends.sequence_features import (
        SequenceFeature,
        load_sequence_map_data,
        save_annotations_multi,
    )

    fasta = tmp_path / "cons.fa"
    fasta.write_text(">cons\n" + ("A" * 200) + "\n")
    feats = [
        SequenceFeature(
            1, 20, "LTR", feature_type="LTR", source="repeatviewer_internal"
        ),
        SequenceFeature(
            10,
            50,
            "NCBI:ACC",
            feature_type="domain",
            source=ONLINE_BLAST_SOURCE,
            pident=90.0,
            subject_id="ACC",
            parent="NCBI:ACC",
        ),
    ]
    save_annotations_multi(
        fasta,
        features_by_id={"cons": feats},
        lengths_by_id={"cons": 200},
        record_ids=["cons"],
    )
    outs = strip_online_blast_from_annotations(fasta)
    assert outs
    data = load_sequence_map_data(fasta)
    assert all(f.source != ONLINE_BLAST_SOURCE for f in data.features)
    assert any(f.name == "LTR" for f in data.features)


def test_online_blast_keeps_all_hsps_for_one_subject(tmp_path: Path) -> None:
    """max_hits caps subjects, not HSP rows — multi-HSP subjects stay linked."""
    from repeatviewer.backends.ncbi_blast_online import (
        ONLINE_BLAST_SOURCE,
        online_blast_features_by_query,
    )

    ob = tmp_path / "q.fa_ob.txt"
    ob.write_text(
        "# columns=qseqid,sseqid,pident,length,mismatch,gapopen,"
        "qstart,qend,sstart,send,sstrand,evalue,bitscore,stitle\n"
        "q\tACC.1\t95.0\t100\t5\t0\t10\t109\t1\t100\tplus\t1e-20\t200\t"
        "Beta vulgaris hit A\n"
        "q\tACC.1\t88.0\t50\t6\t0\t200\t249\t150\t199\tplus\t1e-8\t80\t"
        "Beta vulgaris hit A\n"
        "q\tOTHER.1\t99.0\t40\t1\t0\t1\t40\t1\t40\tplus\t1e-5\t50\t"
        "other subject\n"
    )
    by_query = online_blast_features_by_query(ob, max_hits=1)
    feats = by_query["q"]
    # Only ACC.1 subject kept (best bitscore 200), but BOTH of its HSPs.
    assert len(feats) == 2
    assert {f.subject_id for f in feats} == {"ACC.1"}
    assert {f.parent for f in feats} == {"NCBI:ACC.1"}
    assert {f.start for f in feats} == {10, 200}
    assert {round(f.pident, 1) for f in feats} == {95.0, 88.0}


def test_pack_online_blast_one_row_per_subject() -> None:
    from repeatviewer.backends.ncbi_blast_online import ONLINE_BLAST_SOURCE
    from repeatviewer.backends.sequence_features import (
        SequenceFeature,
        pack_feature_display_rows,
    )

    feats = [
        SequenceFeature(1, 10, "LTR", feature_type="LTR", source="repeatviewer_internal"),
        SequenceFeature(
            10,
            50,
            "hit A",
            feature_type="domain",
            source=ONLINE_BLAST_SOURCE,
            pident=95.0,
            subject_id="ACC.A",
            parent="NCBI:ACC.A",
        ),
        SequenceFeature(
            200,
            240,
            "hit A2",
            feature_type="domain",
            source=ONLINE_BLAST_SOURCE,
            pident=88.0,
            subject_id="ACC.A",
            parent="NCBI:ACC.A",
        ),
        SequenceFeature(
            60,
            90,
            "hit B",
            feature_type="domain",
            source=ONLINE_BLAST_SOURCE,
            pident=70.0,
            subject_id="ACC.B",
            parent="NCBI:ACC.B",
        ),
    ]
    rows, count = pack_feature_display_rows(feats, split_by_source=True)
    # Non-overlapping HSPs for ACC.A stay on one row; ACC.B is another.
    assert rows[1] == rows[2]
    assert rows[1] != rows[3]
    assert count >= 3


def test_pack_online_blast_overlapping_hsps_stack() -> None:
    from repeatviewer.backends.ncbi_blast_online import ONLINE_BLAST_SOURCE
    from repeatviewer.backends.sequence_features import (
        SequenceFeature,
        pack_feature_display_rows,
    )

    feats = [
        SequenceFeature(
            10,
            100,
            "hit A",
            feature_type="domain",
            source=ONLINE_BLAST_SOURCE,
            pident=95.0,
            subject_id="ACC.A",
            parent="NCBI:ACC.A",
        ),
        SequenceFeature(
            50,
            150,
            "hit A2",
            feature_type="domain",
            source=ONLINE_BLAST_SOURCE,
            pident=88.0,
            subject_id="ACC.A",
            parent="NCBI:ACC.A",
        ),
    ]
    rows, count = pack_feature_display_rows(feats, split_by_source=True)
    assert rows[0] != rows[1]
    assert count == 2


def test_filter_features_for_online_blast_display() -> None:
    from repeatviewer.backends.ncbi_blast_online import ONLINE_BLAST_SOURCE
    from repeatviewer.backends.sequence_features import SequenceFeature
    from repeatviewer.widgets.sequence_map_viewer import (
        filter_features_for_online_blast_display,
    )

    feats = [
        SequenceFeature(1, 10, "LTR", feature_type="LTR", source="repeatviewer_internal"),
        SequenceFeature(
            20,
            30,
            "NCBI:ACC",
            feature_type="domain",
            source=ONLINE_BLAST_SOURCE,
            pident=90.0,
            subject_id="ACC",
            parent="NCBI:ACC",
        ),
        SequenceFeature(40, 50, "ORF", feature_type="ORF", source="repeatviewer_internal"),
    ]
    # Overlay is retired: display always hides online BLAST features.
    hidden, mapping = filter_features_for_online_blast_display(
        feats, show_online_blast=False
    )
    assert len(hidden) == 2
    assert mapping == [0, 2]
    assert all(f.source != ONLINE_BLAST_SOURCE for f in hidden)


def test_pident_color_hex_scales() -> None:
    from repeatviewer.backends.annotation_colors import (
        pident_color_hex,
        pident_scale_bounds,
    )

    zmin, zmax = pident_scale_bounds([80.0, 95.0])
    assert zmin == 80.0
    assert zmax == 95.0
    low = pident_color_hex(80.0, zmin, zmax)
    high = pident_color_hex(95.0, zmin, zmax)
    assert low.startswith("#") and high.startswith("#")
    assert low != high
    assert pident_scale_bounds([]) == (50.0, 100.0)


def test_run_ncbi_online_blast_rejects_bad_program(tmp_path: Path) -> None:
    query = tmp_path / "q.fa"
    query.write_text(">q\nACGT\n")
    with pytest.raises(NcbiBlastOnlineError, match="Unsupported"):
        run_ncbi_online_blast(query, tmp_path, program="notablast")


def test_online_blast_databases_cover_ncbi_suite() -> None:
    from repeatviewer.backends.ncbi_blast_online import (
        ONLINE_BLAST_DATABASES,
        default_online_database,
    )
    from repeatviewer.job_options import BLAST_ONLINE_ALGORITHMS

    assert BLAST_ONLINE_ALGORITHMS == (
        "blastn",
        "blastp",
        "blastx",
        "tblastn",
        "tblastx",
    )
    for program in BLAST_ONLINE_ALGORITHMS:
        assert program in ONLINE_BLAST_DATABASES
        assert default_online_database(program)
        keys = {key for key, _ in ONLINE_BLAST_DATABASES[program]}
        assert default_online_database(program) in keys
    nt_keys = {key for key, _ in ONLINE_BLAST_DATABASES["blastn"]}
    for name in ("core_nt", "nt", "refseq_rna", "wgs", "est", "sra"):
        assert name in nt_keys
    prot_keys = {key for key, _ in ONLINE_BLAST_DATABASES["blastp"]}
    for name in ("nr", "swissprot", "refseq_protein", "pdbaa"):
        assert name in prot_keys


def test_online_blast_nt_sets_no_override(tmp_path: Path, monkeypatch) -> None:
    query = tmp_path / "q.fa"
    query.write_text(">q\nACGTACGTACGTACGTACGT\n")
    put_seen: dict[str, str] = {}

    def fake_post(params):
        put_seen.update(params)
        return "    RID = TEST RID\n    RTOE = 1\n"

    def fake_get(params, timeout=60.0):
        if params.get("FORMAT_OBJECT") == "SearchInfo":
            return "Status=READY\n"
        return _SAMPLE_XML

    monkeypatch.setattr(
        "repeatviewer.backends.ncbi_blast_online._ncbi_post", fake_post
    )
    monkeypatch.setattr(
        "repeatviewer.backends.ncbi_blast_online._ncbi_get", fake_get
    )
    monkeypatch.setattr(
        "repeatviewer.backends.ncbi_blast_online.time.sleep", lambda *_a, **_k: None
    )
    run_ncbi_online_blast(
        query, tmp_path, program="blastn", database="nt"
    )
    assert put_seen.get("DATABASE") == "nt"
    assert put_seen.get("NO_DATABASE_OVERRIDE") == "true"


def test_run_batch_online_blast_stop_keeps_finished(tmp_path: Path, monkeypatch) -> None:
    """Stop mid-batch still returns outputs already written for finished queries."""
    from repeatviewer.backends import curation
    from repeatviewer.backends.job_cancel import JobCancellation, JobCancelledError
    from repeatviewer.job_options import BLAST_PROGRAM_ONLINE, BlastOptions
    from repeatviewer.project import RepeatViewerProject

    files = []
    for i in range(3):
        p = tmp_path / f"q{i}.fa"
        p.write_text(f">q{i}\nACGTACGTACGTACGTACGT\n")
        files.append(p)

    calls = {"n": 0}
    ready: list[str] = []

    def fake_run_on_sequence(func_name, seq_path, project, **kwargs):
        calls["n"] += 1
        cancel = kwargs.get("cancel")
        if calls["n"] == 1:
            out = seq_path.parent / f"{seq_path.name}_ob.html"
            out.write_text("<html>ok</html>")
            return [str(out)]
        if cancel is not None:
            cancel.request()
        raise JobCancelledError("Job stopped by user")

    monkeypatch.setattr(curation, "run_on_sequence", fake_run_on_sequence)
    result = curation.run_batch(
        "blast",
        files,
        RepeatViewerProject(),
        lambda *_a, **_k: None,
        blast_options=BlastOptions(program=BLAST_PROGRAM_ONLINE),
        cancel=JobCancellation(),
        on_outputs_ready=lambda paths: ready.extend(paths),
    )
    assert result.get("cancelled") is True
    assert len(result["created_outputs"]) == 1
    assert Path(result["created_outputs"][0]).name == "q0.fa_ob.html"
    assert ready and Path(ready[0]).name == "q0.fa_ob.html"


def test_lineage_merges_shallow_beta_under_eukaryota() -> None:
    """Hits that miss Taxonomy must not create a second root-level Beta."""
    from repeatviewer.backends.online_blast_report import (
        BlastHit,
        BlastHsp,
        BlastIteration,
        _build_lineage_rows,
    )

    def hit(title: str, score: float) -> BlastHit:
        return BlastHit(
            accession=title[:12],
            hit_id="",
            title=title,
            length=100,
            hsps=[
                BlastHsp(
                    bit_score=score,
                    score=score,
                    evalue="0",
                    identity=10,
                    positive=10,
                    gaps=0,
                    align_len=10,
                    qstart=1,
                    qend=10,
                    sstart=1,
                    send=10,
                    sstrand="plus",
                    pident=100.0,
                )
            ],
        )

    it = BlastIteration(
        query_id="q",
        query_def="q",
        query_len=100,
        hits=[
            hit("Beta vulgaris subsp. vulgaris transposon VulMITE I", 200),
            hit("Beta vulgaris subsp. maritima some gene", 80),
        ],
    )
    org_lineages = {
        "Beta vulgaris subsp. vulgaris": [
            "Eukaryota",
            "Viridiplantae",
            "eudicotyledons",
            "Amaranthaceae",
            "Beta",
            "Beta vulgaris",
            "Beta vulgaris subsp. vulgaris",
        ]
        # subsp. maritima intentionally missing → would have been shallow Beta
    }
    rows = _build_lineage_rows(it, org_lineages=org_lineages)
    beta_rows = [r for r in rows if r.name == "Beta"]
    assert len(beta_rows) == 1
    assert beta_rows[0].depth > 0
    assert any(r.name == "Eukaryota" for r in rows)
    assert any(r.name == "Beta vulgaris subsp. maritima" for r in rows)
    maritima = next(r for r in rows if r.name == "Beta vulgaris subsp. maritima")
    # Grafted under the same deep tree, not a separate root branch.
    assert maritima.depth == beta_rows[0].depth + 2


def test_guess_organism_stops_at_binomial() -> None:
    from repeatviewer.backends.online_blast_report import _guess_organism

    assert (
        _guess_organism(
            "Chenopodium quinoa xyloglucan galactosyltransferase XLT2-like, mRNA"
        )
        == "Chenopodium quinoa"
    )
    assert (
        _guess_organism(
            "PREDICTED: Beta vulgaris subsp. vulgaris xyloglucan galactosyltransferase "
            "XLT2 (LOC104906587), mRNA"
        )
        == "Beta vulgaris subsp. vulgaris"
    )
    assert (
        _guess_organism(
            "Beta vulgaris subsp. vulgaris transposon VulMITE II-2, cultivar Rosamona"
        )
        == "Beta vulgaris subsp. vulgaris"
    )


def test_descriptions_cluster_similar_titles() -> None:
    from repeatviewer.backends.online_blast_report import (
        BlastHit,
        BlastHsp,
        BlastIteration,
        _descriptions_table,
    )

    def hit(title: str, score: float) -> BlastHit:
        return BlastHit(
            accession=f"ACC{int(score)}",
            hit_id="",
            title=title,
            length=100,
            hsps=[
                BlastHsp(
                    bit_score=score,
                    score=score,
                    evalue="0",
                    identity=10,
                    positive=10,
                    gaps=0,
                    align_len=10,
                    qstart=1,
                    qend=10,
                    sstart=1,
                    send=10,
                    sstrand="plus",
                    pident=100.0,
                )
            ],
        )

    it = BlastIteration(
        query_id="q",
        query_def="q",
        query_len=100,
        hits=[
            hit("Beta vulgaris subsp. vulgaris transposon VulMITE II-2", 300),
            hit(
                "PREDICTED: Beta vulgaris subsp. vulgaris xyloglucan galactosyltransferase",
                250,
            ),
            hit("Beta vulgaris subsp. vulgaris transposon VulMITE II-1", 200),
            hit("Chenopodium quinoa xyloglucan galactosyltransferase XLT2", 180),
        ],
    )
    html = _descriptions_table(it)
    assert html.count("org-header") == 2
    assert 'id="group-by-organism"' in html
    assert "Group by organism" in html
    assert 'data-group-by-organism="1"' in html
    assert "Beta vulgaris subsp. vulgaris <span class=\"org-count\">(3 hits)</span>" in html
    assert "Chenopodium quinoa <span class=\"org-count\">(1 hit)</span>" in html
    assert "th class=\"sortable\"" in html or 'class="sortable"' in html
    beta_block = html.split("Beta vulgaris subsp. vulgaris <span")[1].split(
        "Chenopodium quinoa <span"
    )[0]
    vul1 = beta_block.find("VulMITE II-1")
    vul2 = beta_block.find("VulMITE II-2")
    gene = beta_block.find("xyloglucan")
    assert vul1 != -1 and vul2 != -1 and gene != -1
    assert (vul1 < gene and vul2 < gene) or (vul1 > gene and vul2 > gene)
