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
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_query-def>LINE_frag description</Iteration_query-def>
      <Iteration_hits>
        <Hit>
          <Hit_id>gi|1|gb|MZ092301.1|</Hit_id>
          <Hit_def>PREDICTED: Beta vulgaris subsp. vulgaris uncharacterized LOC125498962, ncRNA</Hit_def>
          <Hit_accession>MZ092301</Hit_accession>
          <Hit_hsps>
            <Hsp>
              <Hsp_bit-score>200.5</Hsp_bit-score>
              <Hsp_evalue>1e-50</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>100</Hsp_query-to>
              <Hsp_hit-from>50</Hsp_hit-from>
              <Hsp_hit-to>149</Hsp_hit-to>
              <Hsp_identity>95</Hsp_identity>
              <Hsp_gaps>2</Hsp_gaps>
              <Hsp_align-len>100</Hsp_align-len>
            </Hsp>
            <Hsp>
              <Hsp_bit-score>80.0</Hsp_bit-score>
              <Hsp_evalue>1e-10</Hsp_evalue>
              <Hsp_query-from>10</Hsp_query-from>
              <Hsp_query-to>60</Hsp_query-to>
              <Hsp_hit-from>200</Hsp_hit-from>
              <Hsp_hit-to>150</Hsp_hit-to>
              <Hsp_identity>48</Hsp_identity>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>50</Hsp_align-len>
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

    out = run_ncbi_online_blast(
        query,
        tmp_path,
        program="blastn",
        database="nt",
        evalue=0.001,
        hitlist_size=10,
        progress=lambda d, t, label: phases.append(label),
    )
    assert out.name == "q.fa_ob.txt"
    text = out.read_text()
    assert "RID=TEST" in text
    assert "MZ092301" in text
    assert "Beta vulgaris" in text
    assert "gi|" not in text
    assert any("connecting to NCBI" in p for p in phases)
    assert any("submitted" in p for p in phases)
    assert any("downloading" in p for p in phases)
    assert any("2 hit row" in p for p in phases)


def test_merge_online_blast_into_annotations(tmp_path: Path) -> None:
    from repeatviewer.backends.ncbi_blast_online import (
        ONLINE_BLAST_SOURCE,
        merge_online_blast_into_annotations,
    )
    from repeatviewer.backends.sequence_features import load_sequence_map_data

    fasta = tmp_path / "cons.fa"
    fasta.write_text(">cons\n" + ("A" * 200) + "\n")
    ob = tmp_path / "cons.fa_ob.txt"
    ob.write_text(
        "# NCBI online BLAST (blastn vs nt); RID=X\n"
        "# columns=qseqid,sseqid,pident,length,mismatch,gapopen,"
        "qstart,qend,sstart,send,sstrand,evalue,bitscore,stitle\n"
        "cons\tACC.1\t95.0\t100\t5\t0\t10\t109\t1\t100\tplus\t1e-20\t200\t"
        "PREDICTED: Beta vulgaris uncharacterized LOC1, ncRNA\n"
        "cons\tACC2.1\t90.0\t80\t8\t0\t50\t129\t5\t84\tminus\t1e-10\t150\t"
        "Beta vulgaris subsp. vulgaris mRNA\n"
    )
    outs = merge_online_blast_into_annotations(fasta, ob)
    assert outs and outs[0].is_file()
    data = load_sequence_map_data(fasta)
    online = [f for f in data.features if f.source == ONLINE_BLAST_SOURCE]
    assert len(online) == 2
    assert all(f.name.startswith("NCBI:") for f in online)
    assert any("Beta vulgaris" in f.name for f in online)
    assert any("LOC1" in f.name for f in online)
    assert {f.start for f in online} == {10, 50}
    assert {round(f.pident, 1) for f in online} == {95.0, 90.0}
    assert {f.subject_id for f in online} == {"ACC.1", "ACC2.1"}
    assert all(f.parent.startswith("NCBI:") for f in online)

    # JSON round-trip keeps the new fields.
    from repeatviewer.backends.sequence_features import SequenceFeature

    for feat in online:
        restored = SequenceFeature.from_dict(feat.to_dict())
        assert restored.pident == feat.pident
        assert restored.subject_id == feat.subject_id
        assert restored.parent == feat.parent

    # Re-run replaces previous online features (no duplicates).
    merge_online_blast_into_annotations(fasta, ob)
    data2 = load_sequence_map_data(fasta)
    online2 = [f for f in data2.features if f.source == ONLINE_BLAST_SOURCE]
    assert len(online2) == 2


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
    shown, mapping = filter_features_for_online_blast_display(
        feats, show_online_blast=True
    )
    assert len(shown) == 3
    assert mapping == [0, 1, 2]
    hidden, mapping2 = filter_features_for_online_blast_display(
        feats, show_online_blast=False
    )
    assert len(hidden) == 2
    assert mapping2 == [0, 2]
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
