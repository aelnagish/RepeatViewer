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

    # Re-run replaces previous online features (no duplicates).
    merge_online_blast_into_annotations(fasta, ob)
    data2 = load_sequence_map_data(fasta)
    online2 = [f for f in data2.features if f.source == ONLINE_BLAST_SOURCE]
    assert len(online2) == 2


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
