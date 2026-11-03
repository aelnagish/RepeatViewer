"""Tools > Download from NCBI: query building, parsing, downloads (no network)."""

from __future__ import annotations

import io
import zipfile
from pathlib import Path

import pytest

from repeatviewer.backends import ncbi_download as nd

FEATURE_TABLE = (
    ">Feature ref|NC_045512.2|\n"
    "1\t265\t5'UTR\n"
    "1\t29903\tsource\n"
    "\t\t\torganism\tSARS-CoV-2\n"
    "266\t21555\tgene\n"
    "\t\t\tgene\tORF1ab\n"
    "266\t13468\tCDS\n"
    "13468\t21555\n"
    "\t\t\tproduct\tORF1ab polyprotein\n"
    "\t\t\tnote\tpp1ab; frameshift\n"
    "<300\t>200\tmisc_feature\n"
    "\t\t\tnote\treverse\n"
)


def _rec(acc: str, uid: str = "", **kw) -> nd.NcbiRecord:
    return nd.NcbiRecord(uid=uid or acc, accession=acc, **kw)


# ----------------------------------------------------------------- queries


@pytest.mark.parametrize(
    "token, db",
    [
        ("NC_045512.2", nd.DB_NUCLEOTIDE),
        ("MG815809", nd.DB_NUCLEOTIDE),
        ("OZ559852.1", nd.DB_NUCLEOTIDE),
        ("JAQHRD010000001.1", nd.DB_NUCLEOTIDE),
        ("YP_009724389.1", nd.DB_PROTEIN),
        ("QHD43416.1", nd.DB_PROTEIN),
        ("GCA_029339355.1", nd.DB_ASSEMBLY),
        ("GCF_000005845", nd.DB_ASSEMBLY),
        ("CroSat", None),
        ("actin", None),
        ("ITS2", None),
    ],
)
def test_accession_database(token, db):
    assert nd.accession_database(token) == db


def test_looks_like_accessions_needs_every_token():
    assert nd.looks_like_accessions("NC_045512.2, MG815809\nOZ559852") == [
        "NC_045512.2", "MG815809", "OZ559852"]
    assert nd.looks_like_accessions("MG815809 satellite") == []
    assert nd.guess_database("GCA_029339355.1 GCF_000005845.2") == nd.DB_ASSEMBLY
    assert nd.guess_database("NC_045512.2 YP_009724389.1") is None


def test_build_term_accessions_and_filters():
    assert nd.build_term(nd.NcbiQuery(text="NC_045512.2 MG815809")) == (
        "NC_045512.2[accn] OR MG815809[accn]")
    assert nd.build_term(nd.NcbiQuery(database=nd.DB_ASSEMBLY, text="GCA_029339355.1")) == (
        "GCA_029339355.1")
    term = nd.build_term(nd.NcbiQuery(
        text="CroSat", organism="Crocus sativus", min_length=100, max_length=5000,
        molecule=nd.MOLECULE_GENOMIC, refseq_only=True))
    assert term == ('CroSat AND "Crocus sativus"[Organism] AND 100:5000[SLEN] '
                    'AND refseq[filter] AND biomol_genomic[PROP]')
    assert nd.build_term(nd.NcbiQuery(database=nd.DB_ASSEMBLY, organism="Crocus")) == (
        '"Crocus"[Organism] AND latest[filter]')
    with pytest.raises(nd.NcbiError):
        nd.build_term(nd.NcbiQuery())


def test_clipboard_accessions_from_prose():
    text = "Satellites (GenBank MG815805–MG815809.1) and the genome GCA_029339355.1."
    assert nd.clipboard_accessions(text) == ["MG815805", "MG815809.1", "GCA_029339355.1"]


# ----------------------------------------------------------------- records


def test_record_from_nucleotide_summary():
    rec = nd.record_from_summary(nd.DB_NUCLEOTIDE, {
        "uid": "1", "accessionversion": "MG815809.1", "title": "CroSat-6",
        "organism": "Crocus sativus", "slen": 84, "biomol": "genomic", "moltype": "dna",
        "updatedate": "2018/05/01"})
    assert (rec.accession, rec.length, rec.kind, rec.updated) == (
        "MG815809.1", 84, "genomic DNA", "2018-05-01")
    assert rec.web_url(nd.DB_NUCLEOTIDE).endswith("/nuccore/MG815809.1")


def test_record_from_assembly_summary():
    meta = ('<Stats><Stat category="total_length" sequence_tag="all">3011142691</Stat>'
            '<Stat category="scaffold_count" sequence_tag="all">2562722</Stat></Stats>')
    rec = nd.record_from_summary(nd.DB_ASSEMBLY, {
        "uid": "9", "assemblyaccession": "GCA_029339355.1", "assemblyname": "CS_assembly_2",
        "speciesname": "Crocus sativus", "assemblystatus": "Scaffold", "meta": meta,
        "asmreleasedate_genbank": "2023/03/22 00:00"})
    assert rec.length == 3011142691
    assert rec.extra["sequences"] == 2562722
    assert nd.record_file_name(nd.DB_ASSEMBLY, nd.FORMAT_GENOME, rec) == (
        "GCA_029339355.1_CS_assembly_2.fa")


# -------------------------------------------------------------- file names


def test_names_and_planned_outputs(tmp_path):
    recs = [_rec("NC_045512.2"), _rec("MG815809.1")]
    assert nd.default_file_name(nd.DB_NUCLEOTIDE, nd.FORMAT_FASTA, recs, "CroSat satellite") == (
        "CroSat_satellite.fa")
    assert nd.default_file_name(nd.DB_PROTEIN, nd.FORMAT_FASTA, recs[:1], "") == "NC_045512.2.faa"
    assert nd.default_file_name(nd.DB_NUCLEOTIDE, nd.FORMAT_GENBANK, recs,
                                "NC_045512.2 MG815809") == "NC_045512.2_and_1_more.gb"
    one = nd.planned_outputs(nd.DB_NUCLEOTIDE, nd.FORMAT_FASTA, recs, one_file=True,
                             file_name="x.fa", out_dir=tmp_path)
    assert one == [tmp_path / "x.fa"]
    each = nd.planned_outputs(nd.DB_NUCLEOTIDE, nd.FORMAT_FASTA, recs, one_file=False,
                              file_name="", out_dir=tmp_path)
    assert [p.name for p in each] == ["NC_045512.2.fa", "MG815809.1.fa"]
    asm = _rec("GCF_000005845.2", extra={"assembly_name": "ASM584v2"})
    genome = nd.planned_outputs(nd.DB_ASSEMBLY, nd.FORMAT_GENOME_GFF, [asm], one_file=True,
                                file_name="", out_dir=tmp_path)
    assert [p.name for p in genome] == ["GCF_000005845.2_ASM584v2.fa",
                                        "GCF_000005845.2_ASM584v2.gff3"]
    assert nd.safe_file_stem("Crocus sativus [ITS2]") == "Crocus_sativus_ITS2"


# ------------------------------------------------------------ feature table


def test_feature_table_to_gff3():
    gff = nd.feature_table_to_gff3(FEATURE_TABLE).splitlines()
    assert gff[0] == "##gff-version 3"
    rows = [line.split("\t") for line in gff[1:]]
    assert all(r[0] == "NC_045512.2" for r in rows)
    assert "source" not in {r[2] for r in rows}
    cds = [r for r in rows if r[2] == "CDS"]
    assert [(r[3], r[4]) for r in cds] == [("266", "13468"), ("13468", "21555")]
    assert cds[0][8] == cds[1][8]
    assert "Name=ORF1ab polyprotein" in cds[0][8]
    assert "note=pp1ab%3B frameshift" in cds[0][8]
    misc = next(r for r in rows if r[2] == "misc_feature")
    assert (misc[3], misc[4], misc[6]) == ("200", "300", "-")


# ---------------------------------------------------------------- downloads


class _FakeResp(io.BytesIO):
    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


class _FakeClient(nd.NcbiClient):
    def __init__(self):
        super().__init__()
        self.calls: list[dict] = []

    def eutils(self, endpoint, params, *, post=False):
        self.calls.append({"endpoint": endpoint, **params})
        ids = params["id"].split(",")
        if params.get("rettype") == "ft":
            return _FakeResp(FEATURE_TABLE.encode())
        return _FakeResp("".join(f">{i} test\nACGT\n" for i in ids).encode())


def test_download_one_file_and_per_record(tmp_path):
    recs = [_rec(f"AB{n:06d}.1", uid=str(n)) for n in range(1, 26)]
    client = _FakeClient()
    seen = []
    res = nd.download_records(client, nd.DB_NUCLEOTIDE, recs, tmp_path, fmt=nd.FORMAT_FASTA,
                              one_file=True, file_name="all.fa",
                              progress=lambda d, t, s: seen.append((d, t)))
    assert res.files == [tmp_path / "all.fa"]
    assert (tmp_path / "all.fa").read_text().count(">") == 25
    assert len(client.calls) == 2  # batches of 20
    assert seen[-1] == (25, 25)
    assert not list(tmp_path.glob("*.part"))

    res = nd.download_records(client, nd.DB_PROTEIN, recs[:2], tmp_path, fmt=nd.FORMAT_FASTA,
                              one_file=False,
                              targets={str(tmp_path / "AB000002.1.faa"): tmp_path / "renamed.faa"})
    assert [p.name for p in res.files] == ["AB000001.1.faa", "renamed.faa"]


def test_download_with_features_returns_gff3(tmp_path):
    res = nd.download_records(_FakeClient(), nd.DB_NUCLEOTIDE, [_rec("NC_045512.2", uid="7")],
                              tmp_path, fmt=nd.FORMAT_FASTA_FEATURES, one_file=False)
    fasta = res.files[0]
    assert res.fasta_with_features == [fasta]
    assert "ORF1ab" in res.feature_tables[str(fasta)]


def test_download_cancel_leaves_no_partial_file(tmp_path):
    with pytest.raises(nd.NcbiCancelled):
        nd.download_records(_FakeClient(), nd.DB_NUCLEOTIDE, [_rec("AB000001.1")], tmp_path,
                            one_file=False, cancel=lambda: True)
    assert list(tmp_path.iterdir()) == []


def test_extract_assembly_zip(tmp_path):
    zpath = tmp_path / "d.zip"
    with zipfile.ZipFile(zpath, "w") as z:
        z.writestr("ncbi_dataset/data/GCF_1/GCF_1_genomic.fna", ">chr1\nACGT\n")
        z.writestr("ncbi_dataset/data/GCF_1/genomic.gff", "##gff-version 3\n")
        z.writestr("README.md", "x")
    fa, gff = tmp_path / "g.fa", tmp_path / "g.gff3"
    assert nd.extract_assembly_zip(zpath, fa, gff) == (True, True)
    assert fa.read_text() == ">chr1\nACGT\n"
    assert gff.read_text().startswith("##gff-version 3")


# ------------------------------------------------------------------- window


def test_dialog_ticks_names_and_options(qapp, tmp_path):
    from PySide6.QtCore import Qt
    from PySide6.QtWidgets import QGroupBox

    from repeatviewer.job_options import (
        DIALOG_TAB_OPTIONS,
        NCBI_RIGHT_CLICK_CLIPBOARD,
        NcbiDownloadOptions,
    )
    from repeatviewer.views.ncbi_download_dialog import NcbiDownloadDialog, missing_accessions

    dlg = NcbiDownloadDialog(initial=NcbiDownloadOptions(), folder=tmp_path)
    labels = [dlg._db_combo.itemText(i) for i in range(dlg._db_combo.count())]
    assert "Genome Assembly" in labels
    assert "..." not in "".join(labels)
    assert dlg._db_combo.view().textElideMode() == Qt.TextElideMode.ElideNone
    # Options: no right-click chooser; window shrinks to content.
    search_h = dlg.height()
    dlg._tab_bar.setCurrentIndex(dlg._TAB_OPTIONS)
    qapp.processEvents()
    assert dlg.height() < search_h
    assert not any(
        "Right-click" in (g.title() or "")
        for g in dlg._options_page.findChildren(QGroupBox)
    )
    assert dlg.options().right_click_action == NCBI_RIGHT_CLICK_CLIPBOARD
    assert dlg.current_dialog_tab() == DIALOG_TAB_OPTIONS
    dlg._tab_bar.setCurrentIndex(dlg._TAB_SEARCH)
    dlg._search_edit.setText("CroSat")
    dlg.add_records([_rec("MG815809.1", uid="1", length=84), _rec("MG815808.1", uid="2")])
    assert dlg._table.rowCount() == 2
    assert dlg.chosen_records() == []
    dlg._on_tick_all(True)
    assert [r.accession for r in dlg.chosen_records()] == ["MG815809.1", "MG815808.1"]
    assert dlg.file_name() == "CroSat.fa"
    dlg._per_record.setChecked(True)
    assert not dlg.one_file()
    dlg._toggle_row(1)
    assert [r.accession for r in dlg.chosen_records()] == ["MG815809.1"]
    # Pasting a protein accession switches the database.
    dlg._search_edit.setText("YP_009724389.1")
    assert dlg.database() == nd.DB_PROTEIN
    sub = tmp_path / "sub"
    sub.mkdir()
    assert dlg.apply_path_from_panel(sub / "x.fa")
    assert dlg.out_dir() == sub
    opts = dlg.options()
    assert opts.database == nd.DB_PROTEIN and opts.one_file is False
    assert missing_accessions(["MG815809", "XX999999"], [_rec("MG815809.1")]) == ["XX999999"]
    dlg.reject()


def test_progress_detail_for_the_status_bar():
    from repeatviewer.ncbi_download_actions import split_progress_detail

    assert split_progress_detail("GCF_1.1: 161.0 MB downloaded (0 done, 3 left)") == (
        "GCF_1.1", "161.0 MB")
    assert split_progress_detail("NC_045512.2 (0 done, 1 left)") == ("NC_045512.2", "")
    assert split_progress_detail("GCF_1.1: unpacking (1 done, 2 left)") == ("GCF_1.1", "unpacking")
