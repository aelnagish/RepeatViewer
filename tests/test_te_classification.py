from pathlib import Path

import pandas as pd

from teviewer.backends.blast_genome_plot import prepare_hits, read_repeatmasker_out


def test_blast_query_name_classification_levels() -> None:
    hits = prepare_hits(
        pd.DataFrame(
            [
                {
                    "qseqid": "BvmPIF-II-Ea_Beta_vulgaris#TIR/PIF_Harbinger",
                    "sseqid": "chr1",
                    "pident": 97.0,
                    "sstart": 10,
                    "send": 100,
                    "length": 91,
                }
            ]
        )
    )

    assert hits.loc[0, "order"] == "TIR"
    assert hits.loc[0, "superfamily"] == "PIF_Harbinger"
    assert hits.loc[0, "family"] == "BvmPIF-II-Ea_Beta_vulgaris"


def test_blast_query_id_uses_family_order_superfamily_segments() -> None:
    hits = prepare_hits(
        pd.DataFrame(
            [
                {
                    "qseqid": "AmaS-Au_Beta_patula#Class_I/SINE",
                    "sseqid": "ZR_T717_hq-v2.chr9",
                    "pident": 99.448,
                    "sstart": 34279419,
                    "send": 34279239,
                    "length": 181,
                }
            ]
        )
    )

    assert hits.loc[0, "family"] == "AmaS-Au_Beta_patula"
    assert hits.loc[0, "order"] == "Class_I"
    assert hits.loc[0, "superfamily"] == "SINE"


def test_repeatmasker_classification_levels(tmp_path: Path) -> None:
    output = tmp_path / "sample.out"
    output.write_text(
        "   SW  perc perc perc  query       position in query     matching repeat\n"
        "score  div. del. ins.  sequence     begin end (left)   repeat class/family\n"
        "  100  3.0  0.0  0.0  chr1          10 100 (0) + BvmPIF-II-Ea_Beta_vulgaris TIR/PIF_Harbinger 1 91 (0) 1\n",
        encoding="utf-8",
    )

    hits = read_repeatmasker_out(output)
    assert hits.iloc[0]["order"] == "TIR"
    assert hits.iloc[0]["superfamily"] == "PIF_Harbinger"
    assert hits.iloc[0]["family"] == "BvmPIF-II-Ea_Beta_vulgaris"


def test_repeatmasker_class_without_slash_is_order_and_superfamily(tmp_path: Path) -> None:
    output = tmp_path / "no_slash.out"
    output.write_text(
        "   SW  perc perc perc  query       position in query     matching repeat\n"
        "score  div. del. ins.  sequence     begin end (left)   repeat class/family\n"
        "  100  3.0  0.0  0.0  chr1          10 100 (0) + Example_TE LINE 1 91 (0) 1\n",
        encoding="utf-8",
    )

    hits = read_repeatmasker_out(output)
    assert hits.iloc[0]["order"] == "LINE"
    assert hits.iloc[0]["superfamily"] == "LINE"
    assert hits.iloc[0]["family"] == "Example_TE"
