"""Parse nhmmscan --tblout used for DNA HMM / Dfam similarity annotation."""

from __future__ import annotations

from repeatviewer.backends.pyhmmer_search import parse_nhmmscan_tblout


_SAMPLE = """\
# target name                   accession     query name                               accession  hmmfrom hmm to alifrom  ali to envfrom  env to  modlen strand   E-value  score  bias  description of target
#           -------------------    ----------                     -------------------- ---------- ------- ------- ------- ------- ------- ------- ------- ------ --------- ------ ----- ---------------------
Copia-44_rCom                   DF003584236.2 Tork0001485_Beta_vulgaris#Ty1_copia/Tork -              901    4611     975    4605     950    4626    5219    +    4.5e-235  785.4 168.5  Repetitive sequence family
Copia-58_rCom                   DF003584251.2 Tork0001485_Beta_vulgaris#Ty1_copia/Tork -             1775    4648    1700    4603    1680    4623    5470    +    2.8e-202  676.3 123.2  Repetitive sequence family
#
# Program:         nhmmscan
"""


def test_parse_nhmmscan_tblout_env_coords_and_scores() -> None:
    hits = parse_nhmmscan_tblout(_SAMPLE)
    assert len(hits) == 2
    score, evalue, query_id, hmm_name, env_from, env_to = hits[0]
    assert hmm_name == "Copia-44_rCom"
    assert query_id == "Tork0001485_Beta_vulgaris#Ty1_copia/Tork"
    assert env_from == 950
    assert env_to == 4626
    assert score == 785.4
    assert evalue is not None and evalue < 1e-200
    assert hits[1][3] == "Copia-58_rCom"


def test_parse_nhmmscan_tblout_ignores_junk() -> None:
    assert parse_nhmmscan_tblout("# only comments\n\n") == []
    assert parse_nhmmscan_tblout("too few fields\n") == []
