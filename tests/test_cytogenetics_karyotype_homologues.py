"""Homologues that differ, B chromosomes, missing homologues and band presets in a karyotype."""

from pathlib import Path

from repeatviewer.backends.cytogenetics import karyotype as kt


def _k() -> kt.Karyotype:
    k = kt.new_karyotype("test", 2, 3)
    k.chromosomes[0].length = 20
    k.chromosomes[0].centromere_index = 37.5
    k.chromosomes[0].bands.append(kt.Band("DAPI+ heterochromatin", 0.41, 0.59))
    return k


def test_separate_homologue_keeps_own_values():
    k = _k()
    c = k.chromosomes[0]
    kt.set_separate(c, 1, True)
    kt.target(c, 1).length = 21
    assert kt.homologue_view(c, 1).length == 21
    assert kt.homologue_view(c, 0).length == 20
    assert kt.differs(c)
    kt.set_separate(c, 1, False)
    assert not kt.differs(c) and c.homologues == {}


def test_round_trip_and_formula():
    k = _k()
    c = k.chromosomes[0]
    kt.set_separate(c, 1, True)
    kt.target(c, 1).bands = []
    kt.set_missing(c, 2, True)
    kt.add_b_chromosome(k, c, 0)
    k2 = kt.Karyotype.from_dict(k.to_dict())
    c2 = k2.chromosomes[0]
    assert kt.homologue_view(c2, 1).bands == []
    assert kt.is_missing(c2, 2)
    assert k2.chromosomes[-1].b_chromosome and k2.chromosomes[-1].name == "B1"
    assert k2.formula == "2n = 3x - 1 + 1B = 6"
    assert k2.homologue_labels(c2) == ["1a", "1b"]
    k2.set_basic_number(3)
    assert [x.name for x in k2.chromosomes] == ["1", "2", "3", "B1"]


def test_mirror_and_copy_bands():
    k = _k()
    c = k.chromosomes[0]
    c.bands = [kt.Band("CMA+ band", 0.0, 0.1)]
    kt.mirror_bands(c, 2)
    b = kt.homologue_view(c, 2).bands[0]
    assert (round(b.start, 3), round(b.end, 3)) == (0.9, 1.0)
    kt.copy_bands_to_others(k, c, 2)
    assert c.bands[0].start == b.start and not kt.is_separate(c, 2)


def test_presets(tmp_path: Path):
    k = _k()
    kt.save_preset(tmp_path, "Crocus sativus", "rDNA", k)
    assert kt.list_presets(tmp_path, "Crocus sativus") == ["rDNA"]
    k.chromosomes[0].bands = []
    assert kt.apply_preset(tmp_path, "Crocus sativus", "rDNA", k) == []
    assert k.chromosomes[0].bands[0].kind == "DAPI+ heterochromatin"
