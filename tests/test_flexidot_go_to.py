from repeatviewer.backends.flexidot_go_to import parse_dotplot_go_to_targets


def test_previous_flexidot_x_axis_payload() -> None:
    payload = {
        "x_seq": "seqA",
        "x_start": 227,
        "x_end": 231,
        "y_seq": None,
        "y_start": None,
        "y_end": None,
        "ts": 1,
    }
    assert parse_dotplot_go_to_targets(payload) == [("seqA", 227, 231)]


def test_previous_flexidot_y_axis_payload() -> None:
    payload = {
        "x_seq": None,
        "x_start": None,
        "x_end": None,
        "y_seq": "seqA",
        "y_start": 229,
        "y_end": 233,
        "ts": 2,
    }
    assert parse_dotplot_go_to_targets(payload) == [("seqA", 229, 233)]


def test_jiangzhao_legacy_x_only_payload() -> None:
    payload = {
        "x_seq": "seqA",
        "x_start": 229,
        "x_end": 233,
        "y_seq": None,
        "y_start": 229,
        "y_end": 233,
        "ts": 3,
    }
    assert parse_dotplot_go_to_targets(payload) == [("seqA", 229, 233)]
