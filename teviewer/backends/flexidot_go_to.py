"""Parse Flexidot ``window.__flexidotLastClick`` payloads for TEviewer Go-to."""

from __future__ import annotations

from pathlib import Path


def parse_dotplot_go_to_targets(payload: dict) -> list[tuple[str, int, int]]:
    """
    Return (seq_id, start, end) targets from a Flexidot Go-to click.

    Supports:
    - Previous flexidot (axis-specific: only x_* or only y_* set)
    - Jiangzhao 4-arg buttons (coords stored under x_* for both axes)
    """
    targets: list[tuple[str, int, int]] = []
    seen: set[tuple[str, int, int]] = set()

    for axis in ("x", "y"):
        seq = payload.get(f"{axis}_seq")
        start = payload.get(f"{axis}_start")
        end = payload.get(f"{axis}_end")
        item = _axis_target(seq, start, end)
        if item is None or item in seen:
            continue
        seen.add(item)
        targets.append(item)

    if targets:
        return targets

    # Legacy: x_seq populated but y_seq null (both Go-to buttons wrote x_*).
    item = _axis_target(
        payload.get("x_seq"),
        payload.get("x_start"),
        payload.get("x_end"),
    )
    if item is not None:
        targets.append(item)
    return targets


def _axis_target(seq, start, end) -> tuple[str, int, int] | None:
    if seq in (None, ""):
        return None
    try:
        return str(seq), int(start), int(end)
    except (TypeError, ValueError):
        return None


def fasta_for_flexidot_html(html_path: Path) -> Path | None:
    """Map ``sample_flexidot.html`` / ``sample_flexidot_multi.html`` to ``sample.fa``."""
    name = html_path.name
    for marker in ("_flexidot_multi", "_flexidot"):
        if marker not in name:
            continue
        base = name.split(marker, 1)[0]
        for ext in (".fa", ".fasta", ".fna", ".fas"):
            candidate = html_path.parent / f"{base}{ext}"
            if candidate.is_file():
                return candidate
    return None


def is_flexidot_go_to_console_message(message: str) -> bool:
    return message.startswith("TEVIEWER_FLEXIDOT:")


def payload_from_console_message(message: str) -> dict | None:
    import json

    if not is_flexidot_go_to_console_message(message):
        return None
    raw = message[len("TEVIEWER_FLEXIDOT:") :].strip()
    if not raw:
        return None
    try:
        payload = json.loads(raw)
    except (TypeError, ValueError):
        return None
    return payload if isinstance(payload, dict) else None


FLEXIDOT_GO_TO_BRIDGE_JS = """
(function() {
  if (window.__teviewerFlexidotBridgeInstalled) { return; }
  window.__teviewerFlexidotBridgeInstalled = true;
  window.__teviewerFlexidotLastTs = 0;
  function notify(payload) {
    if (!payload || typeof payload.ts !== 'number') { return; }
    if (payload.ts <= window.__teviewerFlexidotLastTs) { return; }
    window.__teviewerFlexidotLastTs = payload.ts;
    try {
      window.__teviewerFlexidotNotifyJson = JSON.stringify(payload);
    } catch (e) {}
  }
  setInterval(function() {
    if (window.__flexidotLastClick) {
      notify(window.__flexidotLastClick);
    }
  }, 120);
  var plotDiv = document.getElementById('flexidot-dotplot');
  if (plotDiv) {
    plotDiv.on('plotly_click', function() {
      setTimeout(function() {
        var box = document.getElementById('flexidot-dotplot-infobox');
        if (box && box.scrollIntoView) {
          box.scrollIntoView({block: 'nearest', behavior: 'smooth'});
        }
      }, 0);
    });
  }
})();
"""
