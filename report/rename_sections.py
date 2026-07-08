#!/usr/bin/env python3
"""Rename report section headers to Vivli's required labels."""

from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SRC = REPO / "report" / "draft_v2.md"
BAK = REPO / "report" / "draft_v2.md.bak2"

RENAMES = [
    ("## 1. Background and aim", "## 1. Objectives"),
    ("## 4. Discussion and policy implications", "## 4. Impact of the work"),
]


def main():
    text = SRC.read_text()
    BAK.write_text(text)
    print(f"backed up to {BAK.relative_to(REPO)}")

    for old, new in RENAMES:
        if old in text:
            text = text.replace(old, new, 1)
            print(f"renamed: '{old}' -> '{new}'")
        else:
            print(f"NOT FOUND: '{old}' (already renamed or drift; check manually)")

    SRC.write_text(text)
    print(f"wrote {SRC.relative_to(REPO)}")


if __name__ == "__main__":
    main()
