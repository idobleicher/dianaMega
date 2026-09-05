# ZZ- paper

New project, started 2026-09-05. Nothing analysed here yet — this file is the scaffold, and the
sections below are the ones every other project in this repo ended up needing. Fill them in as the
work happens rather than at the end; the point of the README in `ubr3_PG_reviewer/` is that a
reader can check a number without opening a script.

## What this project is

*One paragraph: the question, and who is asking it.*

## Source data

*Which file, which sheet, how many rows, and where it came from. Name the file that is the source
of truth, and say plainly if anything else in `data/` is kept only for provenance and read by
nothing.*

| File | Rows | What it is |
|---|---:|---|
| | | |

## Definitions

*Every term the figures and tables use, defined once. Cutoffs get a justification, not a round
number.*

| Term | Definition |
|---|---|
| | |

## Figures

*One row per figure: the stem, the script that writes it, and what the encoding is. Say what a
reader should conclude from it, and what they should not.*

| Figure | Script | Encoding |
|---|---|---|
| | | |

## What survives correction

*The honest headline. How many tests, how many hits, how many expected by chance, what survives
multiple-testing correction, and what that leaves the paper able to claim.*

## Reproducing

```bash
# every script that writes something into figures/ or data/, in order
```

## Conventions carried over from `ubr3_PG_reviewer`

- **One loader.** Every figure reads the source through a single module, so no two figures can
  disagree about the same number or pick up a stale copy of the data.
- **600 dpi PNG + vector PDF + SVG** for every figure, with live text in the SVG.
- **Colour and annotation defined in one place**, imported by the figure scripts, so a set of
  figures reads as one system.
- **The caption carries the caveat.** If a result does not survive correction, the figure says so
  rather than the reader having to find out elsewhere.
