#!/usr/bin/env python3
"""
Find viral proteins (fragments) generated from polyproteins that start with a required
N-terminus (default: Glycine 'G') after cleavage by 3C/3CL protease or other proteases
(user-configurable motifs).

Input: FASTA file containing one or more polyprotein sequences.
Output: CSV table listing fragments whose N-terminus matches the requested prefix and
        was created by a cleavage event.

Motifs are defined in a JSON config file. Each motif provides:
  - name: protease name (e.g., "3C", "3CLpro")
  - pattern: regex pattern to match the cleavage site window
  - cut_after: integer offset (in residues) from match start after which the cut happens
              (i.e., cut position is match_start + cut_after)

Example:
  pattern "QG", cut_after 1  => cut between Q | G, so downstream fragment begins with G.
  pattern "[LFVM]Q[SGA]", cut_after 2 => cut between Q | (S/G/A)
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
from dataclasses import dataclass
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple


@dataclass(frozen=True)
class Motif:
    name: str
    pattern: str
    cut_after: int


@dataclass(frozen=True)
class CleavageEvent:
    protease: str
    match_start_0: int
    match_end_0: int
    cut_pos_0: int  # cut is between cut_pos_0-1 and cut_pos_0, i.e. downstream starts at cut_pos_0
    match_seq: str


def parse_fasta(path: str) -> Iterator[Tuple[str, str]]:
    """
    Minimal FASTA parser. Yields (record_id, sequence).
    record_id is the first token after '>'.
    """
    header: Optional[str] = None
    seq_chunks: List[str] = []

    def flush():
        nonlocal header, seq_chunks
        if header is None:
            return
        rec_id = header.split()[0]
        seq = "".join(seq_chunks).replace(" ", "").replace("\t", "").upper()
        seq = re.sub(r"[^A-Z]", "", seq)
        yield (rec_id, seq)
        header = None
        seq_chunks = []

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                yield from flush()
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        yield from flush()


def load_motifs(config_path: str) -> List[Motif]:
    with open(config_path, "r", encoding="utf-8") as f:
        raw = json.load(f)

    if isinstance(raw, dict) and "motifs" in raw:
        raw_motifs = raw["motifs"]
    else:
        raw_motifs = raw

    motifs: List[Motif] = []
    if not isinstance(raw_motifs, list):
        raise ValueError("Motif config must be a list or an object with key 'motifs' (list).")

    for i, m in enumerate(raw_motifs):
        if not isinstance(m, dict):
            raise ValueError(f"Motif entry #{i} must be an object.")
        name = str(m.get("name", "")).strip()
        pattern = str(m.get("pattern", "")).strip()
        cut_after = m.get("cut_after", None)
        if not name or not pattern or cut_after is None:
            raise ValueError(f"Motif entry #{i} must include name, pattern, cut_after.")
        if not isinstance(cut_after, int) or cut_after < 1:
            raise ValueError(f"Motif entry #{i} cut_after must be integer >= 1.")
        motifs.append(Motif(name=name, pattern=pattern, cut_after=cut_after))

    return motifs


def find_cleavage_events(seq: str, motif: Motif) -> List[CleavageEvent]:
    """
    Find all cleavage events for a given motif in a single sequence.
    cut_pos_0 points to the FIRST residue of the downstream fragment (0-based index).
    """
    events: List[CleavageEvent] = []
    rx = re.compile(motif.pattern)
    for m in rx.finditer(seq):
        s, e = m.start(), m.end()
        cut = s + motif.cut_after
        # cut must be within match and within sequence bounds
        if cut <= s or cut >= e:
            continue
        if cut < 1 or cut > len(seq) - 1:
            continue
        match_seq = seq[s:e]
        events.append(
            CleavageEvent(
                protease=motif.name,
                match_start_0=s,
                match_end_0=e,
                cut_pos_0=cut,
                match_seq=match_seq,
            )
        )
    return events


def digest_sequence(seq: str, cut_positions_0: Sequence[int]) -> List[Tuple[int, int]]:
    """
    Return list of (start_0, end_0_exclusive) segments after cutting at cut_positions_0.
    cut_positions_0 indicates the start index of downstream segments.
    """
    cuts = sorted(set([c for c in cut_positions_0 if 0 < c < len(seq)]))
    bounds = [0] + cuts + [len(seq)]
    segments: List[Tuple[int, int]] = []
    for a, b in zip(bounds, bounds[1:]):
        if a < b:
            segments.append((a, b))
    return segments


def context_around_cut(seq: str, cut_pos_0: int, flank: int = 6) -> str:
    """
    Provide a readable context string around the cut:
      AAAAAB|CCCCCC
    where | marks the cut.
    """
    left_start = max(0, cut_pos_0 - flank)
    right_end = min(len(seq), cut_pos_0 + flank)
    left = seq[left_start:cut_pos_0]
    right = seq[cut_pos_0:right_end]
    return f"{left}|{right}"


def build_rows_for_record(
    rec_id: str,
    seq: str,
    motifs: Sequence[Motif],
    require_created_by_cleavage: bool,
    flank: int,
    required_nterm: str,
) -> List[Dict[str, object]]:
    """
    For each protease motif, digest independently and return rows for fragments whose
    N-terminus starts with required_nterm (e.g., "G" or "GE"). If require_created_by_cleavage,
    exclude the first segment.
    """
    rows: List[Dict[str, object]] = []

    required_nterm = required_nterm.upper()
    if not required_nterm:
        raise ValueError("required_nterm must not be empty")

    for motif in motifs:
        events = find_cleavage_events(seq, motif)
        cut_positions = [ev.cut_pos_0 for ev in events]
        segments = digest_sequence(seq, cut_positions)

        # Map cut_pos -> list of events (there can be multiple matches pointing to same cut)
        by_cut: Dict[int, List[CleavageEvent]] = {}
        for ev in events:
            by_cut.setdefault(ev.cut_pos_0, []).append(ev)

        for idx, (s0, e0) in enumerate(segments, start=1):
            frag = seq[s0:e0]
            if not frag:
                continue
            created_by_cleavage = s0 != 0
            if require_created_by_cleavage and not created_by_cleavage:
                continue
            if not frag.startswith(required_nterm):
                continue

            upstream_events = by_cut.get(s0, [])
            upstream_match = ";".join([ev.match_seq for ev in upstream_events]) if upstream_events else ""
            upstream_match_span = (
                ";".join([f"{ev.match_start_0+1}-{ev.match_end_0}" for ev in upstream_events])
                if upstream_events
                else ""
            )

            rows.append(
                {
                    "polyprotein_id": rec_id,
                    "protease": motif.name,
                    "fragment_index": idx,
                    "fragment_start_1based": s0 + 1,
                    "fragment_end_1based": e0,
                    "fragment_length": e0 - s0,
                    "n_term": frag[0],
                    "n_term_2aa": frag[:2],
                    "n_term_5aa": frag[:5],
                    "c_term": frag[-1],
                    "upstream_cut_pos_1based": (s0 + 1) if created_by_cleavage else "",
                    "upstream_site_context": context_around_cut(seq, s0, flank=flank)
                    if created_by_cleavage
                    else "",
                    "upstream_motif_match": upstream_match,
                    "upstream_motif_match_span_1based": upstream_match_span,
                }
            )

    return rows


def write_csv(path: str, rows: Sequence[Dict[str, object]]) -> None:
    if not rows:
        # still write header for consistency
        header = [
            "polyprotein_id",
            "protease",
            "fragment_index",
            "fragment_start_1based",
            "fragment_end_1based",
            "fragment_length",
            "n_term",
            "n_term_2aa",
            "n_term_5aa",
            "c_term",
            "upstream_cut_pos_1based",
            "upstream_site_context",
            "upstream_motif_match",
            "upstream_motif_match_span_1based",
        ]
        with open(path, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=header)
            w.writeheader()
        return

    header = list(rows[0].keys())
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=header)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def print_preview(rows: Sequence[Dict[str, object]], limit: int = 30) -> None:
    """
    Best-effort console preview. Uses pandas if available, otherwise a small fixed-width table.
    """
    if not rows:
        print("No matching fragments found (given the motifs and input sequences).")
        return

    show = rows[:limit]
    try:
        import pandas as pd  # type: ignore

        df = pd.DataFrame(show)
        with pd.option_context("display.max_columns", 999, "display.width", 180):
            print(df.to_string(index=False))
        if len(rows) > limit:
            print(f"\n... showing first {limit} rows of {len(rows)} total")
        return
    except Exception:
        pass

    cols = [
        "polyprotein_id",
        "protease",
        "fragment_index",
        "fragment_start_1based",
        "fragment_end_1based",
        "fragment_length",
        "n_term_5aa",
        "upstream_site_context",
    ]
    # compute widths
    widths: Dict[str, int] = {}
    for c in cols:
        widths[c] = max(len(c), max(len(str(r.get(c, ""))) for r in show))
        widths[c] = min(widths[c], 48)

    def fmt(c: str, v: object) -> str:
        s = str(v if v is not None else "")
        if len(s) > widths[c]:
            s = s[: widths[c] - 1] + "…"
        return s.ljust(widths[c])

    print(" | ".join([c.ljust(widths[c]) for c in cols]))
    print("-+-".join(["-" * widths[c] for c in cols]))
    for r in show:
        print(" | ".join([fmt(c, r.get(c, "")) for c in cols]))
    if len(rows) > limit:
        print(f"\n... showing first {limit} rows of {len(rows)} total")


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(
        description="Digest polyproteins by protease motifs and find fragments that start with a required N-terminus after cleavage."
    )
    p.add_argument("--fasta", required=True, help="Input FASTA of polyproteins.")
    p.add_argument(
        "--motifs",
        default="protease_motifs.json",
        help="JSON file with protease motifs (default: protease_motifs.json).",
    )
    p.add_argument(
        "--out",
        default="glycine_fragments.csv",
        help="Output CSV path (default: glycine_fragments.csv).",
    )
    p.add_argument(
        "--include_polyprotein_nterm",
        action="store_true",
        help="Also include the first fragment if the polyprotein itself starts with G (not created by cleavage).",
    )
    p.add_argument(
        "--nterm",
        default="G",
        help="Required N-terminus prefix for reported fragments (default: G). Example: --nterm GE",
    )
    p.add_argument(
        "--flank",
        type=int,
        default=6,
        help="Amino acids to show on each side of the cut in context (default: 6).",
    )
    p.add_argument(
        "--preview",
        type=int,
        default=30,
        help="Number of rows to preview in console (default: 30; 0 disables preview).",
    )

    args = p.parse_args(list(argv) if argv is not None else None)

    required_nterm = str(args.nterm or "").strip().upper()
    if not required_nterm or not re.fullmatch(r"[A-Z]+", required_nterm):
        print("ERROR: --nterm must be a non-empty string of letters A-Z (e.g., G or GE).", file=sys.stderr)
        return 2

    if not os.path.exists(args.fasta):
        print(f"ERROR: FASTA not found: {args.fasta}", file=sys.stderr)
        return 2
    if not os.path.exists(args.motifs):
        print(f"ERROR: motif config not found: {args.motifs}", file=sys.stderr)
        return 2

    motifs = load_motifs(args.motifs)
    if not motifs:
        print("ERROR: no motifs loaded.", file=sys.stderr)
        return 2

    all_rows: List[Dict[str, object]] = []
    for rec_id, seq in parse_fasta(args.fasta):
        if not seq:
            continue
        rows = build_rows_for_record(
            rec_id=rec_id,
            seq=seq,
            motifs=motifs,
            require_created_by_cleavage=(not args.include_polyprotein_nterm),
            flank=max(0, int(args.flank)),
            required_nterm=required_nterm,
        )
        all_rows.extend(rows)

    # deterministic ordering
    all_rows.sort(
        key=lambda r: (
            str(r.get("polyprotein_id", "")),
            str(r.get("protease", "")),
            int(r.get("fragment_start_1based", 0) or 0),
            int(r.get("fragment_end_1based", 0) or 0),
        )
    )

    write_csv(args.out, all_rows)
    print(f"Wrote {len(all_rows)} rows to {args.out}")
    if args.preview and args.preview > 0:
        print_preview(all_rows, limit=int(args.preview))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

