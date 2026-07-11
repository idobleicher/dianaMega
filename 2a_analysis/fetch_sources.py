#!/usr/bin/env python3
"""Fetch full source protein sequences for UniProtKB + UniParc accessions.
Resumable: caches to fasta_cache.json; safe to re-run. Prints progress to stdout."""
import json, time, sys, urllib.request, urllib.error, urllib.parse, os

ROOT = os.path.dirname(os.path.abspath(__file__))
recs = json.load(open(os.path.join(ROOT, "parsed_rows.json")))
CACHE = os.path.join(ROOT, "fasta_cache.json")
cache = json.load(open(CACHE)) if os.path.exists(CACHE) else {}

up_acc = sorted({r["acc"] for r in recs if r["db"] == "UniProtKB"})
uparc_acc = sorted({r["acc"] for r in recs if r["db"] == "UniParc"})

def http_get(url, tries=5):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "2A-motif-analysis/1.0 (research)"})
            with urllib.request.urlopen(req, timeout=90) as resp:
                return resp.read().decode("utf-8", "replace")
        except urllib.error.HTTPError as e:
            if e.code in (429, 500, 502, 503, 504):
                time.sleep(2 * (i + 1)); continue
            if e.code == 400:  # bad batch; signal caller
                return None
            time.sleep(1 + i)
        except Exception:
            time.sleep(2 * (i + 1))
    return None

def parse_fasta(text):
    out = {}; hdr = None; buf = []
    for line in text.splitlines():
        if line.startswith(">"):
            if hdr is not None: out[hdr] = "".join(buf)
            hdr = line[1:]; buf = []
        else:
            buf.append(line.strip())
    if hdr is not None: out[hdr] = "".join(buf)
    return out

def save():
    tmp = CACHE + ".tmp"
    json.dump(cache, open(tmp, "w"))
    os.replace(tmp, CACHE)

# ---- UniProtKB: batch 100 via accessions endpoint ----
todo = [a for a in up_acc if a not in cache]
print(f"[UniProtKB] {len(up_acc)} unique, {len(todo)} to fetch", flush=True)
B = 100
for i in range(0, len(todo), B):
    batch = todo[i:i+B]
    url = "https://rest.uniprot.org/uniprotkb/accessions?accessions=" + ",".join(batch) + "&format=fasta"
    text = http_get(url)
    if text:
        fa = parse_fasta(text)
        for hdr, seq in fa.items():
            # header: db|ACC|NAME ...
            parts = hdr.split("|")
            if len(parts) >= 2:
                cache[parts[1]] = seq
    # mark misses so we don't retry forever
    for a in batch:
        cache.setdefault(a, cache.get(a, ""))
    if (i // B) % 5 == 0:
        save()
        print(f"[UniProtKB] {min(i+B,len(todo))}/{len(todo)}", flush=True)
    time.sleep(0.15)
save()
print(f"[UniProtKB] done", flush=True)

# ---- UniParc: batch 50 via stream OR-query ----
todo = [a for a in uparc_acc if a not in cache]
print(f"[UniParc] {len(uparc_acc)} unique, {len(todo)} to fetch", flush=True)
B = 50
for i in range(0, len(todo), B):
    batch = todo[i:i+B]
    q = " OR ".join("upi:" + a for a in batch)
    url = "https://rest.uniprot.org/uniparc/stream?query=" + urllib.parse.quote(q) + "&format=fasta"
    text = http_get(url)
    if text:
        fa = parse_fasta(text)
        for hdr, seq in fa.items():
            upi = hdr.split()[0]  # ">UPI.... status=active"
            cache[upi] = seq
    for a in batch:
        cache.setdefault(a, cache.get(a, ""))
    if (i // B) % 10 == 0:
        save()
        print(f"[UniParc] {min(i+B,len(todo))}/{len(todo)}  (cache={len(cache)})", flush=True)
    time.sleep(0.15)
save()
got = sum(1 for a in up_acc+uparc_acc if cache.get(a))
print(f"[DONE] fetched non-empty for {got}/{len(up_acc)+len(uparc_acc)} accessions", flush=True)
