#!/usr/bin/env python3
"""Audit bib entries that carry no DOI, by searching CrossRef on the title.

Complements audit_refs.py, which can only check entries that already have a DOI.
Reports volume / first-page / year disagreements against the best title match,
and prints the DOI found so it can be added.
"""
import json, re, subprocess, sys, time
from pathlib import Path

BIB = Path(__file__).resolve().parent.parent.parent / "paper" / "refs.bib"

def parse(text):
    out = []
    for m in re.finditer(r'@(\w+)\s*\{\s*([^,\s]+)\s*,(.*?)(?=\n\s*@|\Z)', text, re.S):
        f = {}
        for fm in re.finditer(r'(\w+)\s*=\s*(\{(?:[^{}]|\{[^{}]*\})*\}|"[^"]*"|[^,\n]+)', m.group(3)):
            f[fm.group(1).lower()] = re.sub(r'\s+', ' ', fm.group(2).strip().strip('{}"').strip())
        out.append({"kind": m.group(1).lower(), "key": m.group(2), **f})
    return out

def clean(t):
    return re.sub(r'[{}\\$]', '', t)

def search(title, author):
    q = clean(title)
    try:
        r = subprocess.run(["curl", "-s", "--max-time", "25", "-G",
            "https://api.crossref.org/works",
            "--data-urlencode", f"query.bibliographic={q}",
            "--data-urlencode", f"query.author={author}",
            "--data-urlencode", "rows=3"],
            capture_output=True, text=True, timeout=40)
        return json.loads(r.stdout)["message"]["items"]
    except Exception:
        return []

def main():
    entries = [e for e in parse(BIB.read_text(encoding="utf-8"))
               if e["kind"] == "article" and not e.get("doi") and e.get("title")]
    print(f"checking {len(entries)} article entries with no DOI\n")
    flagged = 0
    for e in entries:
        fam = re.split(r'\s*(?:and|,)\s*', e.get("author", ""))[0].strip()
        items = search(e["title"], fam)
        if not items:
            print(f"  {e['key']:26s} no registry match"); continue
        best = items[0]
        bt = clean(e["title"]).lower()[:40]
        ct = (best.get("title") or [""])[0].lower()[:40]
        if bt.split()[0] not in ct and ct.split()[0] not in bt:
            print(f"  {e['key']:26s} title mismatch -> '{ct}'"); flagged += 1; continue
        probs = []
        cv, bv = str(best.get("volume", "")), e.get("volume", "")
        if bv and cv and bv != cv: probs.append(f"VOLUME bib {bv} vs {cv}")
        cy = str(best.get("issued", {}).get("date-parts", [[""]])[0][0] or "")
        by = e.get("year", "")
        if by and cy and by != cy: probs.append(f"YEAR bib {by} vs {cy}")
        cp, bp = str(best.get("page", "")), e.get("pages", "")
        bn, cn = re.findall(r'\d+', bp), re.findall(r'\d+', cp)
        if bn and cn and bn[0] != cn[0]: probs.append(f"PAGES bib {bp} vs {cp}")
        if probs:
            flagged += 1
            print(f"  {e['key']:26s} {'; '.join(probs)}")
            print(f"  {'':26s} doi: {best.get('DOI','?')}")
        time.sleep(0.4)
    print(f"\n{flagged} flagged")

main()
