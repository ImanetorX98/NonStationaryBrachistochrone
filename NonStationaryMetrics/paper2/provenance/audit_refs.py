#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
audit_refs.py -- check the bibliography against the DOI registry.

Written after the CQG revision found a wrong page range in a reference we had
been citing for months, and after a referee's own citation of the same paper
turned out to name the wrong author. Neither is catchable by reading the
manuscript; both are catchable here.

For every entry with a DOI, this queries CrossRef -- publisher-deposited metadata
-- and compares author family names, given names, year, volume and page range
against what the .bib says. It reports MISMATCH lines only; silence means the
entry agrees with the registry.

What it does NOT check is whether the cited work says what the manuscript
attributes to it. That needs reading, and is tracked separately.

    python3 audit_refs.py                 # all entries with a DOI
    python3 audit_refs.py --group 3       # the third group of five
    python3 audit_refs.py --key Perlick1991
"""
from __future__ import annotations
import json, re, subprocess, sys, time
from pathlib import Path

HERE = Path(__file__).resolve().parent
BIB = HERE.parent.parent / "paper" / "refs.bib"


def parse_bib(text: str) -> list[dict]:
    out = []
    for m in re.finditer(r'@(\w+)\s*\{\s*([^,\s]+)\s*,(.*?)(?=\n\s*@|\Z)', text, re.S):
        kind, key, body = m.group(1), m.group(2), m.group(3)
        fields = {}
        for fm in re.finditer(r'(\w+)\s*=\s*(\{(?:[^{}]|\{[^{}]*\})*\}|"[^"]*"|[^,\n]+)',
                              body):
            v = fm.group(2).strip().strip('{}"').strip()
            fields[fm.group(1).lower()] = re.sub(r'\s+', ' ', v)
        out.append({"kind": kind, "key": key, **fields})
    return out


def crossref(doi: str) -> dict | None:
    try:
        r = subprocess.run(
            ["curl", "-s", "--max-time", "25",
             "-H", "User-Agent: refs-audit (mailto:iman.rosignoli@gmail.com)",
             f"https://api.crossref.org/works/{doi}"],
            capture_output=True, text=True, timeout=40)
        return json.loads(r.stdout)["message"]
    except Exception:
        return None


def norm(s: str) -> str:
    s = s.lower()
    for a, b in [("\\'e", "e"), ("\\`e", "e"), ("\\\"o", "o"), ("\\\"u", "u"),
                 ("\\^i", "i"), ("\\c{c}", "c"), ("\\o", "o"), ("{", ""), ("}", ""),
                 ("é", "e"), ("è", "e"), ("ö", "o"), ("ü", "u"), ("à", "a"),
                 ("ò", "o"), ("ù", "u"), ("í", "i"), ("ñ", "n"), ("ç", "c"),
                 ("\\", ""), ("-", " "), ("'", ""), (".", "")]:
        s = s.replace(a, b)
    return " ".join(s.split())


def bib_authors(e: dict) -> list[tuple[str, str]]:
    """[(family, given), ...] from a bibtex author field."""
    raw = e.get("author", "")
    out = []
    for part in re.split(r'\s+and\s+', raw):
        part = part.strip()
        if not part:
            continue
        if "," in part:
            fam, giv = part.split(",", 1)
        else:
            bits = part.split()
            fam, giv = (bits[-1], " ".join(bits[:-1])) if len(bits) > 1 else (part, "")
        out.append((norm(fam), norm(giv)))
    return out


def check(e: dict) -> list[str]:
    doi = e.get("doi", "").strip()
    if not doi:
        return []
    m = crossref(doi)
    if m is None:
        return [f"    ? could not reach the registry for {doi}"]
    p = []
    ca = [(norm(a.get("family", "")), norm(a.get("given", "")))
          for a in m.get("author", [])]
    ba = bib_authors(e)
    if ca:
        if len(ca) != len(ba):
            p.append(f"    ! author COUNT: bib {len(ba)}, registry {len(ca)}")
        for i, (bf, bg) in enumerate(ba):
            if i >= len(ca):
                break
            cf, cg = ca[i]
            if bf != cf:
                p.append(f"    ! author {i+1} FAMILY: bib '{bf}' vs registry '{cf}'")
            elif bg and cg and not (bg[0] == cg[0]):
                p.append(f"    ! author {i+1} GIVEN: bib '{bg}' vs registry '{cg}'")
            elif bg and cg and len(bg) > 2 and len(cg) > 2 and bg != cg:
                p.append(f"    ! author {i+1} GIVEN: bib '{bg}' vs registry '{cg}'")
    cy = str(m.get("issued", {}).get("date-parts", [[""]])[0][0] or "")
    by = e.get("year", "").strip()
    if by and cy and by != cy:
        p.append(f"    ! YEAR: bib {by} vs registry {cy}")
    cv, bv = str(m.get("volume", "")), e.get("volume", "")
    if bv and cv and norm(bv) != norm(cv):
        p.append(f"    ! VOLUME: bib {bv} vs registry {cv}")
    cp, bp = str(m.get("page", "")), e.get("pages", "")
    if bp and cp:
        bn = re.findall(r'\d+', bp)
        cn = re.findall(r'\d+', cp)
        if bn and cn and bn[0] != cn[0]:
            p.append(f"    ! PAGES: bib {bp} vs registry {cp}")
    return p


def main(argv):
    entries = [e for e in parse_bib(BIB.read_text(encoding="utf-8"))]
    withdoi = [e for e in entries if e.get("doi")]
    if "--key" in argv:
        k = argv[argv.index("--key") + 1]
        withdoi = [e for e in withdoi if e["key"] == k]
    elif "--group" in argv:
        g = int(argv[argv.index("--group") + 1])
        withdoi = withdoi[(g - 1) * 5: g * 5]
    print(f"{len(entries)} entries, {len([e for e in entries if e.get('doi')])} with a DOI; "
          f"checking {len(withdoi)}")
    bad = 0
    for e in withdoi:
        probs = check(e)
        if probs:
            bad += 1
            print(f"  {e['key']}  ({e.get('doi')})")
            for x in probs:
                print(x)
        time.sleep(0.3)
    print(f"\n{bad} entr{'y' if bad==1 else 'ies'} with a discrepancy "
          f"out of {len(withdoi)} checked")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
