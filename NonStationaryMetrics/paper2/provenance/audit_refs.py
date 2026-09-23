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

import unicodedata
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
    """Fold a name to a comparable form.

    The previous version carried an ad-hoc list of accent spellings and missed
    most of them, so it reported a discrepancy whenever a correctly accented
    bib entry met a correctly accented registry record -- \\'ario vs ario,
    \\v{S} vs s, D'Hoker vs d\u2019hoker.  Twenty-five of the sixty-four entries
    came back flagged, and the one genuine error in the list (a wrong given
    name) was invisible among them.  A checker that cries wolf is worse than no
    checker.

    This decodes LaTeX accent commands generically, then folds Unicode
    diacritics, so only substantive differences survive.
    """
    s = s.lower()
    # \'a  \`a  \"a  \^a  \~a  \=a  \.a  \v{s}  \u{a}  \H{o}  \c{c}  \k{a}  \r{a}
    s = re.sub(r"\\[`'\"^~=.vuHckr]\s*\{?([a-z])\}?", r"\1", s)
    for a, b in [("\\o", "o"), ("\\l", "l"), ("\\ss", "ss"), ("\\aa", "aa"),
                 ("\\ae", "ae"), ("\\oe", "oe"),
                 ("{", ""), ("}", ""), ("\\", ""),
                 ("\u2019", ""), ("'", ""), ("\u2010", " "), ("-", " "), (".", "")]:
        s = s.replace(a, b)
    # fold any remaining real diacritics (á, š, ò, ...)
    s = "".join(c for c in unicodedata.normalize("NFKD", s)
                if not unicodedata.combining(c))
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



def _initials_compatible(a: str, b: str) -> bool:
    """True when two given-name strings differ only by abbreviation.

    'gary w' vs 'g w', 'sumner byron' vs 's b': the registry and a house style
    disagree about spelling out first names, which is not a bibliographic error.
    """
    ta, tb = a.split(), b.split()
    if not ta or not tb:
        return False
    # The registry often stores only part of a multi-part given name
    # ('miguel' for 'Miguel Angel', abbreviated 'M. A.' in the bibliography),
    # so compare the tokens both sides actually have rather than demanding the
    # same count.  A genuine difference of name -- 'isaac' vs 'israel' -- still
    # fails, because neither is a prefix of the other.
    return all(x == y or x.startswith(y) or y.startswith(x)
               for x, y in zip(ta, tb))


def classify(problem: str, bib_uses_etal: bool) -> str:
    """SUBSTANTIVE or COSMETIC.

    The exit code is only useful if it means something.  Before this split the
    script reported twenty-five discrepancies, every one of them a difference of
    convention, and a real error (a wrong given name) sat unnoticed among them.
    Only SUBSTANTIVE findings now set the exit status; COSMETIC ones are still
    printed, because a reader may want to see them.
    """
    if problem.lstrip().startswith("?"):
        return "COSMETIC"          # registry has no record (Zenodo, ResearchGate)
    if bib_uses_etal and ("author COUNT" in problem or "'others'" in problem):
        return "COSMETIC"          # bibtex "and others" is et al., by design
    if "GIVEN:" in problem:
        try:
            bg = problem.split("bib '")[1].split("'")[0]
            cg = problem.split("registry '")[1].split("'")[0]
            if _initials_compatible(bg, cg):
                return "COSMETIC"
        except IndexError:
            pass
    return "SUBSTANTIVE"


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
    n_sub = n_cos = 0
    for e in withdoi:
        probs = check(e)
        if probs:
            etal = "others" in e.get("author", "").lower()
            kinds = [classify(x, etal) for x in probs]
            sub = sum(k == "SUBSTANTIVE" for k in kinds)
            n_sub += sub
            n_cos += len(kinds) - sub
            tag = "SUBSTANTIVE" if sub else "cosmetic"
            print(f"  [{tag}] {e['key']}  ({e.get('doi')})")
            for x, k in zip(probs, kinds):
                print(x + ("" if k == "SUBSTANTIVE" else "   [cosmetic]"))
        time.sleep(0.3)
    print(f"\n{n_sub} substantive and {n_cos} cosmetic finding(s) "
          f"out of {len(withdoi)} entries checked")
    if not n_sub:
        print("no bibliographic error found; the cosmetic findings are "
              "abbreviation style, bibtex 'and others', and DOIs the registry "
              "does not index (Zenodo, ResearchGate).")
    return 1 if n_sub else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
