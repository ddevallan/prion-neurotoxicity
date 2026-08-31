"""Inventory every citation in the repo and rank it by verification risk.

The attribution audit on 2026-08-31 found four wrong attributions in six
checked entries, including an invented author name attached to two unrelated
papers. That rate means an unchecked citation here is unverified, not correct.
This script produces the worklist: what exists, what carries a resolvable
identifier, and which entries appear in material addressed to a patient's
medical team, where a fabricated citation would do real harm.
"""
import re, os, json
from collections import OrderedDict

BASE = os.path.dirname(os.path.abspath(__file__))
BIB = 'docs/references/BIBLIOGRAPHY.md'
OUTBOUND = ['email_lito.md']          # leaves the repo, addressed to clinicians
DOCS = ['README.md', 'docs/MODEL_v5.md', 'docs/THERAPEUTICS.md',
        'docs/CROSS_DISCIPLINARY.md', 'docs/MEMANTINE_DOSSIER.md',
        'docs/GAPS.md', 'docs/EXPERIMENTS.md', 'docs/STRESS_TEST.md']

ID_RE = re.compile(r'(PMID[:\s]*\d+|PMC\d+|doi:\s*\S+|10\.\d{4,}/\S+)', re.I)
YEAR_RE = re.compile(r'\b(19[89]\d|20[0-2]\d)\b')
JOURNAL_RE = re.compile(r'\*([^*]+)\*')
# "Surname AB," or "Surname AB, et al." at the start of an entry
AUTHOR_RE = re.compile(r'^\*{0,2}([A-Z][A-Za-z\'\-]+(?:\s+[A-Z]{1,3})?)[,\s]')

VERIFIED_MARKERS = ['Kloda T', 'Phillis JW', 'do Amaral MJ', 'Yan R, Zhang Y',
                    'Torres MDT', 'Riemer C, Schulz-Schaeffer']


def entries(path):
    """Top-level bullet entries, with their continuation lines attached."""
    out, cur = [], None
    for line in open(os.path.join(BASE, path)):
        if line.startswith('> ') or line.startswith('>|'):
            continue                        # the audit banner, not a citation
        if re.match(r'^- ', line):
            if cur:
                out.append(cur)
            cur = line.rstrip()
        elif cur is not None and re.match(r'^\s+', line) and line.strip():
            cur += ' ' + line.strip()
        elif line.startswith('#') or not line.strip():
            if cur:
                out.append(cur)
            cur = None
    if cur:
        out.append(cur)
    return out


bib = entries(BIB)
print("=" * 78)
print("REFERENCE INVENTORY")
print("=" * 78)
print(f"\nBibliography entries: {len(bib)}")

records = []
for e in bib:
    body = e[2:].strip()
    ids = ID_RE.findall(body)
    m = AUTHOR_RE.match(body)
    rec = {
        'text': body[:150],
        'ids': ids,
        'has_id': bool(ids),
        'year': (YEAR_RE.search(body) or [None])[0] if YEAR_RE.search(body) else None,
        'journal': (JOURNAL_RE.search(body).group(1) if JOURNAL_RE.search(body) else None),
        'named_author': bool(m),
        'verified': any(v in body for v in VERIFIED_MARKERS),
    }
    records.append(rec)

has_id = [r for r in records if r['has_id']]
no_id = [r for r in records if not r['has_id']]
verified = [r for r in records if r['verified']]
named_no_id = [r for r in no_id if r['named_author']]

print(f"  verified against a database : {len(verified)}")
print(f"  carry a PMID/PMC/DOI        : {len(has_id)}")
print(f"  no resolvable identifier    : {len(no_id)}")
print(f"    of those, name an author  : {len(named_no_id)}  <-- highest risk")
print(f"  no journal named            : {sum(1 for r in records if not r['journal'])}")

# --- author name reuse: the signal that caught the fabrication -------------
authors = {}
for r in records:
    m = AUTHOR_RE.match(r['text'])
    if m:
        authors.setdefault(m.group(1), []).append(r)
repeated = {a: rs for a, rs in authors.items() if len(rs) > 1}
print(f"\nAuthor surnames appearing on more than one entry: {len(repeated)}")
for a, rs in sorted(repeated.items(), key=lambda x: -len(x[1])):
    yrs = sorted({r['year'] for r in rs if r['year']})
    span = int(yrs[-1]) - int(yrs[0]) if len(yrs) > 1 else 0
    flag = '  <-- CHECK: wide span, unrelated topics?' if span >= 10 else ''
    print(f"  {a:<22} {len(rs)} entries, years {','.join(yrs) or '?'}{flag}")

# --- what actually leaves the repository -----------------------------------
print(f"\n{'=' * 78}")
print("CITATIONS IN OUTBOUND MATERIAL (highest priority)")
print("=" * 78)
for path in OUTBOUND:
    p = os.path.join(BASE, path)
    if not os.path.exists(p):
        continue
    txt = open(p).read()
    ids = ID_RE.findall(txt)
    years = YEAR_RE.findall(txt)
    journals = JOURNAL_RE.findall(txt)
    # crude author-year detections like "Riemer 2008" or "Sandberg et al."
    cites = re.findall(r'\b([A-Z][a-z]{3,})\s+(?:et al\.?,?\s*)?(19[89]\d|20[0-2]\d)', txt)
    uniq = OrderedDict((f"{a} {y}", None) for a, y in cites)
    print(f"\n{path}:")
    print(f"  resolvable identifiers present : {len(set(ids))}")
    print(f"  author-year mentions           : {len(uniq)}")
    for c in uniq:
        print(f"    - {c}")

with open(os.path.join(BASE, 'docs/references/AUDIT_WORKLIST.json'), 'w') as f:
    json.dump({'total': len(records), 'verified': len(verified),
               'with_id': len(has_id), 'without_id': len(no_id),
               'named_author_no_id': len(named_no_id),
               'repeated_authors': {a: [r['text'][:90] for r in rs]
                                    for a, rs in repeated.items()},
               'unverified': [r['text'] for r in records if not r['verified']]},
              f, indent=2)
print(f"\nWrote docs/references/AUDIT_WORKLIST.json")
