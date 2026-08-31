"""Verify citations against PubMed and CrossRef. Free, no API key.

Why not ask an LLM: the four fabricated attributions found on 2026-08-31 were
produced by an LLM and would be checked by an LLM with the same failure mode.
PubMed E-utilities and CrossRef are the indexes themselves — a title either
resolves to a record or it does not.

Both are free and keyless. NCBI asks for tool= and email= identifiers and
allows 3 requests/second without a key (10/s with one, which is also free);
CrossRef asks for a mailto in the User-Agent for its polite pool.

Usage:
  uv run python verify_citations.py                 # audit the bibliography
  uv run python verify_citations.py email_lito.md   # audit one file
"""
import json, os, re, sys, time, urllib.parse, urllib.request

BASE = os.path.dirname(os.path.abspath(__file__))
EMAIL = "allanmfx@gmail.com"          # NCBI asks callers to identify themselves
TOOL = "prion-neurotoxicity-audit"
UA = f"{TOOL}/1.0 (mailto:{EMAIL})"

ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
CROSSREF = "https://api.crossref.org/works"


def get(url, params, tries=3):
    q = f"{url}?{urllib.parse.urlencode(params)}"
    for i in range(tries):
        try:
            req = urllib.request.Request(q, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=30) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            if i == tries - 1:
                return {"_error": str(e)}
            time.sleep(1.5 * (i + 1))


def pubmed_lookup(title, author=None, year=None):
    """Search PubMed by title words, optionally narrowed by author and year."""
    term = f'{title}[Title]'
    if author:
        term += f' AND {author}[Author]'
    if year:
        term += f' AND {year}[DP]'
    r = get(ESEARCH, {"db": "pubmed", "term": term, "retmode": "json",
                      "retmax": 5, "tool": TOOL, "email": EMAIL})
    ids = (r or {}).get("esearchresult", {}).get("idlist", [])
    if not ids:
        return None
    s = get(ESUMMARY, {"db": "pubmed", "id": ",".join(ids), "retmode": "json",
                       "tool": TOOL, "email": EMAIL})
    out = []
    for pid in ids:
        d = (s or {}).get("result", {}).get(pid)
        if not d:
            continue
        out.append({
            "pmid": pid,
            "title": d.get("title", "").rstrip("."),
            "authors": [a["name"] for a in d.get("authors", [])][:6],
            "journal": d.get("fulljournalname") or d.get("source"),
            "year": (d.get("pubdate") or "")[:4],
            "volume": d.get("volume"), "pages": d.get("pages"),
            "doi": d.get("elocationid", ""),
        })
    return out


def crossref_lookup(title, rows=3):
    r = get(CROSSREF, {"query.bibliographic": title, "rows": rows,
                       "mailto": EMAIL})
    items = (r or {}).get("message", {}).get("items", [])
    return [{
        "doi": i.get("DOI"),
        "title": (i.get("title") or [""])[0],
        "authors": [f"{a.get('family','')} {a.get('given','')}".strip()
                    for a in (i.get("author") or [])][:6],
        "journal": (i.get("container-title") or [""])[0],
        "year": str((i.get("issued", {}).get("date-parts") or [[None]])[0][0]),
    } for i in items]


# Recover a searchable title and a claimed first author from a bibliography line
TITLE_RE = re.compile(r'^(?:\*\*)?(?P<auth>[A-Z][^.*]{0,120}?)(?:\*\*)?\.\s+'
                      r'(?P<title>[^.*]{15,300}?)\.\s')
YEAR_RE = re.compile(r'\b(19[89]\d|20[0-2]\d)\b')
ID_RE = re.compile(r'(?:PMID[:\s]*(\d+))|(?:PMC(\d+))|(?:doi:\s*(\S+))', re.I)


def parse_entry(line):
    body = line.lstrip('- ').strip()
    m = TITLE_RE.match(body)
    ids = ID_RE.search(body)
    y = YEAR_RE.search(body)
    if m:
        auth = m.group('auth').split(',')[0].strip()
        return {'author': auth, 'title': m.group('title').strip(),
                'year': y.group(1) if y else None,
                'has_id': bool(ids), 'raw': body[:130]}
    # No leading author: treat the sentence before the first period as the title
    parts = body.split('. ')
    if parts and len(parts[0]) > 20:
        return {'author': None, 'title': parts[0].strip(' *'),
                'year': y.group(1) if y else None,
                'has_id': bool(ids), 'raw': body[:130]}
    return None


def norm(s):
    return re.sub(r'[^a-z0-9 ]', ' ', (s or '').lower()).split()


def title_similarity(a, b):
    """Jaccard over word sets. A loose search can return an unrelated paper,
    and reporting that as 'the real citation' would accuse a correct entry of
    being fabricated -- as bad as missing a fabricated one. So a candidate has
    to look like the same title before any mismatch is claimed."""
    wa, wb = set(norm(a)), set(norm(b))
    if not wa or not wb:
        return 0.0
    return len(wa & wb) / len(wa | wb)


MATCH_MIN = 0.65      # below this, we do not claim to have found the paper


def fetch_abstract(pmid):
    q = f"{EFETCH}?" + urllib.parse.urlencode(
        {"db": "pubmed", "id": pmid, "rettype": "abstract", "retmode": "text",
         "tool": TOOL, "email": EMAIL})
    try:
        req = urllib.request.Request(q, headers={"User-Agent": UA})
        with urllib.request.urlopen(req, timeout=30) as r:
            return r.read().decode(errors='replace').strip()
    except Exception as e:
        return f"(abstract unavailable: {e})"


def check(entry, want_abstract=False):
    candidates, src = [], 'pubmed'
    for kwargs in ({'year': entry['year']}, {}):
        hits = pubmed_lookup(entry['title'], **kwargs)
        if hits:
            candidates = hits
            break
    if not candidates:
        candidates = crossref_lookup(entry['title'])
        src = 'crossref'
    if not candidates:
        return {'status': 'NOT FOUND', 'source': None, 'best': None,
                'similarity': 0.0}

    scored = sorted(((title_similarity(entry['title'], c['title']), c)
                     for c in candidates), key=lambda x: -x[0])
    sim, best = scored[0]
    if sim < MATCH_MIN:
        # We found records, but none is convincingly this paper. Say so
        # instead of naming an unrelated one as the correction.
        return {'status': 'NO CONFIDENT MATCH', 'source': src, 'best': best,
                'similarity': sim,
                'note': 'search returned records but none matched the title closely'}

    mismatch = []
    if entry['author']:
        surname = entry['author'].split()[0].lower()
        if not any(surname in a.lower() for a in best['authors']):
            mismatch.append(f"author '{entry['author']}' not in {best['authors'][:3]}")
    if entry['year'] and best.get('year') and entry['year'] != best['year']:
        mismatch.append(f"year {entry['year']} vs {best['year']}")
    if want_abstract and best.get('pmid'):
        best = {**best, 'abstract': fetch_abstract(best['pmid'])}
    return {'status': 'MISMATCH' if mismatch else 'OK', 'source': src,
            'best': best, 'mismatch': mismatch, 'similarity': sim}


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else 'docs/references/BIBLIOGRAPHY.md'
    limit = int(os.environ.get('LIMIT', '0'))
    lines = [l for l in open(os.path.join(BASE, path))
             if l.startswith('- ') and not l.startswith('- **[')]
    entries = [e for e in (parse_entry(l) for l in lines) if e]
    if limit:
        entries = entries[:limit]
    print(f"{path}: {len(entries)} parseable entries\n")

    want_abs = os.environ.get('ABSTRACTS') == '1'
    results = {'OK': [], 'MISMATCH': [], 'NOT FOUND': [], 'NO CONFIDENT MATCH': []}
    for i, e in enumerate(entries, 1):
        r = check(e, want_abstract=want_abs)
        results[r['status']].append({'entry': e, 'result': r})
        mark = {'OK': 'ok  ', 'MISMATCH': 'DIFF', 'NOT FOUND': 'MISS',
                'NO CONFIDENT MATCH': '????'}[r['status']]
        print(f"[{i:3d}/{len(entries)}] {mark} {e['title'][:64]}")
        if r['status'] == 'MISMATCH':
            for m in r['mismatch']:
                print(f"            -> {m}")
            b = r['best']
            print(f"            -> real: {', '.join(b['authors'][:3])}. "
                  f"{b.get('journal')} {b.get('year')}. "
                  f"PMID {b.get('pmid') or '-'} {b.get('doi') or ''}")
        time.sleep(0.35)          # stay under the 3 req/s keyless limit

    print(f"\n{'='*70}")
    print(f"OK                 {len(results['OK'])}")
    print(f"MISMATCH           {len(results['MISMATCH'])}   <- wrong author or year")
    print(f"NO CONFIDENT MATCH {len(results['NO CONFIDENT MATCH'])}   <- needs a human look, NOT an accusation")
    print(f"NOT FOUND          {len(results['NOT FOUND'])}   <- nothing indexed under this title")
    out = os.path.join(BASE, 'docs/references/VERIFICATION.json')
    with open(out, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nWrote {out}")


if __name__ == '__main__':
    main()
