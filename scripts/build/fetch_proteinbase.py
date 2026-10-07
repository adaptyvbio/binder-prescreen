#!/usr/bin/env python3
"""Published Proteinbase designs -> FASTA, from the public API. No credentials.

This is the designed-binder prior art: what other people have already published for any
target. A submission matching it closely is an existing design rather than a new one.

Only public designs are reachable here, which is the correct scope for this arm — the
unpublished set is not prior art anyone may be told about.

Usage: fetch_proteinbase.py --out OUT.fasta [--page-size 500] [--base URL]
"""

import argparse
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

DEFAULT_BASE = "https://proteinbase.com/api/proteins"
VALID = set("ACDEFGHIKLMNPQRSTVWYXBZU")


# The API paginates by `offset`, not by `page` — passing `page` silently returns the
# first page every time, which walks forever and yields `limit` unique records. It also
# caps `limit` at 100 regardless of what you ask for.
MAX_LIMIT = 100


def fetch_page(base, offset, limit, retries=4):
    url = f"{base}?{urllib.parse.urlencode({'offset': offset, 'limit': limit})}"
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=120) as fh:
                return json.load(fh)
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
            if attempt == retries - 1:
                raise
            # the API rate-limits a long walk; back off rather than lose the run
            wait = 2 ** attempt
            print(f"  offset {offset} failed ({exc}); retrying in {wait}s", file=sys.stderr)
            time.sleep(wait)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", required=True)
    ap.add_argument("--page-size", type=int, default=MAX_LIMIT)
    ap.add_argument("--base", default=DEFAULT_BASE)
    a = ap.parse_args()

    limit = min(a.page_size, MAX_LIMIT)
    seen: dict[str, str] = {}
    offset, total, skipped = 0, None, 0
    while True:
        payload = fetch_page(a.base, offset, limit)
        rows = payload.get("data") or []
        if total is None:
            total = payload.get("total")
            print(f"proteinbase: {total} designs to walk")
        if not rows:
            break
        for row in rows:
            seq = (row.get("sequence") or "").strip().upper()
            # the API returns only public designs, but assert it rather than assume:
            # a private sequence in this arm would be reported to a competitor as
            # "already known", which is the one thing this corpus must never do.
            if row.get("visibility") != "public":
                skipped += 1
                continue
            if not seq or set(seq) - VALID:
                skipped += 1
                continue
            name = row.get("slug") or row.get("name") or row.get("id")
            seen.setdefault(seq, f"proteinbase_{name}")
        offset += len(rows)
        print(f"\r  {offset}/{total} walked, {len(seen)} unique", end="", file=sys.stderr)
        if total is not None and offset >= total:
            break
        time.sleep(0.2)
    print(file=sys.stderr)

    with open(a.out, "w") as out:
        for seq, name in seen.items():
            out.write(f">{name}\n{seq}\n")
    print(f"proteinbase_public: {len(seen)} unique sequences -> {a.out} ({skipped} skipped)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
