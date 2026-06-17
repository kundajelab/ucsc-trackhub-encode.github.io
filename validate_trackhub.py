#!/usr/bin/env python3
import re
import sys
from collections import Counter
from pathlib import Path
from urllib.parse import urlparse


def main():
    if len(sys.argv) != 2:
        sys.exit("usage: validate_trackhub.py hub_repo_dir")
    repo = Path(sys.argv[1])
    hg = repo / "hg38"
    top = hg / "trackDb.txt"
    errors = []

    includes = []
    for line in top.read_text(encoding="utf-8").splitlines():
        if line.startswith("include "):
            includes.append(line.split(None, 1)[1])
    for name in includes:
        if not (hg / name).exists():
            errors.append(f"missing include: {name}")

    track_ids = []
    urls = []
    for path in sorted(hg.glob("trackDb*.txt")):
        for line in path.read_text(encoding="utf-8").splitlines():
            if line.startswith("track "):
                track_ids.append(line.split(None, 1)[1])
            elif line.startswith("bigDataUrl "):
                urls.append(line.split(None, 1)[1])

    duplicates = [k for k, v in Counter(track_ids).items() if v > 1]
    if duplicates:
        errors.append(f"duplicate track IDs: {duplicates[:10]}")

    bad_urls = [
        u for u in urls
        if not re.search(r"\.(bigWig|bigBed|bw|bb)$", urlparse(u).path, re.I)
    ]
    if bad_urls:
        errors.append(f"non-bigWig/bigBed URLs: {bad_urls[:10]}")

    print(f"includes: {len(includes)}")
    print(f"tracks: {len(track_ids)}")
    print(f"bigDataUrls: {len(urls)}")
    if errors:
        for e in errors:
            print(f"ERROR: {e}")
        sys.exit(1)
    print("validation passed")


if __name__ == "__main__":
    main()
