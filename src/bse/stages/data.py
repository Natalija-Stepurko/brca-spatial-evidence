"""Stage `data`: download every input with sha256 provenance, copy study 2's tables at a pinned commit,
unpack the archives and write the section inventory (DESIGN §3, §6).

Downloads are resumable and skipped when the file is present with the recorded sha256.
"""
import sys
import tarfile
import time
import zipfile
from pathlib import Path

import requests

from bse import config as C
from bse import provenance as P
from bse.sources import SOURCES, STUDY2, STUDY2_FILES

RAW = C.DATA / "raw"
UNPACKED = C.DATA / "unpacked"
OUT = C.RESULTS / "data"
STUDY2_COMMIT = "b5cba3226a2f44124c3a11f6a311e74259aa70da"


def fetch(url: str, dest: Path, size_mb: float, retries=4):
    part = dest.with_suffix(dest.suffix + ".part")
    have = part.stat().st_size if part.exists() else 0
    for attempt in range(retries):
        try:
            headers = {"Range": f"bytes={have}-"} if have else {}
            with requests.get(url, stream=True, timeout=180, headers=headers, allow_redirects=True) as r:
                if r.status_code == 416:
                    break
                r.raise_for_status()
                mode = "ab" if r.status_code == 206 else "wb"
                if mode == "wb": have = 0
                with part.open(mode) as f:
                    last = have
                    for chunk in r.iter_content(1 << 20):
                        f.write(chunk); have += len(chunk)
                        if have - last > 200 << 20:
                            last = have
                            print(f"    {dest.name}: {have / 1e6:,.0f} / ~{size_mb:,.0f} MB", flush=True)
            break
        except (requests.RequestException, OSError) as e:
            if attempt == retries - 1: raise
            print(f"    retry {attempt + 1} after {type(e).__name__}", flush=True)
            time.sleep(15 * (attempt + 1))
    part.rename(dest)


def recorded() -> dict:
    out = {}
    if P.MANIFEST.exists():
        for line in P.MANIFEST.read_text().splitlines():
            digest, rel, *_ = line.split("  ")
            out[rel] = digest
    return out


def download_all():
    RAW.mkdir(parents=True, exist_ok=True)
    rec = recorded()
    for s in SOURCES:
        dest = RAW / s.local
        rel = str(dest.relative_to(C.ROOT))
        if dest.exists() and rec.get(rel) == P.sha256(dest):
            print(f"  ok      {s.local}"); continue
        print(f"  fetch   {s.local}  ({s.release}, ~{s.size_mb:,.0f} MB)", flush=True)
        fetch(s.url, dest, s.size_mb)
        print(f"  sha256  {P.record_file(dest, s.url, s.release)[:16]}…  {s.local}")


def copy_study2():
    """Study 2's result tables at a pinned commit, through the raw GitHub endpoint."""
    d = C.DATA / "study2"
    d.mkdir(parents=True, exist_ok=True)
    for rel in STUDY2_FILES:
        url = f"{STUDY2}/{STUDY2_COMMIT}/{rel}"
        dest = d / rel.replace("/", "__")
        if not dest.exists():
            r = requests.get(url, timeout=60); r.raise_for_status()
            dest.write_bytes(r.content)
        P.record_file(dest, url, f"brca-target-evidence @ {STUDY2_COMMIT[:12]}")
    print(f"  study 2 tables: {len(STUDY2_FILES)} files at {STUDY2_COMMIT[:12]}")


def unpack():
    UNPACKED.mkdir(parents=True, exist_ok=True)
    for name, sub in (("wu2021_filtered_count_matrices.tar.gz", "wu"), ("wu2021_spatial.tar.gz", "wu"),
                      ("wu2021_metadata.tar.gz", "wu"), ("janesick_vis_spatial.tar.gz", "janesick")):
        target = UNPACKED / sub
        marker = target / f".{name}.done"
        if marker.exists(): continue
        target.mkdir(parents=True, exist_ok=True)
        with tarfile.open(RAW / name) as t:
            t.extractall(target, filter="data")
        marker.touch(); print(f"  unpacked {name}")
    for name, sub in (("li2025_spaceranger_output.zip", "li"), ("li2025_images.zip", "li")):
        target = UNPACKED / sub
        marker = target / f".{name}.done"
        if marker.exists(): continue
        target.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(RAW / name) as z:
            z.extractall(target)
        marker.touch(); print(f"  unpacked {name}")


def run():
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    print("downloading", flush=True)
    download_all()
    copy_study2()
    print("unpacking", flush=True)
    unpack()
    from bse.stages.inventory import write_inventory  # separate module: it needs the unpacked layout
    inv = write_inventory()
    P.log_run("data", {"study2_commit": STUDY2_COMMIT, "n_sources": len(SOURCES)},
              [str(p.relative_to(C.ROOT)) for p in sorted(OUT.glob("*"))], time.time() - t0)
    print(f"done in {time.time() - t0:,.0f}s; {len(inv)} sections in the inventory")


if __name__ == "__main__":
    sys.exit(run())
