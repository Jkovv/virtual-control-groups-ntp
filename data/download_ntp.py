"""Download the public NTP/DTT individual-animal data used in this project.

Nothing in this script is specific to the analysis: it only fetches the files
listed in ``MANIFEST`` (the versioned "DTT Integrated High Value Datasets"
published through CEBS) plus the CEBS Study Catalog, which supplies the
laboratory for each study.

Sources
-------
* DTT Integrated High Value Datasets (landing page, versioned downloads):
  https://cebs.niehs.nih.gov/cebs/paper/16015  (DOI 10.22427/NTP-DATA-110-021-001-000-6)
* DTT Data Collections Guided Search (filterable tables, same data):
  https://cebs-ext.niehs.nih.gov/datasets
* CEBS Study Catalog (laboratory, start year, route, vehicle per study):
  https://cebs-ext.niehs.nih.gov/datasets/search/cebs-study-catalog

The two CSV collections are large (in-life body weights 5.8 GB, histopathology
2.5 GB).  Use ``--skip-large`` to download only the four spreadsheets
(~385 MB); the pipeline then skips the body-weight-trajectory and incidence
endpoints and says so.

Usage
-----
    python data/download_ntp.py              # everything into data/raw/
    python data/download_ntp.py --skip-large
    python data/download_ntp.py --dest /other/dir --no-verify
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import requests

BASE = "https://cebs-ext.niehs.nih.gov/cahs/file/download/datasets/"
CATALOG_PAGE = "https://cebs-ext.niehs.nih.gov/datasets/search/cebs-study-catalog"
CATALOG_API = "https://cebs-ext.niehs.nih.gov/datasets/api/dataset/data/download"
CATALOG_FILE_BASE = "https://cebs-ext.niehs.nih.gov/datasets/"
CATALOG_COLUMNS = [
    "STUDY_UID", "DEPOSITOR_STUDY_NUMBER", "PUBLICATION_UID", "STANDARD_STUDY_TYPE",
    "START_YEAR", "STANDARD_TEST_COMPOUND_NAME", "DTXSID", "DTT_CASRN", "STANDARD_SPECIES",
    "STANDARD_STRAIN", "STANDARD_SEX", "STANDARD_ROUTE", "STANDARD_VEHICLE",
    "STANDARD_POSITIVE_CONTROL", "STUDY_TITLE", "GENERATION", "LABORATORY",
]
USER_AGENT = "virtual-control-groups-ntp/0.1 (+https://github.com/Jkovv/virtual-control-groups-ntp)"


@dataclass(frozen=True)
class RemoteFile:
    name: str           # local file name
    url: str
    sha256: str | None  # pinned when known; None = record only
    large: bool = False


# Versions pinned on 2026-10-21 from the DTT landing page.  Checksums were
# computed on the files as downloaded that day; NTP may republish a version in
# place, in which case ``--no-verify`` lets you proceed and the manifest written
# next to the data records what you actually got.
MANIFEST: list[RemoteFile] = [
    RemoteFile("DTT_Integrated_Data_Page_Definitions.xlsx",
               BASE + "Dataset_And_Column_Definitions/DTT_Integrated_Data_Page_Definitions.xlsx",
               "bba1ac2d403707ade701cea01e00d30c016c7b8f3169a77649b9f37e09341bf0"),
    RemoteFile("202602_Clinical_Chemistry_IAD.xlsx",
               BASE + "Clinical_Chemistry_IAD/202602_Clinical_Chemistry_IAD.xlsx",
               "8dc22338d3c0c2192a728a3dc366773a9752056cb808e33c9d89fd6407db93e9"),
    RemoteFile("202602_Hematology_IAD.xlsx",
               BASE + "Hematology_IAD/202602_Hematology_IAD.xlsx",
               "02ef172929f6c28a5f3e4388ca47f534767698f2eef9dbfbc1c89a0b4ee2f25b"),
    RemoteFile("202602_Organ_Weight_IAD.xlsx",
               BASE + "Organ_Weight_IAD/202602_Organ_Weight_IAD.xlsx",
               "6f938d319a816241cc7cd04e555b9f011ca96db6e331bf99e57410f81d765ca8"),
    RemoteFile("202602_Terminal_Bodyweight_IAD.xlsx",
               BASE + "Terminal_Bodyweight_IAD/202602_Terminal_Bodyweight_IAD.xlsx",
               "4e91376209275f66354d8b9f7d0836453a52887bda9166b046daf98b39cd6c8a"),
    RemoteFile("202507_In_Life_Bodyweight_IAD.csv",
               BASE + "In_Life_Bodyweight_IAD/202507_In_Life_Bodyweight_IAD.csv",
               "a3858d767502367e12fa707524d536a77bad39dfcb5acb081dbb5c03a386cc34", large=True),
    RemoteFile("202602_Histopathology_IAD.csv",
               BASE + "Histopathology_IAD/202602_Histopathology_IAD.csv",
               "3a3aca7096a9efd5a91d335f7b7275cbd7de39f61fecff17cf0e4e01ef746294", large=True),
]


def sha256_of(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while block := fh.read(chunk):
            h.update(block)
    return h.hexdigest()


def download(session: requests.Session, url: str, dest: Path, retries: int = 5) -> None:
    """Stream ``url`` to ``dest`` with HTTP range resume and simple retries."""
    tmp = dest.with_suffix(dest.suffix + ".part")
    for attempt in range(1, retries + 1):
        have = tmp.stat().st_size if tmp.exists() else 0
        headers = {"Range": f"bytes={have}-"} if have else {}
        try:
            with session.get(url, headers=headers, stream=True, timeout=(30, 300)) as r:
                if r.status_code == 416:  # range not satisfiable: already complete
                    break
                r.raise_for_status()
                resuming = bool(have) and r.status_code == 206
                total = int(r.headers.get("Content-Length", 0)) + (have if resuming else 0)
                done = have if resuming else 0
                t0 = time.time()
                with tmp.open("ab" if resuming else "wb") as fh:
                    for block in r.iter_content(chunk_size=1 << 20):
                        fh.write(block)
                        done += len(block)
                        if total and time.time() - t0 > 5:
                            t0 = time.time()
                            print(f"    {done / 1e6:9.1f} / {total / 1e6:.1f} MB", end="\r", flush=True)
            print(" " * 40, end="\r")
            break
        except (requests.RequestException, OSError) as exc:
            print(f"    attempt {attempt} failed: {exc}")
            time.sleep(min(60, 5 * attempt))
    else:
        raise RuntimeError(f"could not download {url}")
    tmp.replace(dest)


def fetch_study_catalog(session: requests.Session, dest: Path) -> str:
    """Request a CSV export of the CEBS Study Catalog through the site's own
    download endpoint (the same call the page's Download button makes)."""
    page = session.get(CATALOG_PAGE, timeout=120)
    page.raise_for_status()
    match = re.search(r'name="_token" value="([^"]+)"', page.text)
    if not match:
        raise RuntimeError("could not find the CSRF token on the study catalog page")
    token = match.group(1)
    form = [("_token", token), ("searchFilters", f"_token={token}&query="),
            ("slug", "cebs-study-catalog"), ("previewMode", "false")]
    form += [("columnList[]", c) for c in CATALOG_COLUMNS]
    r = session.post(CATALOG_API, data=form, timeout=600,
                     headers={"X-CSRF-TOKEN": token, "X-Requested-With": "XMLHttpRequest",
                              "Referer": CATALOG_PAGE})
    r.raise_for_status()
    remote_name = re.sub(r"<[^>]+>", "", r.text).strip()
    if not remote_name.endswith(".csv"):
        raise RuntimeError(f"unexpected response from download endpoint: {r.text[:200]!r}")
    download(session, CATALOG_FILE_BASE + remote_name, dest)
    return remote_name


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dest", type=Path, default=Path(__file__).resolve().parent / "raw",
                    help="target directory (default: data/raw)")
    ap.add_argument("--skip-large", action="store_true",
                    help="skip the multi-GB in-life body weight and histopathology CSVs")
    ap.add_argument("--no-verify", action="store_true", help="do not fail on checksum mismatch")
    ap.add_argument("--force", action="store_true", help="re-download files that already exist")
    args = ap.parse_args(argv)

    dest: Path = args.dest
    dest.mkdir(parents=True, exist_ok=True)
    session = requests.Session()
    session.headers["User-Agent"] = USER_AGENT
    record: dict[str, dict] = {}
    status = 0

    for item in MANIFEST:
        if item.large and args.skip_large:
            print(f"skip   {item.name} (--skip-large)")
            continue
        target = dest / item.name
        if target.exists() and not args.force:
            print(f"exists {item.name}")
        else:
            print(f"get    {item.name}")
            download(session, item.url, target)
        digest = sha256_of(target)
        ok = item.sha256 is None or digest == item.sha256
        record[item.name] = {"url": item.url, "sha256": digest, "bytes": target.stat().st_size,
                             "pinned_sha256": item.sha256, "checksum_ok": ok}
        if not ok:
            print(f"  checksum mismatch for {item.name}: got {digest[:12]}..., pinned {item.sha256[:12]}...")
            if not args.no_verify:
                status = 1

    catalog = dest / "cebs_study_catalog.csv"
    if catalog.exists() and not args.force:
        print("exists cebs_study_catalog.csv")
    else:
        print("get    cebs_study_catalog.csv (via datasets API)")
        remote = fetch_study_catalog(session, catalog)
        record["cebs_study_catalog.csv"] = {"url": CATALOG_FILE_BASE + remote,
                                            "sha256": sha256_of(catalog), "bytes": catalog.stat().st_size,
                                            "note": "unversioned export; file name carries the export date"}

    manifest_path = dest / "download_manifest.json"
    existing = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    existing.update(record)
    existing["_downloaded_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    manifest_path.write_text(json.dumps(existing, indent=2))
    print(f"manifest written to {manifest_path}")
    if status:
        print("one or more checksums did not match the pinned values (use --no-verify to accept)")
    return status


if __name__ == "__main__":
    sys.exit(main())
