#!/usr/bin/env python3
"""Download and checksum the compact STARmap visual-cortex input for SFig11."""

from __future__ import annotations

import argparse
import hashlib
import os
import urllib.request
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_OUTPUT = SCRIPT_DIR / "Panel_A_E_data" / "STARmap_mouse_visual_cortex.zip"
SOURCE_URL = (
    "https://zenodo.org/records/10698912/files/"
    "STARmap_mouse_visual_cortex.zip?download=1"
)
EXPECTED_SHA256 = "c1ffc459cf149b25773411fc6c5ef32845104dd8974802f627171fbf394e083c"
EXPECTED_MD5 = "8f8464ab54504d8854a9fb39cdbcb478"
EXPECTED_BYTES = 471_002


def digest(path: Path, algorithm: str) -> str:
    checksum = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            checksum.update(chunk)
    return checksum.hexdigest()


def verify(path: Path) -> None:
    observed = {
        "bytes": path.stat().st_size,
        "sha256": digest(path, "sha256"),
        "md5": digest(path, "md5"),
    }
    expected = {
        "bytes": EXPECTED_BYTES,
        "sha256": EXPECTED_SHA256,
        "md5": EXPECTED_MD5,
    }
    if observed != expected:
        raise RuntimeError(f"STARmap archive mismatch: {observed}; expected {expected}.")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace an existing archive only after the new copy verifies.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists() and not args.overwrite:
        verify(output)
        print(f"Verified existing file: {output}")
        return

    temporary = output.with_suffix(output.suffix + ".part")
    if temporary.exists():
        temporary.unlink()
    request = urllib.request.Request(
        SOURCE_URL,
        headers={"User-Agent": "SimSpace-reproducibility/1.0"},
    )
    try:
        with urllib.request.urlopen(request, timeout=120) as response:
            with temporary.open("wb") as handle:
                while chunk := response.read(1024 * 1024):
                    handle.write(chunk)
        verify(temporary)
        os.replace(temporary, output)
    finally:
        if temporary.exists():
            temporary.unlink()
    print(f"Downloaded and verified {output} ({EXPECTED_BYTES:,} bytes).")


if __name__ == "__main__":
    main()
