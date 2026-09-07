#!/usr/bin/env python3
"""Restore the archived 128-grid checkpoint and verify every recorded hash."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path


def restore(index_path, output):
    index = json.loads(index_path.read_text())
    root = index_path.resolve().parent
    packed = bytearray()
    for part in index["parts"]:
        path = (root / part["path"]).resolve()
        path.relative_to(root)
        data = path.read_bytes()
        if len(data) != part["bytes"] or hashlib.sha256(data).hexdigest() != part["sha256"]:
            raise ValueError("Checkpoint part identity mismatch: " + part["path"])
        packed.extend(data)
    if hashlib.sha256(packed).hexdigest() != index["compressed_sha256"]:
        raise ValueError("Combined checkpoint identity mismatch")
    raw = gzip.decompress(packed)
    if len(raw) != index["uncompressed_bytes"] or hashlib.sha256(raw).hexdigest() != index["uncompressed_sha256"]:
        raise ValueError("Restored checkpoint identity mismatch")
    if raw[:8] != b"NSCONT1\n":
        raise ValueError("Unexpected checkpoint format")
    # Exclusive creation protects an existing checkpoint or other user file.
    with output.open("xb") as destination:
        destination.write(raw)
    return index


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--index", type=Path, default=Path(__file__).with_name("fine128_checkpoint.json"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    metadata = restore(args.index, args.output)
    print("Restored grid", metadata["grid"], "at time", metadata["time"], "to", args.output)
