#!/usr/bin/env python3
"""Pack completed target states; the final audit verifies full restoration.

A completed .08 checkpoint may be packed while the .10 trajectory continues.
That does not mark the whole trajectory or its scientific gates completed.
Original gzip files remain as ignored local copies until final verification.
"""
import hashlib
import json
from pathlib import Path

from analyze import ROOT, INITIAL_SHA, driver, save


def main():
    folder = ROOT / "reference128"
    manifest = json.loads((folder / "manifest.json").read_text())
    if manifest["status"] not in ("running", "completed-trajectory-awaiting-crosschecks"):
        raise ValueError("Unexpected trajectory status; inspect the failed run first")
    if manifest["initial_state_sha256"] != INITIAL_SHA:
        raise ValueError("Initial field differs")
    for checkpoint in manifest["checkpoints"]:
        path = folder / checkpoint["path"]
        index_path = folder / ("t%03d_checkpoint.json" % round(checkpoint["time"] * 1000))
        if index_path.exists():
            index = json.loads(index_path.read_text())
            if index["compressed_sha256"] != checkpoint["sha256"] or index["uncompressed_sha256"] != checkpoint["uncompressed_sha256"]:
                raise ValueError("Existing archive belongs to a different state")
            combined = hashlib.sha256()
            size = 0
            for part in index["parts"]:
                data = (folder / part["path"]).read_bytes()
                if len(data) != part["bytes"] or hashlib.sha256(data).hexdigest() != part["sha256"]:
                    raise ValueError("Existing archive part changed")
                combined.update(data)
                size += len(data)
            if combined.hexdigest() != checkpoint["sha256"] or size != index["compressed_bytes"]:
                raise ValueError("Existing archive stream changed")
            print("Verified existing packed target", checkpoint["time"], flush=True)
            continue
        if driver.sha(path) != checkpoint["sha256"]:
            raise ValueError("Compressed solver state changed")
        packed = path.read_bytes()
        grid, cutoff, time = manifest["settings"]["grid"], 42, checkpoint["time"]
        budget = next(item for item in manifest["budgets"] if item["time"] == time)
        if grid != 128 or budget["cutoff"] != cutoff or budget["checkpoint_sha256"] != checkpoint["uncompressed_sha256"]:
            raise ValueError("Completed target does not match the independent snapshot")
        directory = folder / ("t%03d_state" % round(time * 1000))
        directory.mkdir(exist_ok=False)
        index = {"format": "split-gzip", "checkpoint_format": "NSCONT1", "grid": grid,
                 "cutoff": cutoff, "time": time, "initial_state_sha256": INITIAL_SHA,
                 "original_path": checkpoint["path"], "compressed_bytes": len(packed),
                 "compressed_sha256": checkpoint["sha256"], "uncompressed_bytes": checkpoint["uncompressed_bytes"],
                 "uncompressed_sha256": checkpoint["uncompressed_sha256"],
                 "trajectory_status_when_packed": manifest["status"], "parts": []}
        chunk_size = 3 * 1024 * 1024
        for number, start in enumerate(range(0, len(packed), chunk_size)):
            part = directory / ("part-%03d.gzchunk" % number)
            data = packed[start:start + chunk_size]
            with part.open("xb") as stream:
                stream.write(data)
            index["parts"].append({"path": str(part.relative_to(folder)), "bytes": len(data), "sha256": driver.sha(part)})
        combined = hashlib.sha256()
        for item in index["parts"]:
            combined.update((folder / item["path"]).read_bytes())
        if combined.hexdigest() != checkpoint["sha256"]:
            raise ValueError("Split archive failed compressed-stream reconstruction")
        save(index_path, index)
        print("Packed", grid, time, checkpoint["uncompressed_bytes"], "restored bytes in", len(index["parts"]), "parts", flush=True)


if __name__ == "__main__":
    main()
