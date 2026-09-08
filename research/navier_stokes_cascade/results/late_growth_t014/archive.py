#!/usr/bin/env python3
"""Pack completed extension states without changing any scientific outcome."""
import hashlib
import json

from extend import ROOT, INITIAL_SHA, driver, save


def verify_parts(folder, index):
    combined, size = hashlib.sha256(), 0
    for part in index["parts"]:
        path = (folder / part["path"]).resolve()
        path.relative_to(folder.resolve())
        data = path.read_bytes()
        if len(data) != part["bytes"] or hashlib.sha256(data).hexdigest() != part["sha256"]:
            raise ValueError("Checkpoint part identity changed")
        combined.update(data)
        size += len(data)
    if combined.hexdigest() != index["compressed_sha256"] or size != index["compressed_bytes"]:
        raise ValueError("Compressed checkpoint stream changed")


def main():
    protocol = json.loads((ROOT / "protocol.json").read_text())
    for grid in protocol["evolution_grids"]:
        folder = ROOT / ("reference%d" % grid)
        manifest = json.loads((folder / "manifest.json").read_text())
        if manifest["status"] not in ("running", "completed-trajectory-awaiting-crosschecks"):
            raise ValueError("Inspect the stopped or failed trajectory before archiving")
        if manifest["initial_state_sha256"] != INITIAL_SHA:
            raise ValueError("Different initial field")
        for checkpoint in manifest["checkpoints"]:
            time = checkpoint["time"]
            if time not in protocol["endpoints"]:
                raise ValueError("Undeclared saved target")
            index_path = folder / ("t%03d_checkpoint.json" % round(time * 1000))
            budget = next(b for b in manifest["budgets"] if b["time"] == time)
            identity = {"grid": grid, "cutoff": (grid - 1) // 3, "time": time,
                        "initial_state_sha256": INITIAL_SHA, "compressed_sha256": checkpoint["sha256"],
                        "uncompressed_sha256": checkpoint["uncompressed_sha256"],
                        "uncompressed_bytes": checkpoint["uncompressed_bytes"]}
            if budget["grid"] != grid or budget["cutoff"] != identity["cutoff"] or budget["checkpoint_sha256"] != identity["uncompressed_sha256"]:
                raise ValueError("Saved field differs from its independent budget")
            if index_path.exists():
                index = json.loads(index_path.read_text())
                if any(index[k] != v for k, v in identity.items()):
                    raise ValueError("Existing archive belongs to another field")
                verify_parts(folder, index)
                print("Verified packed state", grid, time, flush=True)
                continue
            path = folder / checkpoint["path"]
            if driver.sha(path) != checkpoint["sha256"]:
                raise ValueError("Saved compressed state changed")
            packed = path.read_bytes()
            directory = folder / ("t%03d_state" % round(time * 1000))
            directory.mkdir(exist_ok=False)
            index = {**identity, "format": "split-gzip", "checkpoint_format": "NSCONT1",
                     "original_path": checkpoint["path"], "compressed_bytes": len(packed),
                     "trajectory_status_when_packed": manifest["status"], "parts": []}
            chunk_size = 3 * 1024 * 1024
            for number, start in enumerate(range(0, len(packed), chunk_size)):
                part = directory / ("part-%03d.gzchunk" % number)
                with part.open("xb") as stream:
                    stream.write(packed[start:start + chunk_size])
                index["parts"].append({"path": str(part.relative_to(folder)),
                                       "bytes": part.stat().st_size, "sha256": driver.sha(part)})
            verify_parts(folder, index)
            save(index_path, index)
            print("Packed", grid, time, "in", len(index["parts"]), "parts", flush=True)


if __name__ == "__main__":
    main()
