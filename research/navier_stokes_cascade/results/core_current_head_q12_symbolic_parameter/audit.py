"""Reproducible Q12 symbolic arithmetic and frozen provenance, not a PDE proof."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import parameter as p

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
RESULTS = HERE.parent
DEPENDENCIES = ("core_current_head_q11_parameter_obstruction",
                "core_current_head_q10_pressure_comparison_join")


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, RESULTS / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def closure():
    seen = set()
    edges = sources = 0

    def visit(path):
        nonlocal edges, sources
        if path in seen:
            return
        seen.add(path)
        record = json.loads(path.read_text())
        for key, value in record.get("source_hashes", {}).items():
            assert digest(PROJECT/key) == value, key
            sources += 1
        for key, value in record.get("dependency_evidence_sha256", {}).items():
            child = PROJECT/key
            assert digest(child) == value, key
            edges += 1
            visit(child)

    for name in DEPENDENCIES:
        visit(RESULTS/name/"evidence.json")
    return dict(nodes=len(seen), edges=edges, source_bindings=sources)


def encode(value):
    if isinstance(value, p.Q):
        return str(value)
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [encode(item) for item in value]
    return value


def generate():
    q11 = load("q12_parent_audit", "core_current_head_q11_parameter_obstruction/audit.py")
    assert q11.generate() == json.loads((q11.HERE/"evidence.json").read_text())
    samples = [p.rational_plan(h, m, N) for h, m, N in (
        (p.Q(1, 200), 0, 0), (p.Q(1, 513), 1, 3),
        (p.Q(1, 10**30), 5, 11), (p.Q(1, 10**100), 0, 0))]
    for row in samples:
        assert row["complete_minimum"] > row["target"]+1
        assert min(row["reserves"].values()) > 0
        assert min(row["request"]["base_strict_margins"].values()) > 0
        assert row["request"]["integrations_before_extra"]["A1"] == 1
    return {
        "status": "root-parameter-arithmetic-repaired-cutoff-identification-open",
        "parent": "efb5d0777f22aa31cd78b6fad0dbd54a33636856",
        "symbolic_root_plans": [p.symbolic_root_plan(m, N) for m, N in ((0, 0), (1, 3), (5, 11))],
        "rational_diagnostics_not_root_instantiations": encode(samples),
        "closure": closure(),
        "dependency_evidence_sha256": {
            f"results/{name}/evidence.json": digest(RESULTS/name/"evidence.json")
            for name in DEPENDENCIES},
        "source_hashes": {str(path.relative_to(PROJECT)): digest(path)
                          for path in (HERE/"README.md", HERE/"parameter.py", HERE/"audit.py", HERE/"test_audit.py")},
        "scope": {"root_evaluated_numerically": False, "historical_evidence_modified": False,
                  "same_field_join_restored": False, "cutoff_sequence_reselected": False,
                  "endpoint_bound_proved": False, "unconditional_residual_or_blowup": False},
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--verify-record", type=Path)
    args = parser.parse_args()
    record = generate()
    if args.verify_record:
        assert record == json.loads(args.verify_record.read_text())
    with args.output.open("x") as out:
        json.dump(record, out, indent=2, sort_keys=True)
        out.write("\n")
    print(json.dumps({"status": record["status"], "closure": record["closure"]}))
