"""Validate workflow references, source coverage, and preserved asset hashes."""
from pathlib import Path
import csv
import hashlib
import json


def main():
    root = Path(__file__).resolve().parents[1]
    references = set()
    manifests = sorted((root / "config/workflows").glob("*.json"))
    for path in manifests:
        config = json.loads(path.read_text(encoding="utf-8"))
        assert path.stem == config["id"], path
        targets = set()
        for item in config["files"]:
            source = (root / item["source"]).resolve()
            target = Path(item["target"])
            assert source.is_relative_to(root) and source.is_file(), source
            assert not target.is_absolute() and ".." not in target.parts, target
            assert str(target).lower() not in targets, (path, target)
            targets.add(str(target).lower())
            references.add(source)
    code = set((root / "src").rglob("*.m"))
    assert code <= references, code - references
    checked = set()
    with (root / "docs/file_inventory.csv").open(encoding="utf-8", newline="") as stream:
        for item in csv.DictReader(stream):
            if item["status"] != "input_asset":
                continue
            path = root / item["destination"]
            if path in checked:
                continue
            assert hashlib.sha256(path.read_bytes()).hexdigest() == item["sha256"], path
            checked.add(path)
    print(f"PASS: {len(manifests)} manifests, {len(code)} MATLAB sources, {len(checked)} unchanged assets.")


if __name__ == "__main__":
    main()
