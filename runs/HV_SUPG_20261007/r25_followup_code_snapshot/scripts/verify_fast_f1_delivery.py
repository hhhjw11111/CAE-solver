"""Verify the packaged hashes and local HTML links without rerunning the solver."""
from pathlib import Path
from html.parser import HTMLParser
import argparse
import hashlib
import json


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.targets = []

    def handle_starttag(self, tag, attrs):
        self.targets.extend(v for k, v in attrs if k in ("src", "href") and v and not v.startswith("#"))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--delivery-root", type=Path, required=True)
    args = parser.parse_args()
    root = args.delivery_root
    manifest = json.loads((root/"delivery_manifest.json").read_text(encoding="utf-8"))
    bad = []
    for name, record in manifest["files"].items():
        with (root/name).open("rb") as stream:
            if hashlib.file_digest(stream, "sha256").hexdigest() != record["sha256"]:
                bad.append(name)
    links = Links()
    links.feed((root/"图像索引.html").read_text(encoding="utf-8"))
    missing = [name for name in links.targets if not (root/name).exists()]
    print(json.dumps({"files": len(manifest["files"]), "bad_hashes": bad,
                      "index_links": len(links.targets), "missing_targets": missing}))
    if bad or missing:
        raise ValueError("Package integrity or index links failed")


if __name__ == "__main__":
    main()
