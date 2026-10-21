#!/usr/bin/env python3
"""Capture a reproducible inventory for a long-read RNA workflow run."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def command_output(command: list[str], cwd: Path | None = None) -> dict[str, object]:
    try:
        result = subprocess.run(
            command,
            cwd=cwd,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=30,
            check=False,
        )
        return {
            "command": command,
            "returncode": result.returncode,
            "output": result.stdout.strip(),
        }
    except (OSError, subprocess.TimeoutExpired) as exc:
        return {"command": command, "error": str(exc)}


def files_under(paths: list[Path], excluded: set[Path]) -> list[Path]:
    files: set[Path] = set()
    for path in paths:
        if path.is_file() or path.is_symlink():
            files.add(path)
        elif path.is_dir():
            files.update(candidate for candidate in path.rglob("*") if candidate.is_file())
    return sorted(path for path in files if path.resolve() not in excluded)


def file_record(path: Path, with_hash: bool) -> dict[str, object]:
    resolved = path.resolve()
    stat = resolved.stat()
    record: dict[str, object] = {
        "path": str(path),
        "resolved_path": str(resolved),
        "size": stat.st_size,
        "mtime": datetime.fromtimestamp(stat.st_mtime, timezone.utc).isoformat(),
        "is_symlink": path.is_symlink(),
    }
    if with_hash:
        record["sha256"] = sha256(resolved)
    return record


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--retention", choices=("standard", "full"), required=True)
    parser.add_argument("--scope", choices=("input", "output", "failure"), required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--git-root", type=Path)
    parser.add_argument("--path", action="append", default=[], type=Path)
    parser.add_argument("--tool", action="append", default=[], help="NAME=EXECUTABLE")
    args = parser.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    excluded = {args.output.resolve()}
    records = [
        file_record(path, with_hash=args.retention == "full")
        for path in files_under(args.path, excluded)
    ]
    manifest: dict[str, object] = {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "sample_id": args.sample_id,
        "scope": args.scope,
        "retention": args.retention,
        "hash_policy": "sha256_all_files" if args.retention == "full" else "metadata_only",
        "files": records,
    }

    if args.scope == "input":
        manifest["runtime"] = {
            "platform": platform.platform(),
            "python": sys.version,
            "executable": sys.executable,
            "cwd": os.getcwd(),
        }
        if args.git_root:
            manifest["git"] = {
                "head": command_output(["git", "rev-parse", "HEAD"], args.git_root),
                "status": command_output(["git", "status", "--short"], args.git_root),
            }
        tools = {}
        for item in args.tool:
            name, separator, executable = item.partition("=")
            if not separator or not name or not executable:
                parser.error("--tool must use NAME=EXECUTABLE")
            tools[name] = command_output([executable, "--version"])
        manifest["tools"] = tools
        manifest["python_packages"] = command_output([sys.executable, "-m", "pip", "freeze"])
        manifest["conda_packages"] = command_output(["conda", "list", "--json"])

    args.output.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"provenance_manifest={args.output}")
    print(f"provenance_files={len(records)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
