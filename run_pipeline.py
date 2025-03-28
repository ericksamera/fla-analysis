#!/usr/bin/env python3

import argparse
import json
import os
from fla_pipeline.config import GlobalConfig, MarkerConfig
from fla_pipeline.pipeline import run_pipeline

def parse_marker_config(json_path: str) -> list[MarkerConfig]:
    with open(json_path) as f:
        raw = json.load(f)

    if not isinstance(raw, list):
        raise ValueError("Marker config JSON must be a list of marker definitions.")

    return [MarkerConfig(**entry) for entry in raw]

def sanitize_for_json(obj):
    """Recursively sanitize an object for JSON serialization."""
    if isinstance(obj, (int, float, str, bool)) or obj is None:
        return obj
    elif isinstance(obj, (list, tuple)):
        return [sanitize_for_json(i) for i in obj]
    elif isinstance(obj, dict):
        return {k: sanitize_for_json(v) for k, v in obj.items()}
    elif hasattr(obj, "dict"):  # GenotypeResult or similar
        return sanitize_for_json(obj.dict())
    elif hasattr(obj, "__dict__"):  # fallback for dataclasses
        return sanitize_for_json(vars(obj))
    else:
        return str(obj)  # last-resort fallback

def main():
    parser = argparse.ArgumentParser(description="Run FLA diploid pipeline.")
    parser.add_argument("file", type=str, help="Path to the FSA file (.fsa)")
    parser.add_argument("--markers", required=True, help="Path to marker config (JSON)")
    parser.add_argument("--output", help="Optional path to write results as JSON")
    args = parser.parse_args()

    if not os.path.exists(args.file):
        print(f"[ERROR] FSA file not found: {args.file}")
        return

    if not os.path.exists(args.markers):
        print(f"[ERROR] Marker config not found: {args.markers}")
        return

    marker_list = parse_marker_config(args.markers)
    config = GlobalConfig()  # uses defaults

    result = run_pipeline(args.file, config, marker_list)

    sanitized = sanitize_for_json(result)

    if args.output:
        with open(args.output, "w") as out:
            json.dump(sanitized, out, indent=4)
        print(f"[✓] Output written to: {args.output}")
    else:
        print(json.dumps(sanitized, indent=4))

if __name__ == "__main__":
    main()
