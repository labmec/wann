#!/usr/bin/env python3

"""Run the geoNew Gmsh workflow using geometry parameters from JSON."""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path


PARAM_KEYS = {
    "Lw": ("WellboreData", "length"),
    "Rw": ("WellboreData", "radius"),
    "Hw": ("WellboreData", "height"),
    "ecc": ("WellboreData", "eccentricity", 2),
    "Hr": ("ReservoirData", "height"),
    "Wr": ("ReservoirData", "width"),
    "Lr": ("ReservoirData", "length"),
}

MESH_PARAMS = {
    "h_div": 5,
    "axial_div": 200,
    "radial_div": 8,
    "p_res": 1.6,
    "p_well": 0.3,
    "size_min_res": 15.0,
    "size_max_res": 300.0,
    "dist_min_res": 80.0,
    "dist_max_res": 200.0,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate geoNew meshes using geometry values stored in JSON."
    )
    parser.add_argument(
        "json_file",
        nargs="?",
        default="../input/penmatcha1999.json",
        help="Path to the JSON input file.",
    )
    parser.add_argument(
        "--gmsh",
        default="gmsh",
        help="Gmsh executable to run.",
    )
    return parser.parse_args()


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        raw_text = handle.read()

    cleaned_text = re.sub(r"//.*$", "", raw_text, flags=re.MULTILINE)
    return json.loads(cleaned_text)


def extract_value(data: dict, path: tuple) -> float:
    value = data
    for key in path:
        value = value[key]
    return float(value)


def extract_optional_value(data: dict, path: tuple) -> float | None:
    value = data
    for key in path:
        if not isinstance(value, dict) or key not in value:
            return None
        value = value[key]
    return float(value)


def build_gmsh_args(values: dict[str, float], geo_file: str, batch_mode: bool) -> list[str]:
    args: list[str] = []
    for name, value in values.items():
        args.extend(["-setnumber", name, f"{value:.16g}"])
    args.append(geo_file)
    if batch_mode:
        args.append("-")
    return args


def run_gmsh(gmsh: str, cwd: Path, args: list[str]) -> None:
    subprocess.run([gmsh, *args], cwd=cwd, check=True)


def build_merge_args(values: dict[str, float], geo_file: str, output_file: str) -> list[str]:
    args = build_gmsh_args(values, geo_file, True)
    return args


def resolve_output_path(json_path: Path, output_file: str) -> Path:
    output_path = Path(output_file)
    if output_path.is_absolute():
        return output_path
    return (json_path.parent / output_path).resolve()


def main() -> int:
    args = parse_args()

    script_dir = Path(__file__).resolve().parent
    json_path = Path(args.json_file)
    if not json_path.is_absolute():
        json_path = (script_dir / json_path).resolve()

    if shutil.which(args.gmsh) is None:
        print(f"Gmsh executable not found: {args.gmsh}", file=sys.stderr)
        return 1

    if not json_path.exists():
        print(f"JSON file not found: {json_path}", file=sys.stderr)
        return 1

    json_data = load_json(json_path)

    values = {name: extract_value(json_data, path) for name, path in PARAM_KEYS.items() if name != "Hw"}
    hw = extract_optional_value(json_data, PARAM_KEYS["Hw"])
    values["Hw"] = hw if hw is not None else values["Hr"] / 2.0
    values.update(MESH_PARAMS)

    reservoir_format = json_data["ReservoirData"]["format"]
    if reservoir_format not in {"box", "pill", "ball"}:
        print(f"Unsupported ReservoirData.format value: {reservoir_format}", file=sys.stderr)
        return 1

    output_file = json_data["MeshData"]["file"]
    if reservoir_format == "box":
        reservoir_geo = "reservoirBox.geo"
    elif reservoir_format == "pill":
        reservoir_geo = "reservoirPill.geo"
    else:
        reservoir_geo = "reservoirBall.geo"
    final_mesh_path = resolve_output_path(json_path, output_file)

    run_gmsh(args.gmsh, script_dir, build_gmsh_args(values, reservoir_geo, True))
    run_gmsh(args.gmsh, script_dir, build_gmsh_args(values, "nearWell.geo", True))
    run_gmsh(args.gmsh, script_dir, build_merge_args(values, "combineMeshes.geo", output_file))

    merged_candidates = [
        script_dir / "combined.msh",
        json_path.parent / "combined.msh",
    ]
    for merged_mesh_path in merged_candidates:
        if merged_mesh_path.exists():
            final_mesh_path.parent.mkdir(parents=True, exist_ok=True)
            if final_mesh_path.exists():
                final_mesh_path.unlink()
            merged_mesh_path.replace(final_mesh_path)
            break

    return 0


if __name__ == "__main__":
    raise SystemExit(main())