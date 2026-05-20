#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

read -r -p "Select reservoir geometry (Box/Pill): " reservoir_geometry

case "$reservoir_geometry" in
	Box|box)
		reservoir_geo="reservoirBox.geo"
		;;
	Pill|pill)
		reservoir_geo="reservoirPill.geo"
		;;
	*)
		echo "Invalid reservoir geometry: $reservoir_geometry" >&2
		exit 1
		;;
esac

gmsh "$reservoir_geo" -
gmsh nearWell.geo -
gmsh combineMeshes.geo
