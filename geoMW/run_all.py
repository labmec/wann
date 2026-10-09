#!/usr/bin/env python3
"""Two parallel wells in a shared reservoir, using geoNew's near-well topology."""
from __future__ import annotations

import argparse
import importlib.util
import math
from pathlib import Path
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parent
# Reuse JSON parsing, parameter defaults, and material ID validation.
spec = importlib.util.spec_from_file_location('geo_new', ROOT.parent / 'geoNew/run_all.py')
geo_new = importlib.util.module_from_spec(spec)
spec.loader.exec_module(geo_new)


def parameters(data):
    wells = data['WellboreData']
    if not isinstance(wells, list) or len(wells) != 2:
        raise ValueError('geoMW requires exactly two wells in WellboreData')
    reservoir = data['ReservoirData']
    if reservoir['format'] not in {'box', 'pill'}:
        raise ValueError('geoMW supports box and pill reservoirs')
    hr, width, length = (float(reservoir[k]) for k in ('height', 'width', 'length'))
    if not all(math.isfinite(v) and v > 0 for v in (hr, width, length)):
        raise ValueError('Reservoir dimensions must be finite and positive')
    values = []
    used = set()
    for i, well in enumerate(wells):
        v = geo_new.build_parameters({**data, 'WellboreData': well})
        if not all(math.isfinite(v[k]) and v[k] > 0 for k in ('Lw', 'Rw')):
            raise ValueError('Well lengths and radii must be finite and positive')
        if 2*v['Rw'] >= hr:
            raise ValueError('Well diameters must be smaller than reservoir height')
        # Ignore the old height/eccentricity convention for this placement example.
        v.update(Hw=hr/2, ecc=0, well_y=(-1 if i == 0 else 1)*width/4)
        ids = {value for key, value in v.items() if key.startswith('id_') and key not in
               {'id_reservoir', 'id_farfield', 'id_cap_rock'}}
        if used & ids:
            raise ValueError('The two wells must use distinct material IDs')
        used.update(ids)
        if math.floor(v['h_div']*v['Lw']/(3*hr)) < 2:
            raise ValueError('Well too short for the near-well axial discretization')
        values.append(v)
    if width <= 2*hr:
        raise ValueError('Reservoir width must exceed twice its height for separate near-well boxes')
    max_lw = max(v['Lw'] for v in values)
    if reservoir['format'] == 'box' and length <= max_lw + hr:
        raise ValueError('Box length must exceed longest well length plus reservoir height')
    return values


def reservoir_geometry(values, shape, size_scale=2.0):
    """Extrude a plane with two holes matching the structured near-well boxes."""
    v = values[0]
    hr, width, lr = v['Hr'], v['Wr'], v['Lr']
    lw = max(w['Lw'] for w in values)
    lines = ['SetFactory("OpenCASCADE");', 'Mesh.MshFileVersion = 2.2;',
             'General.NumThreads = 1;', 'Geometry.Tolerance = 1e-10;',
             'Geometry.MatchMeshTolerance = 1e-8;', 'Mesh.ToleranceReferenceElement = 1e-10;']
    def point(tag, x, y):
        lines.append(f'Point({tag}) = {{{x:.16g}, {y:.16g}, {-hr/2:.16g}, 1}};')
    if shape == 'box':
        for tag, x, y in [(1, -(lr-lw)/2, -width/2), (2, (lr+lw)/2, -width/2),
                          (3, (lr+lw)/2, width/2), (4, -(lr-lw)/2, width/2)]:
            point(tag, x, y)
        for tag, a, b in [(1,1,2),(2,2,3),(3,3,4),(4,4,1)]:
            lines.append(f'Line({tag}) = {{{a},{b}}};')
    else:
        # geoNew's pill is a stadium: length of straight section = longest well.
        for tag, x, y in [(1,0,-width/2),(2,lw,-width/2),(3,lw,width/2),
                          (4,0,width/2),(5,0,0),(6,lw,0)]:
            point(tag, x, y)
        lines += ['Line(1) = {1,2};', 'Circle(2) = {3,6,2};',
                  'Line(3) = {3,4};', 'Circle(4) = {1,5,4};']
    lines.append('Curve Loop(1) = {1,2,3,4};' if shape == 'box' else
                 'Curve Loop(1) = {1,-2,3,-4};')
    hole_curves = []
    for i, w in enumerate(values):
        base = 10 + 10*i
        y = w['well_y']
        for j, x, yy in [(0,-hr/2,y-hr/2),(1,w['Lw']+hr/2,y-hr/2),
                          (2,w['Lw']+hr/2,y+hr/2),(3,-hr/2,y+hr/2)]:
            point(base+j,x,yy)
        for j in range(4):
            lines.append(f'Line({base+j}) = {{{base+j},{base+(j+1)%4}}};')
        curves = ','.join(str(base+j) for j in range(4))
        hole_curves.extend(range(base,base+4))
        lines += [f'Curve Loop({i+2}) = {{{curves}}};',
                  f'Transfinite Curve {{{base},{base+2}}} = {math.floor(w["h_div"]*w["Lw"]/(3*hr))};',
                  f'Transfinite Curve {{{base+1},{base+3}}} = {int(w["h_div"])+1};']
    lines += ['Plane Surface(1) = {1,2,3};',
              'Field[1] = Distance;',
              'Field[1].CurvesList = {' + ','.join(map(str,hole_curves)) + '};',
              'Field[1].Sampling = 100;', 'Field[2] = Threshold;',
              'Field[2].InField = 1;',
              f'Field[2].SizeMin = {size_scale*hr};', f'Field[2].SizeMax = {size_scale*lr/30};',
              f'Field[2].DistMin = {hr};', f'Field[2].DistMax = {width/3};',
              'Background Field = 2;', 'Mesh.MeshSizeFromPoints = 0;',
              'Mesh.MeshSizeFromCurvature = 0;', 'Mesh.MeshSizeExtendFromBoundary = 0;',
              'Mesh.Algorithm = 5;', 'Mesh.RecombinationAlgorithm = 1;',
              'Mesh.RecombineMinimumQuality = 0.01;', 'Mesh.Smoothing = 100;', 'Mesh.RecombineOptimizeTopology = 100;',
              'Recombine Surface {1};',
              f'v[] = Extrude {{0,0,{hr}}} {{ Surface{{1}}; Layers{{{int(v["h_div"])}}}; Recombine; }};',
              f'Physical Surface("surface_farfield", {v["id_farfield"]}) = {{v[2],v[3],v[4],v[5]}};',
              f'Physical Surface("surface_cap_rock", {v["id_cap_rock"]}) = {{1,v[0]}};',
              f'Physical Volume("volume_reservoir", {v["id_reservoir"]}) = {{v[1]}};',
              'Mesh 3;', 'Save "reservoir.msh";']
    return '\n'.join(lines) + '\n'


def combine_meshes(paths, output):
    """Merge ASCII MSH2 files, sharing interface nodes and preserving physical IDs.

    Entity tags are offset between files: identical local entity IDs must not
    conflate different wells. Physical tags deliberately remain unchanged.
    """
    nodes, elements, names, coordinates = [], [], {}, {}
    entity_offset = 0
    for index, path in enumerate(paths):
        text = path.read_text().splitlines()
        def section(name):
            start = text.index('$'+name)+1
            end = text.index('$End'+name)
            return text[start:end]
        for row in section('PhysicalNames')[1:]:
            dim, tag, name = row.split(maxsplit=2)
            # Well group names are unique, while reservoir groups are shared.
            if index:
                name = name[:-1] + f'_{index}"' if int(tag) not in names else name
            names[int(tag)] = (dim, name)
        mapping = {}
        for row in section('Nodes')[1:]:
            tag, *xyz = row.split()
            xyz = tuple(map(float, xyz))
            key = tuple(round(x/1e-8) for x in xyz)
            if key not in coordinates:
                coordinates[key] = len(nodes)+1
                nodes.append(xyz)
            mapping[int(tag)] = coordinates[key]
        max_entity = 0
        for row in section('Elements')[1:]:
            fields = list(map(int,row.split()))
            _, kind, ntags = fields[:3]
            tags = fields[3:3+ntags]
            if ntags >= 2:
                max_entity = max(max_entity,tags[1])
                tags[1] += entity_offset
            conn = [mapping[n] for n in fields[3+ntags:]]
            elements.append([kind,ntags,*tags,*conn])
        entity_offset += max_entity+1
    with output.open('w') as f:
        f.write('$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$PhysicalNames\n')
        f.write(f'{len(names)}\n')
        for tag, (dim,name) in sorted(names.items()):
            f.write(f'{dim} {tag} {name}\n')
        f.write(f'$EndPhysicalNames\n$Nodes\n{len(nodes)}\n')
        for i, xyz in enumerate(nodes,1):
            f.write(f'{i} ' + ' '.join(f'{x:.16g}' for x in xyz)+'\n')
        f.write(f'$EndNodes\n$Elements\n{len(elements)}\n')
        for i, element in enumerate(elements,1):
            f.write(f'{i} ' + ' '.join(map(str,element))+'\n')
        f.write('$EndElements\n')
    return len(nodes), len(elements)


def validate_interfaces(mesh_path, values):
    """Check hex corner Jacobians and require conforming well-box interfaces.

    Check connectivity before MSH4 export. Matching coordinates alone would
    miss unmerged nodes or a quad adjoining multiple smaller faces.
    """
    from collections import defaultdict

    lines = mesh_path.read_text().splitlines()
    start = lines.index('$Nodes') + 2
    nodes = {int(row.split()[0]): tuple(map(float, row.split()[1:]))
             for row in lines[start:lines.index('$EndNodes')]}
    faces = defaultdict(list)
    layouts = {
        5: ((0,1,2,3), (4,5,6,7), (0,1,5,4), (1,2,6,5), (2,3,7,6), (3,0,4,7)),
        6: ((0,1,2), (3,4,5), (0,1,4,3), (1,2,5,4), (2,0,3,5)),
    }
    start = lines.index('$Elements') + 2
    for row in lines[start:lines.index('$EndElements')]:
        fields = list(map(int, row.split()))
        kind, ntags = fields[1:3]
        if fields[3] != values[0]['id_reservoir']:
            continue
        if kind not in layouts:
            raise ValueError(f'Unsupported reservoir element type for interface validation: {kind}')
        conn = fields[3+ntags:]
        if kind == 5:
            # A hex can pass integration-point checks yet have a zero Jacobian
            # at a corner (e.g. a recombined quad with three collinear vertices).
            signs = ((-1,-1,-1), (1,-1,-1), (1,1,-1), (-1,1,-1),
                     (-1,-1,1), (1,-1,1), (1,1,1), (-1,1,1))
            neighbors = ((1,0,3,2,5,4,7,6), (3,2,1,0,7,6,5,4), (4,5,6,7,0,1,2,3))
            for corner, sign in enumerate(signs):
                origin = nodes[conn[corner]]
                columns = [tuple((nodes[conn[neighbors[axis][corner]]][j]-origin[j])/(-2*sign[axis])
                                 for j in range(3)) for axis in range(3)]
                a, b, c = columns
                determinant = (a[0]*(b[1]*c[2]-b[2]*c[1])
                               - a[1]*(b[0]*c[2]-b[2]*c[0])
                               + a[2]*(b[0]*c[1]-b[1]*c[0]))
                scale = math.prod(math.sqrt(sum(x*x for x in column)) for column in columns)
                if scale == 0 or determinant <= 1e-10*scale:
                    raise ValueError(f'Degenerate/inverted reservoir hex {fields[0]} at corner {corner}')
        center = tuple(sum(nodes[n][axis] for n in conn)/len(conn) for axis in range(3))
        for face in layouts[kind]:
            faces[tuple(sorted(conn[i] for i in face))].append(center)

    tolerance = 1e-7
    for index, well in enumerate(values):
        hr, y = well['Hr'], well['well_y']
        bounds = ((-hr/2, well['Lw']+hr/2), (y-hr/2, y+hr/2), (-hr/2, hr/2))
        def inside(center):
            return all(lo < coord < hi for coord, (lo, hi) in zip(center, bounds))
        count = 0
        for face, owners in faces.items():
            coords = [nodes[n] for n in face]
            if not all(all(lo-tolerance <= coord <= hi+tolerance
                           for coord, (lo, hi) in zip(point, bounds)) for point in coords):
                continue
            if not any(all(abs(point[axis]-plane) < tolerance for point in coords)
                       for axis in (0, 1) for plane in bounds[axis]):
                continue
            if len(owners) != 2 or sum(inside(center) for center in owners) != 1:
                raise ValueError(f'Nonconforming well {index} interface face: {face}')
            count += 1
        axial = math.floor(well['h_div']*well['Lw']/(3*hr)) - 1
        height = int(well['h_div'])
        expected = 2*height*(axial+height)
        if count != expected:
            raise ValueError(f'Well {index} interface has {count} shared faces; expected {expected}')
        print(f'Well {index}: {count} conforming interface faces')


def run_gmsh(executable, work, args):
    # Some Gmsh versions return zero even after geometry/meshing errors.
    result = subprocess.run([executable, *args, '-v', '3'], cwd=work,
                            capture_output=True, text=True)
    log = result.stdout + result.stderr
    if result.returncode or 'Error' in log:
        raise ValueError(f'Gmsh failed:\n{log}')
    if log.strip():
        print(log, end='')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('json_file', nargs='?', default='../input/twoWell.json')
    parser.add_argument('--gmsh', default='gmsh')
    parser.add_argument('--output', type=Path, help='Override MeshData.file (relative to current directory)')
    parser.add_argument('--reservoir-size-scale', type=float, default=2.0,
                        help='Outer reservoir mesh size multiplier (default: 2; original: 1)')
    args = parser.parse_args()
    try:
        if not math.isfinite(args.reservoir_size_scale) or args.reservoir_size_scale <= 0:
            raise ValueError('--reservoir-size-scale must be finite and positive')
        if shutil.which(args.gmsh) is None:
            raise ValueError(f'Gmsh executable not found: {args.gmsh}')
        path = Path(args.json_file)
        path = path if path.is_absolute() else (ROOT/path).resolve()
        data = geo_new.load_json(path)
        values = parameters(data)
        output = args.output.resolve() if args.output else geo_new.resolve_output_path(path,data['MeshData']['file'])
        with tempfile.TemporaryDirectory(prefix='geoMW-') as tmp:
            work = Path(tmp)
            for name in ('params.geo','nearWell.geo'):
                shutil.copyfile(ROOT/name, work/name)
            (work/'reservoir.geo').write_text(reservoir_geometry(values,data['ReservoirData']['format'], args.reservoir_size_scale))
            run_gmsh(args.gmsh, work, ['reservoir.geo', '-'])
            meshes = [work/'reservoir.msh']
            for i, v in enumerate(values,1):
                run_gmsh(args.gmsh, work, geo_new.build_gmsh_args(v, 'nearWell.geo', True))
                mesh = work/f'well{i}.msh'
                (work/'nearWell.msh').rename(mesh)
                meshes.append(mesh)
            merged = work/'combined.msh'
            nn, ne = combine_meshes(meshes,merged)
            validate_interfaces(merged, values)
            # NeoPZ's TPZGmshReader requires MSH 4.x. Keep the simple MSH2
            # intermediates for merging, then let Gmsh serialize the final mesh
            # with entity/physical-group metadata in ASCII MSH 4.1. No remeshing.
            (work/'export.geo').write_text(
                'Merge "combined.msh";\nMesh.MshFileVersion = 4.1;\n'
                'Mesh.Binary = 0;\nSave "final.msh";\n')
            run_gmsh(args.gmsh, work, ['export.geo', '-'])
            final_mesh = work/'final.msh'
            with final_mesh.open() as handle:
                if handle.readline().strip() != '$MeshFormat' or handle.readline().strip() != '4.1 0 8':
                    raise ValueError('Gmsh did not export ASCII MSH 4.1')
            output.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(final_mesh,output)
        print(f'Wrote {output}: {nn} nodes, {ne} elements')
        return 0
    except (OSError, KeyError, TypeError, ValueError, subprocess.CalledProcessError) as error:
        parser.exit(1,f'Mesh generation failed: {error}\n')


if __name__ == '__main__':
    raise SystemExit(main())
