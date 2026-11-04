#pragma once

#include "ProblemData.h"
class TPZGeoMesh;
class TPZCompMesh;

/** Reservoir-only adaptation of TSFApproxCreator::BuildTransportCmesh and
 * CreateInterfaceElements. Returns an owned-by-caller P0 mesh on existing geometry.
 * Requires a conforming, unrefined 2D/3D reservoir mesh with tagged boundaries.
 * Adds transport interface geometry; call only once on a given geometric mesh.
 */
class TPZWannTransportMesh {
public:
    static TPZCompMesh* Build(TPZGeoMesh* geometry, const ProblemData& data);
};
