#pragma once

#include "TPZWannAlgebraicTransport.h"
class TPZCompMesh;
class TPZMultiphysicsCompMesh;

/** Geometry-to-cell/interface correspondence adapted from TSFDataTransfer.
 * Unlike subFlow, accepts uncondensed WANN mixed reservoir elements and integrates
 * their physical normal flux by quadrature (no low-order coefficient-sum assumption).
 * Does not modify the Darcy solution. Reinitialize after changing meshes or flow.
 */
class TPZWannTransportDataTransfer {
public:
    static void Initialize(TPZMultiphysicsCompMesh* flow, TPZCompMesh* transport,
        const ProblemData& data, TPZWannAlgebraicTransport& algebraic);
    static REAL MaxCellFluxImbalance(const TPZWannAlgebraicTransport& algebraic);
};
