#include <fstream>
#include <iostream>
#include <string>

#include "TPZGmshReader.h"
#include "TPZVTKGeoMesh.h"
#include "TPZWannGeometryTools.h"

namespace {
bool WriteMesh(TPZGeoMesh& mesh, const std::string& filename) {
    std::ofstream output(filename);
    if (!output) {
        std::cerr << "Could not open " << filename << " for writing\n";
        return false;
    }
    TPZVTKGeoMesh::PrintGMeshVTK(&mesh, output);
    return static_cast<bool>(output);
}
}

int main(int argc, char* argv[]) {
    if (argc > 2) {
        std::cerr << "Usage: " << argv[0] << " [input.json]\n";
        return 1;
    }
    if (argc == 2) {
        std::string input = argv[1];
        if (input.find(".json") == std::string::npos) input += ".json";
        ProblemData data;
        data.ReadJson(input);
        TPZGeoMesh* mesh = TPZWannGeometryTools::ReadMeshFromGmsh(&data);
        const bool written = WriteMesh(*mesh, "reservoirMesh.vtk");
        delete mesh;
        return written ? 0 : 1;
    }

    // The bundled mesh parts already carry their WANN material IDs.
    TPZGmshReader reader;
    TPZGeoMesh nearWell, outer;
    reader.GeometricGmshMesh(std::string(MESHESDIRNEW) + "/nearWell.msh", &nearWell);
    reader.GeometricGmshMesh(std::string(MESHESDIRNEW) + "/reservoir.msh", &outer);
    return WriteMesh(nearWell, "nearWellMesh.vtk") && WriteMesh(outer, "outerMesh.vtk") ? 0 : 1;
}
