#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

// NeoPZ includes
#include <iostream>
#include <TPZVTKGeoMesh.h>
#include <TPZGmshReader.h>

#include "ProblemData.h"

// =====================
// Global exact solution
// =====================

ProblemData SimData;

// =============
// Main function
// =============

int main(int argc, char *argv[]) {
  std::string jsonfile = "ozkan1999.json";

  if (argc > 2) {
    std::cout << argv[0] << " being called with too many arguments." << std::endl;
    DebugStop();
  } else if (argc == 2) {
    jsonfile = argv[1];
  }

  if (jsonfile.find(".json") == std::string::npos) {
    jsonfile += ".json";
  }

  // Near-well mesh
  {
    std::string path(std::string(MESHESDIRNEW) + "/" + "nearWell.msh");
    TPZGeoMesh* gmeshN = new TPZGeoMesh();
    TPZGmshReader reader;
    TPZManVector<std::map<std::string, int>, 4> stringtoint(4);
    stringtoint[3]["volume_reservoir"] = SimData.EDomain;
    SimData.m_Reservoir.matid = SimData.EDomain;

    stringtoint[2]["surface_wellbore_cylinder"] = SimData.ESurfWellCyl;
    stringtoint[2]["surface_wellbore_heel"] = SimData.ESurfHeel;
    stringtoint[2]["surface_wellbore_toe"] = SimData.ESurfToe;
    // stringtoint[2]["surface_farfield"] = SimData.EFarField;
    // SetBC(SimData, "surface_farfield", SimData.EFarField);
    stringtoint[2]["surface_cap_rock"] = SimData.ECapRock;
    
    stringtoint[1]["curve_wellbore"] = SimData.ECurveWell;
    stringtoint[1]["curve_heel"] = SimData.ECurveHeel;
    stringtoint[1]["curve_toe"] = SimData.ECurveToe;
    SimData.m_Wellbore.matid = SimData.ECurveWell;

    stringtoint[0]["point_heel"] = SimData.EPointHeel;
    stringtoint[0]["point_toe"] = SimData.EPointToe;
    
    reader.SetDimNamePhysical(stringtoint);
    reader.GeometricGmshMesh(path, gmeshN);

    std::ofstream out("nearWellMesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmeshN, out);
  }

  // Outer mesh
  {
    std::string path(std::string(MESHESDIRNEW) + "/" + "reservoir.msh");
    TPZGeoMesh* gmeshO = new TPZGeoMesh();
    TPZGmshReader reader;
    TPZManVector<std::map<std::string, int>, 4> stringtoint(4);
    stringtoint[3]["volume_reservoir"] = SimData.EDomain;
    SimData.m_Reservoir.matid = SimData.EDomain;

    stringtoint[2]["surface_farfield"] = SimData.EFarField;
    stringtoint[2]["surface_cap_rock"] = SimData.ECapRock;

    reader.SetDimNamePhysical(stringtoint);
    reader.GeometricGmshMesh(path, gmeshO);

    std::ofstream out("outerMesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmeshO, out);
  }
}