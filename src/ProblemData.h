#ifndef PROBLEMDATA_H
#define PROBLEMDATA_H

#include <iostream>
#include "json.hpp"
#include <pzfmatrix.h>
#include <pzvec.h>
#include <string>
#include "input_config.h"

// declaration of simulation data class.
// all the data herein used are storaged in a .json file. It can be called and storaged using ReadJson

class ProblemData
{
  enum EMatid
  {
    ENone,
    EDomain,
    EFarField,
    ESurfWellCyl,
    ESurfHeel,
    ESurfToe, // 5
    ECurveWell,
    ECurveHeel,
    ECurveToe,
    ESurfWellCylNonLin,
    ECapRock, // 10
    EPressure2DSkin,
    EPressureInterface,
    EPointHeel,
    EPointToe,
    EHDivBoundInterface
  };

  struct BoundaryData
  {
    std::string name = "none"; // name of the bc
    int matid = 0;             // bc material ID
    int type = 0;              // bc type 0: direct, 1: neumann
    REAL value = 0;            // bc value
  };

  struct DomainData
  {
    std::string name = "none"; // name of the domain
    int matid = 0;             // domain material ID
    REAL perm = 1.;            // domain permeability
    REAL pOrder = 1;           // polynomial approximation order for flux
    TPZManVector<BoundaryData, 5> BCs; // vector containing all the bcs info
  };

  struct WellboreData : public DomainData
  {
    REAL radius = 0.1; // domain radius
    REAL length = 0.1; // domain length
  };

  struct ReservoirData : public DomainData
  {
    REAL height = 0.1;
    REAL width = 0.1;
    REAL length = 0.1;
  };

public:
  using json = nlohmann::json; // declaration of json class

  std::string m_MeshName = "none"; // name of the mesh file

  int m_NumUniformRef = 0; // number of uniform refinements

  int m_NumDirRef = 0; // number of directional refinements

  int m_Dimension = 3; // dimension of the problem

  WellboreData m_Wellbore; // wellbore data
  
  ReservoirData m_Reservoir; // reservoir data

  int fResolution = 0; // VTK mesh resolution

public:
  ProblemData();

  ~ProblemData();

  void ReadJson(std::string jsonfile);

  void Print(std::ostream &out = std::cout);
};

#endif