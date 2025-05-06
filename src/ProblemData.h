#ifndef PROBLEMDATA_H
#define PROBLEMDATA_H

#include <iostream>
#include <unordered_map>
#include "json.hpp"
#include <pzfmatrix.h>
#include <pzvec.h>
#include <string>
#include "input_config.h"

// declaration of simulation data class.
// all the data herein used are storaged in a .json file. It can be called and storaged using ReadJson

class ProblemData
{
  struct BoundaryData
  {
    int matid;             // bc material ID
    int type;              // bc type 0: direct, 1: neumann
    REAL value;            // bc value
  };

  struct DomainData
  {
    std::string name; // name of the domain
    int matid;             // domain material ID
    REAL perm;            // domain permeability
    REAL pOrder;           // polynomial approximation order for flux
    std::unordered_map<std::string, BoundaryData> BCs; // map containing all the bcs info
  };

  struct WellboreData : public DomainData
  {
    REAL radius; // domain radius
    REAL length; // domain length
    REAL excentricity; // domain excentricity
  };

  struct ReservoirData : public DomainData
  {
    REAL height;
    REAL width;
    REAL length;
    REAL porosity;
  };

  struct FluidData
  {
    std::string name;
    REAL viscosity;
    REAL density;
  };

  struct MeshData
  {
    std::string file;
    int NumUniformRef;
    int NumDirRef;
    int Resolution;
    bool ToCylindrical;
  };

public:
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

  using json = nlohmann::json; // declaration of json class

  WellboreData m_Wellbore; // wellbore data
  
  ReservoirData m_Reservoir; // reservoir data
  
  FluidData m_Fluid; // fluid data
  
  MeshData m_Mesh; // mesh data
  
  int m_VerbosityLevel; // verbosity level

public:
  ProblemData();

  ~ProblemData();

  void ReadJson(std::string jsonfile);

  void Print(std::ostream &out = std::cout);
};

#endif