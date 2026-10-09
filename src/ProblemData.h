#ifndef PROBLEMDATA_H
#define PROBLEMDATA_H

#include <iostream>
#include <unordered_map>
#include "json.hpp"
#include <pzfmatrix.h>
#include <pzvec.h>
#include <string>
#include "dirs_config.h"

// declaration of simulation data class.
// all the data herein used are storaged in a .json file. It can be called and storaged using ReadJson

class ProblemData
{
  struct BoundaryData
  {
    int matid = -1; // bc material ID
    int type;       // bc type 0: direct, 1: neumann
    REAL value;     // bc value
  };

  struct DomainData
  {
    std::string name;                                  // name of the domain
    int matid = -1;                                    // domain material ID
    TPZManVector<REAL, 3> perm;                        // domain permeability
    std::unordered_map<std::string, BoundaryData> BCs; // map containing all the bcs info
  };

  struct WellboreData : public DomainData
  {
    int matidSurf;                      // wellbore surface material ID
    int matidToeSurf;                   // wellbore toe surface material ID
    int matidHeelSurf;                  // wellbore heel surface material ID
    REAL radius;                        // domain radius
    REAL length;                        // domain length
    TPZManVector<REAL, 3> eccentricity; // domain excentricity
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
    // Maybe add how the relative permeability is computed...
  };

  struct MeshData
  {
    std::string file;
    int customRefinement;
    int NumUniformRef;
    int NumDirRef;
    int ToCylindrical;
  };

  struct PostProcData
  {
    std::string wellbore_vtk;
    std::string reservoir_vtk;
    std::string training_data;
    int vtk_resolution;
    int training_resolution;
    int verbosityLevel;
    int nthreads;
  };

  struct NumericsData
  {
    int nthreads;
    int maxIterations;
    int reservoirPorder;
    int wellPorder;
    REAL refPressure;
    REAL pressureScale;
    REAL res_tol;
    REAL corr_tol;
  };

public:
  using json = nlohmann::json; // declaration of json class

  TPZVec<WellboreData> m_Wellbore; // Data for all the wellbores in the simulation

  ReservoirData m_Reservoir; // reservoir data

  TPZVec<FluidData> m_Fluid; // Fluid data for all phases

  MeshData m_Mesh; // mesh data

  PostProcData m_PostProc; // post process data

  NumericsData m_Numerics; // numerics data

  // Auxiliary material IDs for the coupling elements
  // They are initialized to -1, and will be set in the ReadJson function
  int EPressure2DSkin = -1; // material ID for the 2D pressure skin elements
  int EPressureInterface = -1; // material ID for the pressure interface elements
  int EHDivBoundInterface = -1; // material ID for the HDiv boundary interface elements

public:
  ProblemData();

  ~ProblemData();

  void ReadJson(std::string jsonfile);

  void UpdateAuxiliaryMaterialIds();

  void UpdateReservoirBCsForWellbores();

  void ApplyPressureScaling();

  void Print(std::ostream &out = std::cout);

  ProblemData(const ProblemData&) = default;
  ProblemData& operator=(const ProblemData&) = default;
};

#endif