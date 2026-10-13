/**
 * @file
 * @brief Class and method contracts for ProblemData.
 */

#pragma once

#include <iostream>
#include <unordered_map>
#include "json.hpp"
#include <pzfmatrix.h>
#include <pzvec.h>
#include <string>
#include "dirs_config.h"

/**
 * @brief Simulation configuration shared by mesh builders, analysis, and exporters.
 * ReadJson() reads a file relative to INPUTDIR, derives coupling material IDs,
 * adds impermeable well end surfaces, and applies pressure scaling in place.
 * Constructing this object alone does not initialize all simulation fields.
 */
class ProblemData
{
  /**
   * @brief Boundary material identifier, type, and stored (possibly scaled) value.
   */
  struct BoundaryData
  {
    int matid = -1; // bc material ID
    int type;       // bc type 0: direct, 1: neumann
    REAL value;     // bc value
  };

  /**
   * @brief Common domain name/material ID, permeability storage, and named boundaries.
   */
  struct DomainData
  {
    std::string name;                                  // name of the domain
    int matid = -1;                                    // domain material ID
    TPZManVector<REAL, 3> perm;                        // domain permeability
    std::unordered_map<std::string, BoundaryData> BCs; // map containing all the bcs info
  };

  /**
   * @brief Well line/surface IDs, geometry, and endpoint boundary conditions.
   * The current JSON reader does not populate the inherited perm vector.
   */
  struct WellboreData : public DomainData
  {
    int matidSurf;                      // wellbore surface material ID
    int matidToeSurf;                   // wellbore toe surface material ID
    int matidHeelSurf;                  // wellbore heel surface material ID
    REAL radius;                        // domain radius
    REAL length;                        // domain length
    TPZManVector<REAL, 3> eccentricity; // domain excentricity
  };

  /**
   * @brief Reservoir permeability, geometry metadata, porosity, and boundaries.
   * The coupled builders use the diagonal perm vector; porosity is stored but
   * not used by the current stationary assembly.
   */
  struct ReservoirData : public DomainData
  {
    // Coupled mesh builders do not consume these reservoir dimensions.
    // TODO: consider removing them from this struct.
    // They shall remain in the json file for the meshing generation process.
    REAL height;
    REAL width;
    REAL length;
    REAL porosity;
  };

  /**
   * @brief Fluid properties; viscosity/density are scaled in place by ReadJson().
   */
  struct FluidData
  {
    std::string name;
    REAL viscosity;
    REAL density;
    // Maybe add how the relative permeability is computed...
  };

  /**
   * @brief Mesh filename and requested processing flags.
   * CreateGeoMesh() currently does not apply the refinement/cylindrical flags.
   */
  struct MeshData
  {
    std::string file;
    int customRefinement;
    int NumUniformRef;
    int NumDirRef; // TODO: current version does not support directional refinement.
    int ToCylindrical;
  };

  /**
   * @brief Output names, sampling resolution, verbosity, and export thread count.
   */
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

  /**
   * @brief Assembly/Newton settings, approximation orders, and pressure scaling.
   */
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

  // Core data structures
  /// Configured wells, indexed from zero.
  TPZVec<WellboreData> m_Wellbore;
  /// Reservoir domain and boundary configuration.
  ReservoirData m_Reservoir;
  /// Configured fluids; current mesh builders use entry zero.
  TPZVec<FluidData> m_Fluid;
  /// Mesh import and requested processing options.
  MeshData m_Mesh;
  /// Export and diagnostic options.
  PostProcData m_PostProc;
  /// Solver, order, and pressure-scaling options.
  NumericsData m_Numerics;

  // Auxiliary material IDs for the coupling elements
  // They are initialized to -1, and will be set in the ReadJson function
  int EPressure2DSkin = -1; // material ID for the 2D pressure skin elements
  int EPressureInterface = -1; // material ID for the pressure interface elements
  int EHDivBoundInterface = -1; // material ID for the HDiv boundary interface elements

public:
  /**
   * @brief Construct an object; configure simulation/material data before use.
   */
  ProblemData();
  /**
   * @brief Destroy this object using its default/base-owned resources.
   */
  ~ProblemData();

  /**
   * @brief Construct a copy using the existing base/member copy semantics.
   */
  ProblemData(const ProblemData&) = default;
  /**
   * @brief Assign configuration using the existing base/member copy semantics.
   */
  ProblemData& operator=(const ProblemData&) = default;

  /**
   * @brief Read configuration, derive coupling IDs/BCs, and scale physical data.
   * @param[in] jsonfile Filename relative to the configured INPUTDIR.
   * WellboreData and FluidData accept one object or an array. Required fields
   * are checked with DebugStop(); JSON parsing/conversion may also throw.
   * NumericsData must exist, although individual numerical values have defaults.
   * @note This mutates stored data, including viscosity, density, and pressure BCs.
   */
  void ReadJson(std::string jsonfile);

  /**
   * @brief Assign the three coupling IDs above the largest configured material ID.
   * Call after domain and boundary IDs have been populated.
   */
  void UpdateAuxiliaryMaterialIds();

  /**
   * @brief Insert/replace zero-Neumann reservoir BCs on each well toe and heel surface.
   */
  void UpdateReservoirBCsForWellbores();

  /**
   * @brief Scale viscosity/density by pressureScale and pressure BCs in place.
   * Type-0 BC values become (value-refPressure)*pressureScale; other BC types
   * are unchanged. This transformation is not idempotent; ReadJson() calls it once.
   */
  void ApplyPressureScaling();

  /**
   * @brief Write the stored configuration or material summary to out.
   */
  void Print(std::ostream &out = std::cout);
};