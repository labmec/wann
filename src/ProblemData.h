/**
 * @file
 * @brief Class and method contracts for ProblemData.
 */

#pragma once

#include <iostream>
#include <map>
#include <optional>
#include <unordered_map>
#include <utility>
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
public:
  enum class RelativePermeability { Linear, Quadratic };

  struct NewtonData
  {
    REAL res_tol;
    REAL corr_tol;
    int maxIterations;
    bool printNewtonProgress;
  };

  struct TransportData
  {
    REAL initialWaterSaturation = 0.;
    /// Exterior water saturation by reservoir BC name or lateral well name.
    std::map<std::string, REAL> boundaryWaterSaturation;
    RelativePermeability relativePermeability = RelativePermeability::Linear;
    std::string waterFluid = "water";
    std::string oilFluid = "oil";
    REAL timeStep = 0.01;
    int numSteps = 1;
  };

private:
  REAL TransportPhaseViscosity(const std::string &name) const;

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
    REAL viscosity = 0.;
    REAL density = 0.;
    REAL physicalViscosity = viscosity; ///< Viscosity before flow pressure scaling.
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
    NewtonData NewtonFlow{1.e-6, 1.e-6, 10, true};
    NewtonData NewtonTransport{1.e-11, 1.e-6, 30, false};
    int reservoirPorder;
    int wellPorder;
    REAL refPressure;
    REAL pressureScale;
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
  /// Optional settings for reservoir transport at frozen total flow.
  std::optional<TransportData> m_Transport;

  // Auxiliary material IDs for the coupling elements
  // They are initialized to -1, and will be set in the ReadJson function
  int EPressure2DSkin = -1; // material ID for the 2D pressure skin elements
  int EPressureInterface = -1; // material ID for the pressure interface elements
  int ETransportInterface = -1; // material ID for internal transport faces
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
   * @brief Assign the auxiliary IDs above the largest configured material ID.
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

  /// Check the reservoir, phase, and time integration settings for transport.
  void ValidateTransport() const;
  /// Resolve exterior water saturations by physical boundary material ID.
  std::map<int, REAL> TransportBoundarySaturations() const;
  /// Return water fractional flow and its derivative with respect to Sw.
  std::pair<REAL, REAL> FractionalFlow(REAL sw) const;
};