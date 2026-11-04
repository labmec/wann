/**
 * @file
 * @brief JSON input loading, derived coupling IDs/BCs, and pressure scaling.
 * See the matching header for function contracts and docs/ for the architecture.
 */


#include <cmath>
#include <iostream>
#include <string>
#include <fstream>
#include <pzerror.h>

#include "ProblemData.h"

using namespace std;

// constructor
ProblemData::ProblemData()
{
    m_Reservoir.perm.resize(3);
    m_Reservoir.BCs.reserve(3);
}

// deconstructor
ProblemData::~ProblemData() {}

// readjson function. takes a json function as parameter and completes the required simulation data
void ProblemData::ReadJson(std::string file)
{
    std::string path(std::string(INPUTDIR) + "/" + file);
    std::ifstream filejson(path);
    json input = json::parse(filejson, nullptr, true, true); // to ignore comments in json file

    if (input.find("MeshData") == input.end())
        DebugStop();
    json meshdata = input["MeshData"];
    if (meshdata.find("file") == meshdata.end())
        DebugStop();
    m_Mesh.file = meshdata["file"];
    if (meshdata.find("customRefinement") == meshdata.end())
        DebugStop();
    m_Mesh.customRefinement = meshdata["customRefinement"];
    if (meshdata.find("NumUniformRef") == meshdata.end())
        DebugStop();
    m_Mesh.NumUniformRef = meshdata["NumUniformRef"];
    if (meshdata.find("NumDirRef") == meshdata.end())
        DebugStop();
    m_Mesh.NumDirRef = meshdata["NumDirRef"];
    if (meshdata.find("ToCylindrical") == meshdata.end())
        DebugStop();
    m_Mesh.ToCylindrical = meshdata["ToCylindrical"];

    if (input.find("WellboreData") == input.end())
        DebugStop();
    auto readWellbore = [&](const json &wellboreJson) {
        WellboreData wellbore;
        if (wellboreJson.find("name") == wellboreJson.end())
            DebugStop();
        wellbore.name = wellboreJson["name"];
        if (wellboreJson.find("matid") == wellboreJson.end())
            DebugStop();
        wellbore.matid = wellboreJson["matid"];
        if (wellboreJson.find("matidSurf") == wellboreJson.end())
            DebugStop();
        wellbore.matidSurf = wellboreJson["matidSurf"];
        if (wellboreJson.find("matidToeSurf") == wellboreJson.end())
            DebugStop();
        wellbore.matidToeSurf = wellboreJson["matidToeSurf"];
        if (wellboreJson.find("matidHeelSurf") == wellboreJson.end())
            DebugStop();
        wellbore.matidHeelSurf = wellboreJson["matidHeelSurf"];
        if (wellboreJson.find("radius") == wellboreJson.end())
            DebugStop();
        wellbore.radius = wellboreJson["radius"];
        if (wellboreJson.find("length") == wellboreJson.end())
            DebugStop();
        wellbore.length = wellboreJson["length"];
        if (wellboreJson.find("eccentricity") == wellboreJson.end())
            DebugStop();
        const json &eccentricity = wellboreJson["eccentricity"];
        if (!eccentricity.is_array() || eccentricity.size() != 3)
            DebugStop();
        wellbore.eccentricity.Resize(3);
        for (int i = 0; i < 3; i++)
        {
            if (eccentricity[i].is_null())
                DebugStop();
            wellbore.eccentricity[i] = eccentricity[i];
        }
        if (wellboreJson.find("BCs") == wellboreJson.end())
            DebugStop();
        const json &bcs = wellboreJson["BCs"];
        if (!bcs.is_array())
            DebugStop();
        for (int i = 0; i < bcs.size(); i++)
        {
            if (bcs[i].find("name") == bcs[i].end())
                DebugStop();
            std::pair<std::string, BoundaryData> bcpair;
            bcpair.first = bcs[i]["name"];
            if (bcs[i].find("matid") == bcs[i].end())
                DebugStop();
            bcpair.second.matid = bcs[i]["matid"];
            if (bcs[i].find("type") == bcs[i].end())
                DebugStop();
            bcpair.second.type = bcs[i]["type"];
            if (bcs[i].find("value") == bcs[i].end())
                DebugStop();
            bcpair.second.value = bcs[i]["value"];
            wellbore.BCs.insert(bcpair);
        }
        m_Wellbore.push_back(wellbore);
    };
    const json &wellboreData = input["WellboreData"];
    m_Wellbore.Resize(0);
    if (wellboreData.is_array())
    {
        for (const auto &wellboreJson : wellboreData)
            readWellbore(wellboreJson);
    }
    else if (wellboreData.is_object())
    {
        readWellbore(wellboreData);
    }
    else
    {
        DebugStop();
    }

    if (input.find("ReservoirData") == input.end())
        DebugStop();
    json reservoir = input["ReservoirData"];
    if (reservoir.find("name") == reservoir.end())
        DebugStop();
    m_Reservoir.name = reservoir["name"];
    if (reservoir.find("matid") == reservoir.end())
        DebugStop();
    m_Reservoir.matid = reservoir["matid"];
    if (reservoir.find("perm") == reservoir.end())
        DebugStop();
    const json &perm = reservoir["perm"];
    if (perm.is_number())
    {
        m_Reservoir.perm[0] = perm;
        m_Reservoir.perm[1] = perm;
        m_Reservoir.perm[2] = perm;
    }
    else if (perm.is_array())
    {
        if (perm.size() != 3)
            DebugStop();
        for (int i = 0; i < 3; i++)
        {
            if (perm[i].is_null())
                DebugStop();
            m_Reservoir.perm[i] = perm[i];
        }
    }
    else
    {
        DebugStop();
    }
    if (reservoir.find("porosity") == reservoir.end())
        DebugStop();
    m_Reservoir.porosity = reservoir["porosity"];
    if (reservoir.find("height") == reservoir.end())
        DebugStop();
    m_Reservoir.height = reservoir["height"];
    if (reservoir.find("width") == reservoir.end())
        DebugStop();
    m_Reservoir.width = reservoir["width"];
    if (reservoir.find("length") == reservoir.end())
        DebugStop();
    m_Reservoir.length = reservoir["length"];
    if (reservoir.find("BCs") == reservoir.end())
        DebugStop();
    m_Reservoir.BCs.clear();
    json bcsres = reservoir["BCs"];
    for (int i = 0; i < bcsres.size(); i++)
    {
        if (bcsres[i].find("name") == bcsres[i].end())
            DebugStop();
        std::pair<std::string, BoundaryData> bcpair;
        bcpair.first = bcsres[i]["name"];
        if (bcsres[i].find("matid") == bcsres[i].end())
            DebugStop();
        bcpair.second.matid = bcsres[i]["matid"];
        if (bcsres[i].find("type") == bcsres[i].end())
            DebugStop();
        bcpair.second.type = bcsres[i]["type"];
        if (bcsres[i].find("value") == bcsres[i].end())
            DebugStop();
        bcpair.second.value = bcsres[i]["value"];
        m_Reservoir.BCs.insert(bcpair);
    }

    if (input.find("FluidData") == input.end())
        DebugStop();
    auto readFluid = [&](const json &fluidJson) {
        if (!fluidJson.contains("name") || !fluidJson.contains("viscosity") ||
            !fluidJson.contains("density"))
            DebugStop();
        m_Fluid.push_back({fluidJson.at("name").get<std::string>(),
                           fluidJson.at("viscosity").get<REAL>(),
                           fluidJson.at("density").get<REAL>()});
    };
    const json &fluidData = input["FluidData"];
    m_Fluid.Resize(0);
    if (fluidData.is_array())
    {
        for (const auto &fluidJson : fluidData)
            readFluid(fluidJson);
    }
    else if (fluidData.is_object())
    {
        readFluid(fluidData);
    }
    else
    {
        DebugStop();
    }

    if (input.find("PostProcData") == input.end())
        DebugStop();
    json postproc = input["PostProcData"];
    if (postproc.find("wellbore_vtk") == postproc.end())
        DebugStop();
    m_PostProc.wellbore_vtk = postproc["wellbore_vtk"];
    if (postproc.find("reservoir_vtk") == postproc.end())
        DebugStop();
    m_PostProc.reservoir_vtk = postproc["reservoir_vtk"];
    if (postproc.find("training_data") == postproc.end())
        DebugStop();
    m_PostProc.training_data = postproc["training_data"];
    if (postproc.find("vtk_resolution") == postproc.end())
        DebugStop();
    m_PostProc.vtk_resolution = postproc["vtk_resolution"];
    if (postproc.find("training_resolution") == postproc.end())
        DebugStop();
    m_PostProc.training_resolution = postproc["training_resolution"];
    if (postproc.find("nthreads") != postproc.end())
        m_PostProc.nthreads = postproc["nthreads"];
    else
        m_PostProc.nthreads = 0;
    if (postproc.find("verbosityLevel") != postproc.end())
        m_PostProc.verbosityLevel = postproc["verbosityLevel"];
    else
        m_PostProc.verbosityLevel = 0;

    // Default numerics data. Overwrite if present in json file
    m_Numerics.nthreads = 0;
    m_Numerics.NewtonFlow = {1.e-6, 1.e-6, 10, true};
    m_Numerics.NewtonTransport = {1.e-11, 1.e-6, 30, false};
    m_Numerics.reservoirPorder = 1;
    m_Numerics.wellPorder = 2;
    m_Numerics.refPressure = 1e5;
    m_Numerics.pressureScale = 1e5;

    if (input.find("NumericsData") == input.end())
        DebugStop();
    json numerics = input["NumericsData"];
    if (numerics.find("nthreads") != numerics.end())
        m_Numerics.nthreads = numerics["nthreads"];
    auto readNewton = [](const json &settings, NewtonData &data) {
        data.maxIterations = settings.value("maxIterations", data.maxIterations);
        data.res_tol = settings.value("res_tol", data.res_tol);
        data.corr_tol = settings.value("corr_tol", data.corr_tol);
        data.printNewtonProgress = settings.value("printNewtonProgress", data.printNewtonProgress);
    };
    readNewton(numerics.contains("NewtonFlow") ? numerics.at("NewtonFlow") : numerics,
               m_Numerics.NewtonFlow);
    if (numerics.contains("NewtonTransport"))
        readNewton(numerics.at("NewtonTransport"), m_Numerics.NewtonTransport);
    if (numerics.find("reservoirPorder") != numerics.end())
        m_Numerics.reservoirPorder = numerics["reservoirPorder"];
    if (numerics.find("wellPorder") != numerics.end())
        m_Numerics.wellPorder = numerics["wellPorder"];
    if (numerics.find("refPressure") != numerics.end())
        m_Numerics.refPressure = numerics["refPressure"];
    if (numerics.find("pressureScale") != numerics.end())
        m_Numerics.pressureScale = numerics["pressureScale"];

    m_Transport.reset();
    if (input.contains("TransportData")) {
        const auto &settings = input.at("TransportData");
        auto &data = m_Transport.emplace();
        data.initialWaterSaturation = settings.value("initialWaterSaturation", data.initialWaterSaturation);
        data.boundaryWaterSaturation = settings.value("boundaryWaterSaturation", data.boundaryWaterSaturation);
        data.waterFluid = settings.value("waterFluid", data.waterFluid);
        data.oilFluid = settings.value("oilFluid", data.oilFluid);
        data.timeStep = settings.value("timeStep", data.timeStep);
        data.numSteps = settings.value("numSteps", data.numSteps);
        const auto law = settings.value("relativePermeability", std::string("linear"));
        if (law == "quadratic")
            data.relativePermeability = RelativePermeability::Quadratic;
        else if (law != "linear") {
            std::cerr << "Unknown transport relativePermeability: " << law << std::endl;
            DebugStop();
        }
    }

    // After reading the json file, set the auxiliary material IDs
    UpdateAuxiliaryMaterialIds();

    // After reading the json file, we add some additional BCs in the Reservoir Data.
    // This BCs account for the no flux condition on the toe and heel surfaces of each wellbore.
    UpdateReservoirBCsForWellbores();

    // Apply pressure scaling factors
    ApplyPressureScaling();
}

void ProblemData::UpdateAuxiliaryMaterialIds()
{
    int maxID = m_Reservoir.matid;
    for (auto &bcs : m_Reservoir.BCs) {
        if (bcs.second.matid > maxID) maxID = bcs.second.matid;
    }

    for (auto &wellbore : m_Wellbore) {
        if (wellbore.matid > maxID) maxID = wellbore.matid;
        if (wellbore.matidSurf > maxID) maxID = wellbore.matidSurf;
        if (wellbore.matidToeSurf > maxID) maxID = wellbore.matidToeSurf;
        if (wellbore.matidHeelSurf > maxID) maxID = wellbore.matidHeelSurf;
        for (auto &bc : wellbore.BCs) {
            if (bc.second.matid > maxID) maxID = bc.second.matid;
        }
    }

    // Set the auxiliary material IDs to be greater than the maximum existing ID
    EPressure2DSkin = maxID + 1;
    EPressureInterface = maxID + 2;
    EHDivBoundInterface = maxID + 3;
    ETransportInterface = maxID + 4;
}

void ProblemData::UpdateReservoirBCsForWellbores()
{
    for (const auto &wellbore : m_Wellbore) {
        // Add no flux BC for the toe surface
        std::string toeBCName = "no_flux_toe_" + wellbore.name;
        BoundaryData toeBC;
        toeBC.matid = wellbore.matidToeSurf;
        toeBC.type = 1; // Neumann BC
        toeBC.value = 0.0; // No flux
        m_Reservoir.BCs[toeBCName] = toeBC;

        // Add no flux BC for the heel surface
        std::string heelBCName = "no_flux_heel_" + wellbore.name;
        BoundaryData heelBC;
        heelBC.matid = wellbore.matidHeelSurf;
        heelBC.type = 1; // Neumann BC
        heelBC.value = 0.0; // No flux
        m_Reservoir.BCs[heelBCName] = heelBC;
    }
}

void ProblemData::ApplyPressureScaling()
{
    // Update viscosity and density of each fluid based on the pressure scaling factor
    for (auto &fluid : m_Fluid) {
        fluid.viscosity = fluid.viscosity * m_Numerics.pressureScale;
        fluid.density = fluid.density * m_Numerics.pressureScale;
    }

    // Update the pressure boundary conditions
    for (auto &bc : m_Reservoir.BCs) {
        if (bc.second.type == 0) { // Direct BC
            bc.second.value = (bc.second.value - m_Numerics.refPressure) * m_Numerics.pressureScale;
        }
    }
    for (auto &wellbore : m_Wellbore) {
        for (auto &bc : wellbore.BCs) {
            if (bc.second.type == 0) { // Direct BC
                bc.second.value = (bc.second.value - m_Numerics.refPressure) * m_Numerics.pressureScale;
            }
        }
    }
}

void ProblemData::Print(std::ostream &out)
{
    out << "\nSimulation inputs: \n\n";
    out << "Mesh Data:\n";
    out << "File: " << m_Mesh.file << std::endl;
    out << "Number of uniform refinements: " << m_Mesh.NumUniformRef << std::endl;
    out << "Number of directional refinements: " << m_Mesh.NumDirRef << std::endl;
    out << "Cylindrical map: " << (m_Mesh.ToCylindrical ? "Yes" : "No") << std::endl
        << std::endl;

    out << "Wellbore Data:\n";
    int wellboreIndex = 0;
    for (const auto &wellbore : m_Wellbore)
    {
        out << "Wellbore " << wellboreIndex++ << ":\n";
        out << "  Name: " << wellbore.name << std::endl;
        out << "  Material ID: " << wellbore.matid << std::endl;
        out << "  Radius: " << wellbore.radius << std::endl;
        out << "  Length: " << wellbore.length << std::endl;
        out << "  Eccentricity: " << wellbore.eccentricity << std::endl;
        out << "  Boundary conditions:\n";
        for (const auto &bc : wellbore.BCs)
        {
            out << "    Name: " << bc.first << std::endl;
            out << "    Material ID: " << bc.second.matid << std::endl;
            out << "    Type: " << bc.second.type << std::endl;
            out << "    Value: " << bc.second.value << std::endl;
        }
    }
    out << std::endl;

    out << "Reservoir Data:\n";
    out << "Name: " << m_Reservoir.name << std::endl;
    out << "Material ID: " << m_Reservoir.matid << std::endl;
    out << "Permeability: " << m_Reservoir.perm << std::endl;
    out << "Porosity: " << m_Reservoir.porosity << std::endl;
    out << "Height: " << m_Reservoir.height << std::endl;
    out << "Width: " << m_Reservoir.width << std::endl;
    out << "Length: " << m_Reservoir.length << std::endl;
    out << "Boundary conditions:\n";
    for (const auto &bc : m_Reservoir.BCs)
    {
        out << "  Name: " << bc.first << std::endl;
        out << "  Material ID: " << bc.second.matid << std::endl;
        out << "  Type: " << bc.second.type << std::endl;
        out << "  Value: " << bc.second.value << std::endl;
    }
    out << std::endl;

    out << "Fluid Data:\n";
    int fluidIndex = 0;
    for (const auto &fluid : m_Fluid)
    {
        out << "Fluid " << fluidIndex++ << ":\n";
        out << "  Name: " << fluid.name << std::endl;
        out << "  Viscosity: " << fluid.viscosity << std::endl;
        out << "  Density: " << fluid.density << std::endl;
    }
    out << std::endl
        << std::endl;

    out << "Post Processing Data:\n";
    out << "Wellbore VTK: " << m_PostProc.wellbore_vtk << std::endl;
    out << "Reservoir VTK: " << m_PostProc.reservoir_vtk << std::endl;
    out << "Training Data: " << m_PostProc.training_data << std::endl;
    out << "VTK Resolution: " << m_PostProc.vtk_resolution << std::endl;
    out << "Training Data Resolution: " << m_PostProc.training_resolution << std::endl;
    out << "Number of threads: " << m_PostProc.nthreads << std::endl;
    out << "Verbosity Level: " << m_PostProc.verbosityLevel << std::endl
        << std::endl;

    out << "Numerics Data:\n";
    out << "Number of threads: " << m_Numerics.nthreads << std::endl;
    out << "Maximum Iterations: " << m_Numerics.NewtonFlow.maxIterations << std::endl;
    out << "Residual Tolerance: " << m_Numerics.NewtonFlow.res_tol << std::endl;
    out << "Correction Tolerance: " << m_Numerics.NewtonFlow.corr_tol << std::endl;
}

REAL ProblemData::TransportPhaseViscosity(const std::string &name) const {
    const FluidData *phase = nullptr;
    for (const auto &fluid : m_Fluid) {
        if (fluid.name != name) continue;
        if (phase) {
            std::cerr << "Duplicate transport fluid: " << name << '\n';
            DebugStop();
        }
        phase = &fluid;
    }
    if (!phase || !(phase->physicalViscosity > 0.) || !std::isfinite(phase->physicalViscosity)) {
        std::cerr << "Missing or invalid transport fluid: " << name << '\n';
        DebugStop();
    }
    return phase->physicalViscosity;
}

void ProblemData::ValidateTransport() const {
    if (!m_Transport) {
        std::cerr << "Missing TransportData\n";
        DebugStop();
    }
    const auto &settings = *m_Transport;
    const auto &newton = m_Numerics.NewtonTransport;
    if (!(settings.initialWaterSaturation >= 0. && settings.initialWaterSaturation <= 1.) ||
        !(settings.timeStep > 0.) || !std::isfinite(settings.timeStep) || settings.numSteps < 1 ||
        !(m_Reservoir.porosity > 0. && m_Reservoir.porosity <= 1.) ||
        ETransportInterface < 0 || m_Reservoir.matid == ETransportInterface ||
        settings.waterFluid == settings.oilFluid ||
        newton.maxIterations < 1 || !(newton.res_tol > 0.) || !std::isfinite(newton.res_tol)) {
        std::cerr << "Invalid transport parameters\n";
        DebugStop();
    }
    TransportPhaseViscosity(settings.waterFluid);
    TransportPhaseViscosity(settings.oilFluid);
}

std::map<int, REAL> ProblemData::TransportBoundarySaturations() const {
    const auto &overrides = m_Transport->boundaryWaterSaturation;
    std::map<int, REAL> saturations;
    std::map<std::string, int> boundaryIds;
    auto addBoundary = [&](const std::string &name, int matid) {
        const auto [entry, inserted] = boundaryIds.emplace(name, matid);
        if ((!inserted && entry->second != matid && overrides.count(name)) ||
            matid == ETransportInterface || matid == m_Reservoir.matid) {
            std::cerr << "Invalid or ambiguous transport boundary: " << name << '\n';
            DebugStop();
        }
        saturations.emplace(matid, 0.);
    };
    for (const auto &[name, bc] : m_Reservoir.BCs)
        addBoundary(name, bc.matid);
    for (const auto &well : m_Wellbore)
        addBoundary(well.name, well.matidSurf);

    std::map<int, REAL> specified;
    for (const auto &[name, saturation] : overrides) {
        const auto boundary = boundaryIds.find(name);
        if (boundary == boundaryIds.end() || !(saturation >= 0. && saturation <= 1.)) {
            std::cerr << "Unknown transport boundary or invalid saturation: " << name << '\n';
            DebugStop();
        }
        const auto [entry, inserted] = specified.emplace(boundary->second, saturation);
        if (!inserted && entry->second != saturation) {
            std::cerr << "Conflicting transport saturation for boundary: " << name << '\n';
            DebugStop();
        }
        saturations.at(boundary->second) = saturation;
    }
    return saturations;
}

std::pair<REAL, REAL> ProblemData::FractionalFlow(REAL sw) const {
    // Adapted from TSFProblemData::TPetroPhysics (subFlow, MIT; see license).
    const bool quadratic = m_Transport->relativePermeability == RelativePermeability::Quadratic;
    const REAL waterViscosity = TransportPhaseViscosity(m_Transport->waterFluid);
    const REAL oilViscosity = TransportPhaseViscosity(m_Transport->oilFluid);
    const REAL so = 1. - sw;
    const REAL lw = (quadratic ? sw * sw : sw) / waterViscosity;
    const REAL lo = (quadratic ? so * so : so) / oilViscosity;
    const REAL dlw = (quadratic ? 2. * sw : 1.) / waterViscosity;
    const REAL dlo = (quadratic ? -2. * so : -1.) / oilViscosity;
    const REAL total = lw + lo;
    return {lw / total, (dlw * total - lw * (dlw + dlo)) / (total * total)};
}
