
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
        FluidData fluid;
        if (fluidJson.find("name") == fluidJson.end())
            DebugStop();
        fluid.name = fluidJson["name"];
        if (fluidJson.find("viscosity") == fluidJson.end())
            DebugStop();
        fluid.viscosity = fluidJson["viscosity"];
        if (fluidJson.find("density") == fluidJson.end())
            DebugStop();
        fluid.density = fluidJson["density"];
        m_Fluid.push_back(fluid);
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
    m_Numerics.maxIterations = 10;
    m_Numerics.res_tol = 1e-6;
    m_Numerics.corr_tol = 1e-6;
    m_Numerics.reservoirPorder = 1;
    m_Numerics.wellPorder = 2;
    m_Numerics.refPressure = 1e5;
    m_Numerics.pressureScale = 1e5;

    if (input.find("NumericsData") == input.end())
        DebugStop();
    json numerics = input["NumericsData"];
    if (numerics.find("nthreads") != numerics.end())
        m_Numerics.nthreads = numerics["nthreads"];
    if (numerics.find("maxIterations") != numerics.end())
        m_Numerics.maxIterations = numerics["maxIterations"];
    if (numerics.find("res_tol") != numerics.end())
        m_Numerics.res_tol = numerics["res_tol"];
    if (numerics.find("corr_tol") != numerics.end())
        m_Numerics.corr_tol = numerics["corr_tol"];
    if (numerics.find("reservoirPorder") != numerics.end())
        m_Numerics.reservoirPorder = numerics["reservoirPorder"];
    if (numerics.find("wellPorder") != numerics.end())
        m_Numerics.wellPorder = numerics["wellPorder"];
    if (numerics.find("refPressure") != numerics.end())
        m_Numerics.refPressure = numerics["refPressure"];
    if (numerics.find("pressureScale") != numerics.end())
        m_Numerics.pressureScale = numerics["pressureScale"];

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
    out << "Maximum Iterations: " << m_Numerics.maxIterations << std::endl;
    out << "Residual Tolerance: " << m_Numerics.res_tol << std::endl;
    out << "Correction Tolerance: " << m_Numerics.corr_tol << std::endl;
}
