
#include <iostream>
#include <string>
#include <fstream>
#include <pzerror.h>

#include "ProblemData.h"

using namespace std;

// constructor
ProblemData::ProblemData() :
    m_MeshName(""),
    m_NumUniformRef(0),
    m_NumDirRef(0),
    m_Dimension(3),
    m_Resolution(0),
    m_ToCylindrical(false),
    m_Verbose(0)
    
    {
      m_Wellbore.BCs.resize(2);
      m_Reservoir.BCs.resize(3);
    }


// deconstructor
ProblemData::~ProblemData() {}

// readjson function. takes a json function as parameter and completes the required simulation data
void ProblemData::ReadJson(std::string file){
    std::string path(std::string(INPUTDIR) + "/" + file);
    std::ifstream filejson(path);
    json input = json::parse(filejson,nullptr,true,true); // to ignore comments in json file
    
    if(input.find("MeshName") == input.end()) DebugStop();
    m_MeshName = input["MeshName"];
    if(input.find("Dimension") == input.end()) DebugStop();
    m_Dimension = input["Dimension"];
    if(input.find("NumUniformRef") == input.end()) DebugStop();
    m_NumUniformRef = input["NumUniformRef"];
    if(input.find("NumDirRef") == input.end()) DebugStop();
    m_NumDirRef = input["NumDirRef"];
    if(input.find("Resolution") == input.end()) DebugStop();
    m_Resolution = input["Resolution"];
    if(input.find("ToCylindrical") == input.end()) DebugStop();
    m_ToCylindrical = input["ToCylindrical"];

    if(input.find("WellboreData") == input.end()) DebugStop();
    json wellbore = input["WellboreData"];
    if(wellbore.find("name") == wellbore.end()) DebugStop();
    m_Wellbore.name = wellbore["name"];
    if(wellbore.find("matid") == wellbore.end()) DebugStop();
    m_Wellbore.matid = wellbore["matid"];
    if(wellbore.find("perm") == wellbore.end()) DebugStop();
    m_Wellbore.perm = wellbore["perm"];
    if(wellbore.find("pOrder") == wellbore.end()) DebugStop();
    m_Wellbore.pOrder = wellbore["pOrder"];
    if(wellbore.find("BCs") == wellbore.end()) DebugStop();
    json bcs = wellbore["BCs"];
    for(int i = 0; i < bcs.size(); i++){
        if(bcs[i].find("name") == bcs[i].end()) DebugStop();
        m_Wellbore.BCs[i].name = bcs[i]["name"];
        if(bcs[i].find("matid") == bcs[i].end()) DebugStop();
        m_Wellbore.BCs[i].matid = bcs[i]["matid"];
        if(bcs[i].find("type") == bcs[i].end()) DebugStop();
        m_Wellbore.BCs[i].type = bcs[i]["type"];
        if(bcs[i].find("value") == bcs[i].end()) DebugStop();
        m_Wellbore.BCs[i].value = bcs[i]["value"];
    }
    if(wellbore.find("radius") == wellbore.end()) DebugStop();
    m_Wellbore.radius = wellbore["radius"];
    if(wellbore.find("length") == wellbore.end()) DebugStop();
    m_Wellbore.length = wellbore["length"];

    if(input.find("ReservoirData") == input.end()) DebugStop();
    json reservoir = input["ReservoirData"];
    if(reservoir.find("name") == reservoir.end()) DebugStop();
    m_Reservoir.name = reservoir["name"];
    if(reservoir.find("matid") == reservoir.end()) DebugStop();
    m_Reservoir.matid = reservoir["matid"];
    if(reservoir.find("perm") == reservoir.end()) DebugStop();
    m_Reservoir.perm = reservoir["perm"];
    if(reservoir.find("pOrder") == reservoir.end()) DebugStop();
    m_Reservoir.pOrder = reservoir["pOrder"];
    if(reservoir.find("BCs") == reservoir.end()) DebugStop();
    json bcsres = reservoir["BCs"];
    for(int i = 0; i < bcsres.size(); i++){
        if(bcsres[i].find("name") == bcsres[i].end()) DebugStop();
        m_Reservoir.BCs[i].name = bcsres[i]["name"];
        if(bcsres[i].find("matid") == bcsres[i].end()) DebugStop();
        m_Reservoir.BCs[i].matid = bcsres[i]["matid"];
        if(bcsres[i].find("type") == bcsres[i].end()) DebugStop();
        m_Reservoir.BCs[i].type = bcsres[i]["type"];
        if(bcsres[i].find("value") == bcsres[i].end()) DebugStop();
        m_Reservoir.BCs[i].value = bcsres[i]["value"];
    }
    if(reservoir.find("height") == reservoir.end()) DebugStop();
    m_Reservoir.height = reservoir["height"];
    if(reservoir.find("width") == reservoir.end()) DebugStop();
    m_Reservoir.width = reservoir["width"];
    if(reservoir.find("length") == reservoir.end()) DebugStop();
    m_Reservoir.length = reservoir["length"];
}

void ProblemData::Print(std::ostream& out){
    out << "\nSimulation inputs: \n\n";
    out << "Mesh Name: " << m_MeshName << std::endl << std::endl;
    out << "Dimension: " << m_Dimension << std::endl << std::endl;
    out << "Number of uniform refinements: " << m_NumUniformRef << std::endl << std::endl;
    out << "Number of directional refinements: " << m_NumDirRef << std::endl << std::endl;
    out << "VTK mesh resolution: " << m_Resolution << std::endl << std::endl;
    out << "Cylindrical map: " << m_ToCylindrical << std::endl << std::endl;
    out << "Wellbore data: \n";
    out << "Name: " << m_Wellbore.name << std::endl;
    out << "Material ID: " << m_Wellbore.matid << std::endl;
    out << "Permeability: " << m_Wellbore.perm << std::endl;
    out << "Polynomial order: " << m_Wellbore.pOrder << std::endl;
    out << "Radius: " << m_Wellbore.radius << std::endl;
    out << "Length: " << m_Wellbore.length << std::endl;
    out << "Boundary conditions: \n";
    for(int i = 0; i < m_Wellbore.BCs.size(); i++){
        out << "Name: " << m_Wellbore.BCs[i].name << std::endl;
        out << "Material ID: " << m_Wellbore.BCs[i].matid << std::endl;
        out << "Type: " << m_Wellbore.BCs[i].type << std::endl;
        out << "Value: " << m_Wellbore.BCs[i].value << std::endl;
    }
    out << "Reservoir data: \n";
    out << "Name: " << m_Reservoir.name << std::endl;
    out << "Material ID: " << m_Reservoir.matid << std::endl;
    out << "Permeability: " << m_Reservoir.perm << std::endl;
    out << "Polynomial order: " << m_Reservoir.pOrder << std::endl;
    out << "Height: " << m_Reservoir.height << std::endl;
    out << "Width: " << m_Reservoir.width << std::endl;
    out << "Length: " << m_Reservoir.length << std::endl;
    out << "Boundary conditions: \n";
    for(int i = 0; i < m_Reservoir.BCs.size(); i++){
        out << "Name: " << m_Reservoir.BCs[i].name << std::endl;
        out << "Material ID: " << m_Reservoir.BCs[i].matid << std::endl;
        out << "Type: " << m_Reservoir.BCs[i].type << std::endl;
        out << "Value: " << m_Reservoir.BCs[i].value << std::endl;
    }
}