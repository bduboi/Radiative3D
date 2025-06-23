#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "typedefs.hpp"
#include <grid.hpp>
#include <map>

using namespace std;

void SphereEarthModel(Grid & gr, const string & file_path) {


    using Elastic::Velocity;
    using Elastic::VpVs;
    using Elastic::Q;
    using Elastic::QmQk;
    using Elastic::HetSpec;
    using Elastic::HSneak;

    Real sc_nu  = 0.8;    // Default scat params -
    Real sc_eps = 0.005;  // we'll scan args for user values
    Real sc_a   = 4.00;   //
    Real sc_k = 0.8;      //
    Real q = 2000;        //

    HetSpec HSCr = HSneak(sc_nu, sc_eps, sc_a, sc_k); // Crust
    HetSpec HSMa = HSneak(sc_nu, sc_eps, sc_a, sc_k); // Mantle
    HetSpec HSCo = HSneak(sc_nu, sc_eps, 2*sc_a, sc_k); // Core Outer
    HetSpec HSCi = HSneak(sc_nu, sc_eps, sc_a, sc_k); // Core Inner


    bool uHeader=true;
    vector<string> uParamNames = {"Depth", "Vp","Vs","Density","Qmu","Qkappa","nNodes"}; // U stands for User provided
    // Map to store vectors corresponding to each parameter
    
    string homeDir =  "/Volumes/2024_Dubois/Project_ANU/Project_radiative3D/Radiative3D-spherical-shells/Models";
    //string file_path = "/MoonModels/SimplfiedModels/SimplfiedMoonWeber_Science_2011.csv"; // Path to your CSV file
    string full_path  = homeDir + file_path;
    ifstream file(full_path);

    // Check if filed opened
    if (!file.is_open()) {
         cerr << "Error: Could not open the file!" <<  endl;
        return;
    }
    
    map<string, vector<double> > DataMap;
    for (const auto& param : uParamNames) { 
        DataMap[param] = {}; // Initialize an empty vector for each parameter
    }

    // Parse the CSV file into DataMap assuming columns match the order of `uParam_names`
    string line;

    // Skip the header
    if (uHeader == true) {
        // Read and parse the header
        if (!getline(file, line)) {
            cerr << "Error: File is empty or header is missing!" << endl;
            return;
        }

        stringstream header_stream(line);
        string param;
        uParamNames.clear();
        while (getline(header_stream, param, ',')) {
            uParamNames.push_back(param); // Add each header column name to uParamNames
        }
    } 
    else {
        // Check if file is empty
        if (!getline(file, line)) {
            cerr << "Error: File is empty!" << endl;
            return;
        }
    }

    while (getline(file, line)) {
        stringstream ss(line);
        string cell;
        int column_index = 0;

        while (getline(ss, cell, ',')) { // Assuming columns are comma-separated
            if (!cell.empty() && column_index < uParamNames.size()) {
                const string& param_name = uParamNames[column_index];
                DataMap[param_name].push_back(stod(cell)); // Convert and store the value
            }
            column_index++;
        }

        // Ensure no extra columns exist beyond what the user provided
        if (column_index > uParamNames.size()) {
            cerr << "Warning: Row has more columns than expected based on parameter list!" << endl;
        }
    }

    file.close();


    //Print parsed data
    cout << "Parsed Data:\n";
    for (Index i =0; i< uParamNames.size()-1; ++i){
        string ParamName =  uParamNames[i];
        if (ParamName == "nNodes"){
            cout << "\n" << ParamName << " = " << DataMap[ParamName][0]<< endl;
            int nNode = DataMap[ParamName][0]
        }
        else {
        cout << "\n" << ParamName << " Size = " << DataMap[ParamName].size()<< endl;
        for (const auto& val : DataMap[ParamName]) cout << val << " ";
        }
        cout << "\n";
    }


    vector<double> depth = DataMap["Depth"];
    vector<double> Vp = DataMap["Vp"];
    vector<double> Vs = DataMap["Vs"];
    vector<double> Density = DataMap["Density"];
    vector<double> Qm = DataMap["Qmu"];
    vector<double> Qk = DataMap["Qkappa"];

    //vector<double> depth = DataMap["depth"];

    gr.SetSize(1,1,nNode);           // Sets index bounds
    gr.SetIndexBase(0);           // When addressing nodes, use base 0
    gr.SetMapping(Grid::GC_RAE, Grid::GC_SPHERICAL);

    // Loop over the parameters to parse the model
    Index NodeIdx =0;
    for (Index i = 0; i < depth.size() ; ++i) {
        if (i==0){
        // Read the first value and parse the first node. 
        // Update the CurrentNode counter as well as index counter.
        /*
            cout <<  "------- first iter ------- Node number : " << NodeIdx+1 << endl;
            cout <<  "i =  " << i <<  " depth =  " << depth[i] << endl;
            cout <<  "i =  " << i <<  " Vp =  " << Vp[i] << endl;
            cout <<  "i =  " << i <<  " Vs =  " << Vs[i] << endl;
        */
            gr.WNode(0,0,NodeIdx).SetLocation ( 0, 0, depth);
            gr.WNode(0,0,NodeIdx).SetAttributes( VpVs( Vp[i], Vs[i]), Density[i], QmQk( Qm[i], Qk[i]), HSCr );

            }
        // Detect the discontinuity (1) and (2) prevents from parsing the last element twice.
        if (depth[i + 1] == depth[i] && (i<depth.size() -1 )) { 
        // Read both the values for the Node and update the Node counter.
        // The index counter needs to be adjusted so that it does not overlap.
        
        /* Print statements to check for correct parsing
            cout <<  "------ intermediate iter ------- Node number : " << NodeIdx+1 << endl;
            cout <<  "i =  " << i <<  " depth =  " << depth[i] << endl;
            cout <<  "i =  " << i <<  " Vp =  " << Vp[i] << endl;
            cout <<  "i =  " << i <<  " Vs =  " << Vs[i] << endl;
            cout <<  "i+1 =  " << i+1 <<  " depth =  " << depth[i+1] << endl;
            cout <<  "i+1 =  " << i+1 <<  " Vp =  " << Vp[i+1] << endl;
            cout <<  "i+1 =  " << i+1 <<  " Vs =  " << Vs[i+1] << endl;
        */
            gr.WNode(0,0,NodeIdx).SetLocation ( 0, 0, -depth[i] );   // Moho
            gr.WNode(0,0,NodeIdx).SetAttributes( VpVs( Vp[i], Vs[i]), Density[i], QmQk( Qm[i], Qk[i]), HSCr );
            i++; // Increment the model index to skip the double layer.
            cout <<  "i after increment " << i <<  " depth =  " << depth[i] << endl;
            cout <<  "------- end intermediate iter ------- \n" << endl;

        }
        // Normal statement in the case there is no discontinuity (last element)

        if (depth[i]==*std::max_element(depth.begin(), depth.end())){ // Read the last values and parse last node 
            /*
            cout <<  "------- last iter ------- Node number : " << NodeIdx+1 << endl;
            cout <<  " i =  " << i <<  " depth =  " << depth[i] << endl;       
            cout <<  " i =  " << i <<  " Vs =  " << Vp[i] << endl;       
            cout <<  " i =  " << i <<  " Vs =  " << Vs[i] << endl;      
            */
            gr.WNode(0,0,NodeIdx).SetLocation ( 0, 0, -depth[i] );  // Center 
            gr.WNode(0,0,NodeIdx).SetAttributes( VpVs( Vp[i], Vs[i]), Density[i], QmQk( Qm[i], Qk[i]), HSCr );
        }
        NodeIdx++; // Update the node index to parse next layer

    }
    return 0;
}
    