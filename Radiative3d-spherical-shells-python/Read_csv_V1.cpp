#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "typedefs.hpp"
#include "grid.hpp"

using namespace std;

int main() {
     string homeDir =  "/Volumes/2024_Dubois/Project_ANU/Project_radiative3D/Radiative3D-spherical-shells/Models";
     string file_path = "/MoonModels/SimplfiedModels/SimplfiedMoonWeber_Science_2011.csv"; // Path to your CSV file
     string full_path  = homeDir + file_path;
     ifstream file(full_path);

    if (!file.is_open()) {
         cerr << "Error: Could not open the file!" <<  endl;
        return 0;
    }

     string line;

     vector< string> headers;
     vector<double> depth, density, vp, vs, qmu, qkappa;

    // Read the first line to get the headers
    if ( getline(file, line)) {
         stringstream ss(line);
         string header;
        while ( getline(ss, header,  ',')) { // Assuming headers are separated by spaces
            headers.push_back(header);
        }
    }

    // Read the rest of the file to parse data
    Index nNode;
    int currentNode=0;
    bool is_first_row = true;
    while ( getline(file, line)) {
        stringstream ss(line);
         vector<double> row;
         string cell;
        while ( getline(ss, cell, ',')) { // Assuming values are separated by spaces
            if (!cell.empty()) {
                row.push_back( stod(cell)); // Convert each cell to a double
            }
        }

        // Store the values in their respective vectors
        if (row.size() >= 1) { // Ensure there are at least 10 columns
            depth.push_back(row[0]);
            density.push_back(row[1]);
            vp.push_back(row[2]);
            vs.push_back(row[3]);        
            qmu.push_back(row[4]);
            qkappa.push_back(row[5]);
            if (is_first_row) {
                nNode = row[6];
                is_first_row=false;
            }

            //// READ THE FIRST VALUE AND PARSE

            
        }
    }

    file.close();

    //Print parsed data
    cout << "Parsed Data:\n";
    cout << "\nDepth: ";
    for (const auto& val : depth) cout << val << " ";
    cout << "\nDensity: ";
    for (const auto& val : density) cout << val << " ";
    cout << "\nnNode : " << nNode <<  endl;
    cout << "\n";


    // Loop over the parameters to parse the model
    Index NodeIdx =0;
    for (Index i = 0; i < depth.size() ; ++i) {
        if (i==0){
        // Read the first value and parse the first node. 
        // Update the CurrentNode counter as well as index counter.
            cout <<  "------- first iter ------- Node number : " << NodeIdx+1 << endl;
            cout <<  "i =  " << i <<  " depth =  " << depth[i] << endl;

            
        }

        // Detect the discontinuity

        if (depth[i + 1] == depth[i]) {
        // Read both the values for the Node and update the Node counter.
        // The index counter needs to be adjusted so that it does not overlap.
            cout <<  "------ intermediate iter ------- Node number : " << NodeIdx+1 << endl;
            cout <<  "i =  " << i <<  " depth =  " << depth[i] << endl;
            cout <<  "i+1 =  " << i+1 <<  " depth =  " << depth[i+1] << endl;
            i++; // Increment the model index to skip the double layer.
            cout <<  "i after increment " << i <<  " depth =  " << depth[i] << endl;
            cout <<  "------- end intermediate iter ------- \n" << endl;

        }
        // Normal statement in the case there is no discontinuity

        if (depth[i]==*std::max_element(depth.begin(), depth.end())){ // Read the last values and parse last node 
            cout <<  "------- last iter ------- Node number : " << NodeIdx+1 << endl;
            cout <<  " i =  " << i <<  " depth =  " << depth[i] << endl;       
        }
        NodeIdx++; // Update the node index to parse next layer

    }
    return 0;
}
    