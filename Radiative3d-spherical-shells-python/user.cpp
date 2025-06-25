// user.cpp
//
// All compiled code which is intended to be written by the end user
// is, as best can be achieved, stored in this file.  This means a
// user customizing Radiative3D to their own purposes can, hopefully,
// minimize the breadth of changes that they need to make.
//
// An example of "user code" to be stored here is the
// Grid::ConstructGridManual() function, in which a grid can be
// hard-coded (in case a file importer for the desired grid format has
// not yet been written, or else for testing purposes).  This function
// is defined here.
//
// Subordinate user_*_inc.cpp files:
//
// To facilitate organization of multiple custom models, the Makefile
// is coded to trigger on modification of ANY file of pattern
// user_*_inc.cpp.  OK to have one of these per model, and #include
// the subordinate file at a relevant location in this file. Remember
// though to still prototype the function here if needed.
//
#include <iostream>
//#include <filesystem>
#include <cstdlib>      /* exit()     */
#include <algorithm>    /* std::max() */
#include "grid.hpp"
#include <fstream>   // for std::ifstream
#include <map>       // for std::map
#include <vector>    // for std::vector
#include <string>    // for std::string (if you're using it)

// PROTOTYPES:
void LopNorCylinder(Grid &, const std::vector<Real> &);
void LopNorCylinderMoho(Grid &, const std::vector<Real> &);
void LopNorCylinderMoho2(Grid &, const std::vector<Real> &);
void HalfspaceCylinder(Grid &, const std::vector<Real> &);
void ScatParamsStudy(Grid &, const std::vector<Real> &);
void CrustPinchWCG(Grid &, const std::vector<Real> &);
void CrustUpthrustWCG(Grid &, const std::vector<Real> &);
void SphereEarth(Grid &, const std::vector<Real> &);
void ToySphere(Grid &, const std::vector<Real> &);
#include "user_LopNorCyl_inc.cpp"
#include "user_LopNorCylMoho_inc.cpp"
#include "user_LopNorCylMoho2_inc.cpp"
#include "user_Halfspace_inc.cpp"
#include "user_NSCP_inc.cpp"
#include "user_Upthrust_inc.cpp"
#include "user_SphereEarth_inc.cpp"
#include "user_ToySphere_inc.cpp"

//////
// METHOD:  Grid :: ConstructGridManual()
//
//   This user-writable member function is called when a grid source
//   of Grid_COMPILED is requested. Is intended to give the user a way
//   to code a custom grid without too much hassle.  Right now, this
//   is the dominant way to grid models. In the future, a way to read
//   model definition files will be needed for long-term usability.
//
//   Here we have coded this method as a dispatcher function, which
//   just calls out to one of the purpose-motivated custom
//   grid-building functions defined below.
//
void Grid::ConstructGridManual(int Selection, const std::vector<Real> &args) {

  Text LineHead = "|  ConstructGridManual: ";

  if (Selection==0) {
    Selection=5;
    std::cout << LineHead << "Changing selection 0 (no selection) to selection "
              << Selection << ".\n";
  }

  switch (Selection) {
  case 1:
  case 2:                       // Officially ID 1, but set up some aliases in
  case 3:                       // case I want to differentiate in the future
  case 4:
    if (args.size() >= 20) {
      std::cout << LineHead << "Selected Lop Nor Moho Model (Layered).\n";
      LopNorCylinderMoho(*this, args);
    } else {
      std::cout << LineHead << "Selected Lop Nor Baseline Model (Layered).\n";
      LopNorCylinder(*this, args);
    }
    break;

  case 16:
    std::cout << LineHead << "Selected Spherical Earth Model (Spherical).\n";
    SphereEarth(*this, args);
    break;

  case 21:
    std::cout << LineHead << "Selected Lop Nor Moho Model Alt2 (Layered).\n";
    LopNorCylinderMoho2(*this, args);
    break;

  case 5:
  case 6:
  case 7:
    std::cout << LineHead
              << "Selected North Sea Crust Pinch Model (Tetra WCG).\n";
    CrustPinchWCG(*this, args);
    break;

  case 8:
    std::cout << LineHead << "Selected Crust Upthrust Model (Tetra WCG).\n";
    CrustUpthrustWCG(*this, args);
    break;

  case 30:
    std::cout << LineHead << "Selected Spherical Toy Model (Spherical).\n";
    ToySphere(*this, args);
    break;

  case 40:
    std::cout << LineHead << "Selected Halfspace Model (Layered).\n";
    HalfspaceCylinder(*this, args);
    break;

  case 128:
    std::cout << LineHead << "Selected Scatter Params Study (Layered).\n";
    ScatParamsStudy(*this, args);
    break;

  case 666:
    std::cout << "Not using a User hard-Coded grid !" << std::endl;
    // Intended to cause trouble; No grid initialized!
    // Remove in production code.
    break;

  default:
    std::cout << LineHead << "Model selection code not known.\n";
    TextStream errtxt;
    errtxt << "Unknown compiled grid selection ID " << Selection << ".";
    throw(Runtime(errtxt.str()));
    break;
  }

}

// Some grid-making notes:

        // Grids in Radiative3D are, in essence, just a NX-by-NY-by-NZ
        // array of "nodes".  A node (GridNode object) is a point in
        // space at which material properties are "known".  A single
        // node can have attributes (material properties) specified
        // either once or twice.  If specified twice, the node sits on
        // a discontinuity, or "sudden" jump in material properties.

        // Grids are formatted or arranged in different ways, based on
        // the type of "model" one is trying to build.  E.g., a
        // tetrahedral model implies certain expectations about how
        // the grid nodes will be arranged, and a cylinder model has
        // different expectations about how the grid nodes will be
        // arranged.

        // Grid-making for CYLINDER models:

        // A Cylinder grid always has nx=3, ny=1.  The three nodes
        // spanning the x-dimension define a triangle, which establishes
        // the plane boundary between cylinder cells of successive depth.
        // This allows the planar interfaces between cells to have a
        // slope.  It is important to specify the corner locations of the
        // triangle in a counter-clockwise fashion, as seen looking
        // downwards from the +z axis.  The presumption for indexing is
        // that as the z-index increases, the spatial z-coordinates
        // becomes more negative.  The physical attributes of the model
        // for a particular z-level are specified on the ix=0 node.  Any
        // attributes specified on the ix=1 or ix=2 nodes are ignored.
        // MediumCells of one of the Cylinder types span the volume
        // between layers of grid-nodes.  If linear gradients are
        // supported for a particular model attribute (e.g. v_P or v_S),
        // then the top and bottom of the gradient are taken from the
        // node-above, and the node-below, respectively.  If gradients are
        // NOT supported, then the attribute for the cell is taken from
        // the node ABOVE.  (Note: at present gradients are NOT
        // supported for cylinder cells.)
        //
        // It is possible to specify the attributes for each node as
        // many as two times.  If specified only once, the attributes
        // apply to the MediumCells both above and below the node. If
        // specified twice, the first specification applys to the cell
        // above, and the second specification to the cell below.
        // This also signals that full R/T treatement should be
        // applied when a phonon interacts with this interface.

//
//




//////////
// GridFromFile()
//
// Allows to read a grid from a 1D spherical shells grid file (e.g., PREM, AK135, or other planetary models).
// The user can refer to the python/jupyter notebook to create the simplified grids from existing published models. 
// 
//

void Grid::GridFromFile(const Text & file_path) { 
    using namespace std;

    std::cout<< "Entering Function GridFromFile ! " << std::endl;
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
    //Real q = 2000;        //

    //HetSpec HSCr = HSneak(sc_nu, sc_eps, sc_a, sc_k); // Crust
    HetSpec HSMa = HSneak(sc_nu, sc_eps, sc_a, sc_k); // Mantle
    //HetSpec HSCo = HSneak(sc_nu, sc_eps, 2*sc_a, sc_k); // Core Outer
    //HetSpec HSCi = HSneak(sc_nu, sc_eps, sc_a, sc_k); // Core Inner


    bool uHeader=true;
    //vector<Text> uParamNames = {"Depth", "Vp","Vs","Density","Qkappa","Qmu","nNodes"}; // U stands for User provided. Make sure parameters are in this order in the file !!!
    vector<Text> uParamNames = {"Depth", "Vp","Vs","Density","Qkappa","Qmu","nNodes","epsilon","kappa","a","nu"};

    // Map to store vectors corresponding to each parameter
    
   // Text homeDir =  std::filesystem::current_path();
    Text full_path  = file_path;
    std::cout << "the file path: " << full_path << std::endl;\
    ifstream file(full_path);

    // Check if filed opened
    if (!file.is_open()) {
         cerr << "Error: Could not open the file!" <<  endl;
        return;
    }
    
    map<Text, vector<double> > DataMap;
    for (const auto& param : uParamNames) { 
        DataMap[param] = {}; // Initialize an empty vector for each parameter
    }

    // Parse the CSV file into DataMap assuming columns match the order of `uParam_names`
    Text line;

    // Skip the header
    if (uHeader == true) {
        // Read and parse the header
        if (!getline(file, line)) {
            cerr << "Error: File is empty or header is missing!" << endl;
            return;
        }

        stringstream header_stream(line);
        Text param;
        uParamNames.clear();
        while (getline(header_stream, param, ',')) {
            uParamNames.push_back(param);           // Add each header column name to uParamNames
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
   // if (line.empty()) continue;

        stringstream ss(line);
        Text cell;
        int column_index = 0;

        while (getline(ss, cell, ',')) { // Assuming columns are comma-separated
		if (!cell.empty() && column_index < uParamNames.size()) {
			 const Text& param_name = uParamNames[column_index];
			 try {
			        double value = stod(cell);
			        DataMap[param_name].push_back(value);
			    } catch (const std::exception& e) {
			        cerr << "Error converting '" << cell << "' to double at column " << column_index << " (" << param_name << "): " << e.what() << endl;
				    }
			}
            
column_index++;
        }

                                                // Ensure no extra columns exist beyond what the user provided
        if (column_index > uParamNames.size()) {
            cerr << "Warning: Row has more columns than expected based on parameter list!" << endl;
        }
    }

    file.close();

    //std::cout << "DataMap size: " << DataMap.size() << std::endl;
    // Get the data

    vector<double> depth = DataMap["Depth"];
    vector<double> Vp = DataMap["Vp"];
    vector<double> Vs = DataMap["Vs"];
    vector<double> Density = DataMap["Density"];
    vector<double> Qm = DataMap["Qmu"];
    vector<double> Qk = DataMap["Qkappa"];
    int nNode = DataMap["nNodes"][0];
    vector<double> epsilon = DataMap["epsilon"];
    vector<double> kappa   = DataMap["kappa"];
    vector<double> a       = DataMap["a"];
    vector<double> nu      = DataMap["nu"];

    
    // assign a very small value to Vs if it is 0 :
    for (int i = 0; i < Vs.size(); i++) {
        if (Vs[i] == 0) {
            Vs[i] = 0.00001;
        }
    }
    
    //print the data
    /*
    for (int i = 0; i < depth.size(); i++) {
        cout << depth[i] << " " << Vp[i] << " " << Vs[i] << " " << Density[i] << " " << Qk[i] << " " << Qm[i] << endl;
    }
    */
    this->SetSize(1,1,nNode);         // Sets index bounds. Using this-> because already GridFromFile is a memeber function of the Grid Object.
    this->SetIndexBase(0);           // When addressing nodes, use base 0
    this->SetMapping(Grid::GC_RAE, Grid::GC_SPHERICAL);

    
    Index NodeIdx =0;                          // Loop over the parameters to parse the model
    for (Index i = 0; i < depth.size() ; i+=2) {// Read the first value and parse the first node. 
        //if (i==0){                             // Update the CurrentNode counter as well as index counter.
                                               // Need to evaluate separately at 0 and end because there is not detectable discontinuity.

      //      std::cout << "First iteration TOP" << std::endl;                                
//            std::cout << "NodeIdx: " << NodeIdx << std::endl;
  //          std::cout << "Depth: " << depth[i] << std::endl;
    //        std::cout << "Vp : " << Vp[i] << std::endl;


            // Define the heterogeneity spectrum of the layer
            HetSpec HS_top = HSneak(nu[i], epsilon[i], a[i], kappa[i]);
            HetSpec HS_bot = HSneak(nu[i+1], epsilon[i+1], a[i+1], kappa[i+1]);
            this->WNode(0,0,NodeIdx).SetLocation ( 0, 0, depth[i]);                                                             //Set the node for the first layer
            this->WNode(0,0,NodeIdx).SetAttributes( VpVs( Vp[i], Vs[i]), Density[i], QmQk( Qm[i], Qk[i]), HS_top );              // Set the parameters for the first layer
            
            // Print statement to check the values
            std::cout << "Scattering params at depth " << depth[i] << ": epsilon=" << epsilon[i] 
            << " kappa=" << kappa[i] << " a=" << a[i] << " nu=" << nu[i] << std::endl;
  
            NodeIdx++; //update the node index
            this->WNode(0,0,NodeIdx).SetLocation ( 0, 0, -depth[i+1] );                                                       //Set the node for the bottom of the first layer
            this->WNode(0,0,NodeIdx).SetAttributes( VpVs( Vp[i+1], Vs[i+1]), Density[i+1], QmQk( Qm[i+1], Qk[i+1]), HS_bot ); // Set the parameters for the bottom of the first layer

            // Print statement to check the values
        //    std::cout << "First iteration BOTTOM" << std::endl;
          //  std::cout << "NodeIdx: " << NodeIdx << std::endl;
           // std::cout << "Depth: " << depth[i+1] << std::endl;
            //std::cout << "Vp : " << Vp[i+1] << std::endl;

            //std::cout << "Scattering params at depth " << depth[i+1] << ": epsilon=" << epsilon[i+1] 
         // << " kappa=" << kappa[i+1] << " a=" << a[i+1] << " nu=" << nu[i+1] << std::endl;
    }
    std::cout<< "Exiting Function GridFromFile ! " << std::endl;
    //print the node values:
     
    return ; // Exit the function after reading the file and setting the grid nodes.
}
    

//////
// ScatParamsStudy()
//
//   This builds a grid intended to study the effects of varying
//   scattering parameters on mean-free-paths (MFPs) and/or other
//   scattering metrics.
//
//   This grid isn't intended for actual phonon spraying, but rather
//   just to get info from the scatterer-object dump that can be
//   requested at runtime. (Lists all scatterers, their paramters, and
//   MFP's)
//
//   Builds grid as a series of layers in which all parameters but one
//   are held fixed, while one parameter varies.
//
//   We assume a CYLINDER model and the grid-formatting conventions
//   thereof, as this is the easiest way to build a simple layered
//   Earth model.
//
//   The 'args' list is ignored, but could be used to porvide
//   midranges if we wanted to code it that way.
//
void ScatParamsStudy(Grid & gr, const std::vector<Real> & args) {

  using Elastic::HSneak;
  using Elastic::VpVs;
  using Elastic::QmQk;

  Real LayerSpacing = 10;

  // "Mid-Range" Parameters:
  Real nu  = 0.8;
  Real eps = 0.05;
  Real a   = 1.0;
  Real k   = 0.3;
  Real vs  = 4.0;             // El = omega/vs
  Real vp  = vs*1.7321;       // gam0 = vp/vs

  Elastic::HSneak HS(nu,eps,a,k); //
  Elastic::VpVs vpvs(vp,vs);  // 
  Elastic::QmQk Q(1.0/0.0);   // (Infinite. We don't vary this.)
  Real rho = 1.0;             // (We don't vary this either.)

  // Grid Dimension:
  //
  //   If we vary six parameters, with seven data-points each, we need
  //   42 model layers, and 43 grid-sheets to bound them.
  //
  Count ndp = 7;        // Num data points (per parameter)
  Count npar = 6;       // Num parameters that we vary (nu,eps,a,k,el,gam0)
  Count nsh = ndp*npar+1; // Num grid sheets needed
  gr.SetSize(3,1,nsh);  // (Cylinder grid requires NX=3, NY=1)
  gr.SetIndexBase(0);   // Count from 0

  // Set all node LOCATIONS first. (Then we will set attributes.)
  //
  for (Index i = 0; i<nsh; i++) {
    Real depth = -1.0*LayerSpacing*i;
    gr.WNode(0,0,i).SetLocation(0.0, 0.0, depth);
    gr.WNode(1,0,i).SetLocation(1.0, 0.0, depth);
    gr.WNode(2,0,i).SetLocation(0.0, 1.0, depth);
  }

  // Now set grid attributes:
  //
  Index iz = 0;
  //
  //  Varying Nu:
  gr.WNode(0,0,iz++).SetAttributes(vpvs, rho, Q, HSneak(0.3, eps, a, k));
  gr.WNode(0,0,iz++).SetAttributes(vpvs, rho, Q, HSneak(0.5, eps, a, k));
  gr.WNode(0,0,iz++).SetAttributes(vpvs, rho, Q, HSneak(0.7, eps, a, k));
  gr.WNode(0,0,iz++).SetAttributes(vpvs, rho, Q, HSneak(0.9, eps, a, k));
  gr.WNode(0,0,iz++).SetAttributes(vpvs, rho, Q, HSneak(1.1, eps, a, k));
  gr.WNode(0,0,iz++).SetAttributes(vpvs, rho, Q, HSneak(1.3, eps, a, k));
  gr.WNode(0,0,iz++).SetAttributes(vpvs, rho, Q, HSneak(1.5, eps, a, k));
  //
  //  Varying Eps:
  gr.WNode(0,0,iz++).SetAttributes(vpvs, rho, Q, HSneak(nu, 0.01, a, k));
  gr.WNode(0,0,iz++).SetAttributes(vpvs, rho, Q, HSneak(nu, 0.02, a, k));
  gr.WNode(0,0,iz++).SetAttributes(vpvs, rho, Q, HSneak(nu, 0.03, a, k));
  gr.WNode(0,0,iz++).SetAttributes(vpvs, rho, Q, HSneak(nu, 0.04, a, k));
  gr.WNode(0,0,iz++).SetAttributes(vpvs, rho, Q, HSneak(nu, 0.05, a, k));
  gr.WNode(0,0,iz++).SetAttributes(vpvs, rho, Q, HSneak(nu, 0.06, a, k));
  gr.WNode(0,0,iz++).SetAttributes(vpvs, rho, Q, HSneak(nu, 0.07, a, k));
  //
  //  Varying A:
  gr.WNode(0,0,iz++).SetAttributes(vpvs, rho, Q, HSneak(nu, eps, 0.3, k));
  gr.WNode(0,0,iz++).SetAttributes(vpvs, rho, Q, HSneak(nu, eps, 0.6, k));
  gr.WNode(0,0,iz++).SetAttributes(vpvs, rho, Q, HSneak(nu, eps, 0.8, k));
  gr.WNode(0,0,iz++).SetAttributes(vpvs, rho, Q, HSneak(nu, eps, 1.0, k));
  gr.WNode(0,0,iz++).SetAttributes(vpvs, rho, Q, HSneak(nu, eps, 1.2, k));
  gr.WNode(0,0,iz++).SetAttributes(vpvs, rho, Q, HSneak(nu, eps, 1.5, k));
  gr.WNode(0,0,iz++).SetAttributes(vpvs, rho, Q, HSneak(nu, eps, 1.8, k));
  //
  //  Varying kappa:
  gr.WNode(0,0,iz++).SetAttributes(vpvs, rho, Q, HSneak(nu, eps, a, 0.1));
  gr.WNode(0,0,iz++).SetAttributes(vpvs, rho, Q, HSneak(nu, eps, a, 0.2));
  gr.WNode(0,0,iz++).SetAttributes(vpvs, rho, Q, HSneak(nu, eps, a, 0.3));
  gr.WNode(0,0,iz++).SetAttributes(vpvs, rho, Q, HSneak(nu, eps, a, 0.4));
  gr.WNode(0,0,iz++).SetAttributes(vpvs, rho, Q, HSneak(nu, eps, a, 0.5));
  gr.WNode(0,0,iz++).SetAttributes(vpvs, rho, Q, HSneak(nu, eps, a, 0.6));
  gr.WNode(0,0,iz++).SetAttributes(vpvs, rho, Q, HSneak(nu, eps, a, 0.8));
  //
  //  Varying El, keeping gam0 fixed:   (el ~ 1/vs)
  gr.WNode(0,0,iz++).SetAttributes(VpVs(1.5*vp, 1.5*vs), rho, Q, HS);
  gr.WNode(0,0,iz++).SetAttributes(VpVs(1.3*vp, 1.3*vs), rho, Q, HS);
  gr.WNode(0,0,iz++).SetAttributes(VpVs(1.1*vp, 1.1*vs), rho, Q, HS);
  gr.WNode(0,0,iz++).SetAttributes(VpVs(1.0*vp, 1.0*vs), rho, Q, HS);
  gr.WNode(0,0,iz++).SetAttributes(VpVs(0.9*vp, 0.9*vs), rho, Q, HS);
  gr.WNode(0,0,iz++).SetAttributes(VpVs(0.7*vp, 0.7*vs), rho, Q, HS);
  gr.WNode(0,0,iz++).SetAttributes(VpVs(0.5*vp, 0.5*vs), rho, Q, HS);
  //  
  //  Varying gam0:
  gr.WNode(0,0,iz++).SetAttributes(VpVs(1.3321*vs, vs), rho, Q, HS);
  gr.WNode(0,0,iz++).SetAttributes(VpVs(1.5321*vs, vs), rho, Q, HS);
  gr.WNode(0,0,iz++).SetAttributes(VpVs(1.6321*vs, vs), rho, Q, HS);
  gr.WNode(0,0,iz++).SetAttributes(VpVs(1.7321*vs, vs), rho, Q, HS);
  gr.WNode(0,0,iz++).SetAttributes(VpVs(1.8321*vs, vs), rho, Q, HS);
  gr.WNode(0,0,iz++).SetAttributes(VpVs(1.9321*vs, vs), rho, Q, HS);
  gr.WNode(0,0,iz++).SetAttributes(VpVs(2.1321*vs, vs), rho, Q, HS);
  //
  //   Bottom sheet: (values don't matter)
  gr.WNode(0,0,iz++).SetAttributes(vpvs, rho, Q, HS);

}
