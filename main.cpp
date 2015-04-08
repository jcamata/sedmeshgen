/* 
 * File:   main.cpp
 * Author: camata
 *
 * Created on 11 de Agosto de 2014, 11:06
 */

#include <cstdlib>

using namespace std;

#include "petsc.h"
#include "mesh.h"
#include "gmsh.h"
#include "fsi.h"
#include "edgecfd_writer.h"


#ifdef _MESQUITE_

#include "Mesquite.hpp"
#include "MsqError.hpp"
#include "MeshImpl.hpp"
#include "MeshInterface.hpp"
#include "ArrayMesh.hpp"
#include "PlanarDomain.hpp"
#include "ShapeImprover.hpp"
#include "LaplaceWrapper.hpp"
#include "LaplacianSmoother.hpp"
#include "InstructionQueue.hpp"
#include "ConditionNumberQualityMetric.hpp"
#include "EdgeLengthQualityMetric.hpp"
#include "QualityAssessor.hpp"
#include "TerminationCriterion.hpp"


void OptGMesh(Gmsh &gmesh, Mesh &mesh)
{
    
    std::vector<int> fixed(gmesh.getNumberofNodes());
    
    std::fill(fixed.begin(),fixed.end(), 0);
    
    for(int iel = 0; iel < gmesh.getNumberofFaces(); iel++)
    {
        std::vector<int> conn = gmesh.getFaceConnectivity(iel);
        for (int ino = 0; ino < conn.size(); ino++)
            fixed[conn[ino]] = 1;
    }
    
    int dim = 3;
    unsigned long int nnodes = mesh.get_n_nodes();
    unsigned long int nelem  = mesh.get_n_elements();
    unsigned long   *conn  = mesh.get_conn_pointer();
    double *coords = mesh.get_coords_pointer();
    int *fptr = &fixed[0];
    
    
    Mesquite::Vector3D normal(0,0,1);
    Mesquite::Vector3D point(7,6,2.38);

    
    
    Mesquite::ArrayMesh mesq_mesh(dim, nnodes, coords, fptr, nelem, Mesquite::TETRAHEDRON, conn);
    //Mesquite::PlanarDomain domain(normal,point);
    Mesquite::MsqError error;
    //Mesquite::MeshDomainAssoc mesh_domain = Mesquite::MeshDomainAssoc(&mesq_mesh, &domain);
    
  
    //Mesquite::ShapeImprover shape_wrapper; 
    Mesquite::LaplaceWrapper optmizer;
  
    /*

    Mesquite::LaplacianSmoother lapl1;
    Mesquite::InstructionQueue q1;
    Mesquite::ConditionNumberQualityMetric shape_metric;
    Mesquite::EdgeLengthQualityMetric lapl_met;
    lapl_met.set_averaging_method(Mesquite::QualityMetric::RMS);
    Mesquite::QualityAssessor stop_qa = Mesquite::QualityAssessor(&shape_metric);

    stop_qa.add_quality_assessment(&lapl_met);
    
    Mesquite::TerminationCriterion sc2;
    sc2.add_iteration_limit(10);
    lapl1.set_outer_termination_criterion(&sc2);

    
    q1.add_quality_assessor(&stop_qa, error);
    q1.set_master_quality_improver(&lapl1, error);
    q1.add_quality_assessor(&stop_qa, error);
    
    

    q1.run_instructions(&mesq_mesh, error);
     * */
    //shape_wrapper.run_instructions(&mesq_mesh,error);
    optmizer.run_instructions(&mesq_mesh,error);
    if (error)
    {
        std::cout << " Error smoothing mesh: " << std::endl << error <<  std::endl;
    }
    
}

#endif


void TETGenMeshing(Gmsh &gmsh, const char * gts_file);

/*
 * 
 */

int main(int argc, char** argv) {

    
    PetscInitialize(&argc, &argv, 0,0);
    
    Mesh          M;
    Gmsh     gmsh(M);
    
    
    EdgeCFDWriter ecfd(M);
    FluidStructure mesh_moviment(M);
   
    gmsh.read("bacia_se_nova.msh");
    M.extract_bnd_nodes();
    M.compute_quality();
    
    mesh_moviment.set_disp_from_mapping("nova_base.gts",MOVING, 0);
    mesh_moviment.add_dirichlet_surface_value(PNULL,0,0.0);
    mesh_moviment.solve();
     
    ecfd.write("teste");
    gmsh.save("test");
    printf("Done.\n");
    
    PetscFinalize();
  
    return 0;
}





