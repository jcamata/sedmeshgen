


#include "mesh.h"

#ifdef _MESQUITE

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




void Mesh::optimize()
{
    
    std::vector<int> fixed(this->_n_nodes);
    
    std::fill(fixed.begin(),fixed.end(), 0);
    
    for(int s = 0; s < this->_bnd_nodes.surfaces.size(); s++ )
    {
        std::vector<int> *surface = &this->_bnd_nodes.surfaces[s];
        for(std::vector<int>::iterator iter = surface->begin(); iter != surface->end(); iter++) 
        {
            int node_id = *iter;
            fixed[node_id] = 1;
        }
    }
    
    
    int dim = 3;
    unsigned long int nnodes = get_n_nodes();
    unsigned long int nelem  = get_n_elements();
    unsigned long   *conn    = &this->_elements[0];
    double          *coords  = &this->_xyz[0];
    int             *fptr    = &fixed[0];
    
    
   
    
    Mesquite::ArrayMesh mesq_mesh(dim, nnodes, coords, fptr, nelem, Mesquite::TETRAHEDRON, conn);

    Mesquite::MsqError error;

  
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
    //shape_wrapper.run_instructions(&mesq_mesh,error);
    //optmizer.run_instructions(&mesq_mesh,error);
    if (error)
    {
        std::cout << " Error smoothing mesh: " << std::endl << error <<  std::endl;
    }
   
}

#else
void Mesh::optimize()
{
	return;
}
#endif

