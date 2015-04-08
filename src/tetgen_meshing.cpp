

#include "gmsh.h"
#ifdef _TETGEN

#include "tetgen.h"
#include <gts.h>
#include <glib.h>

#include <cstdio>
#include <cstdlib>
#include <vector>
using namespace std;

gdouble  distance(GtsPoint *p, gpointer bounded);

void TETGenMeshing(Gmsh &gmsh, const char * file)
{
    
    double c[3];
    std::vector<int> nodeset;

    FILE      *gts_file;
    GtsSurface *s;
    GtsPoint   *p;
    GtsFile    *fp;
    GNode      *t;
    
    gts_file = fopen(file, "r");
    if(!gts_file) return ;
    fp = gts_file_new(gts_file);
    s  = gts_surface_new(gts_surface_class (),
		         gts_face_class (),
		         gts_edge_class (),
		         gts_vertex_class ());
    if (gts_surface_read (s, fp)) {
        fputs ("file on standard input is not a valid GTS file\n", stderr);
        fprintf (stderr, "stdin:%d:%d: %s\n", fp->line, fp->pos, fp->error);
        return; /* failure */
    }
    
    t = gts_bb_tree_surface(s);
    p = gts_point_new(gts_point_class(), 0.0 , 0.0 ,0.0);
    
    double factor =1000.0;
    GmshPhysicalParameter ph;
    gmsh.getPhysicalIDByName("MOVING",ph);
    gmsh.getNodesMarkedByPhysicalID(ph, nodeset);
    
    std:: cout << nodeset.size() << " nodes located in bottom wall" << endl;
    
    for(int i=0; i < nodeset.size(); i++)
    {
        int id = nodeset[i];
        
        c[0] = gmsh.getX(id);
        c[1] = gmsh.getY(id);
        c[2] = gmsh.getZ(id);
        if(c[2] > 1.0E-3) continue;
         
        c[0]*=factor;
        c[1]*=factor;
        c[2]*=factor;
        
        

        gts_point_set(p, c[0], c[1], c[2]);
        double d = gts_bb_tree_point_distance(t,p,distance,NULL);
        d = (4.6*factor - d)/factor;
        //d = 4.6 - d;
        gmsh.setZ(id, d); 
 
    }
   
    gts_bb_tree_destroy(t, true);   
    
    gmsh.save("moving.msh");
    
    
    cout << "calling tetgen..." << endl;
    tetgenio in, out;
    tetgenio::facet *f;
    tetgenio::polygon *pol;
    
    in.firstnumber    = 1;
    in.numberofpoints = gmsh.getNumberofNodes();
    in.pointlist      = new REAL (in.numberofpoints*3);
    for(int i =0; i < gmsh.getNumberofNodes(); i++)
    {
        in.pointlist[i*3]   = gmsh.getX(i);
        in.pointlist[i*3+1] = gmsh.getY(i);
        in.pointlist[i*3+2] = gmsh.getZ(i);
    }
    
    in.numberoffacets  = gmsh.GetNumberofPhysicalParameters();
    in.facetlist       = new tetgenio::facet[in.numberoffacets];
    in.facetmarkerlist = new int[in.numberoffacets];
    
    for(int i = 0; i < gmsh.GetNumberofPhysicalParameters(); i++)
    {
        GmshPhysicalParameter *pp = gmsh.getPhysicalParameter(i);
        
        in.facetmarkerlist[i] = pp->physical_number;
        f = &in.facetlist[i];
        
        std::vector<int> elmap;
        gmsh.getFacesMarkedByPhysicalID(*pp,elmap);
        
        f->numberofpolygons = elmap.size();
        f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
        f->numberofholes = 0;
        f->holelist = NULL;
   
        for(int j = 0; j < elmap.size(); j++ ) {
            GmshElement *elem = gmsh.getElement(elmap[j]);
            pol = &f->polygonlist[j];
 
            pol->numberofvertices = elem->getNumNodes();
            pol->vertexlist = new int[pol->numberofvertices];
            for(int k=0; k <  pol->numberofvertices; k++)
                pol->vertexlist[k] = elem->node(k);
        }
        
    }
    
  tetrahedralize("pq1.414a0.1", &in, &out);

  // Output mesh to files 'barout.node', 'barout.ele' and 'barout.face'.
  out.save_nodes("barout");
  out.save_elements("barout");
  out.save_faces("barout");
    
    
}

#else
void TETGenMeshing(Gmsh &gmsh, const char * file)
{
	return;
}
#endif


