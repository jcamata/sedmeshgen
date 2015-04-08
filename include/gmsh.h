/*
 * gmsh.h
 *
 *  Created on: Mar 5, 2013
 *      Author: camata
 */

#ifndef GMSH_H_
#define GMSH_H_

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cstring>

#include "mesh.h"
#include "define.h"

using namespace std;


static int elemNumNodes[] = {2,2,3,4,4,8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};

typedef struct
{
    int    physical_dim;
    int    physical_number;
    std::string physical_name;
} GmshPhysicalParameter;


class GmshElement
{
  public:
      GmshElement();
      void read(fstream *efile, int elemtype, int nnodes, int ntags);
      //int  getElementNumTags();
      //int  getElementNodes();
      int  getNumNodes();
      int  getNumTags();
      std::vector<int>& getConnectivity();
      int node(int ith);
      int tag(int ith);
      int type();
      int getPhysicalID();
      GmshElement& operator=(const GmshElement& gelem);
      ~GmshElement();

  private:
      int element_id;
      int element_type;
      int num_tags;
      int num_nodes;
      int physicalID;
      int elementary;
      vector<int>   nodeIDs;
      vector<int>   tags;

};




class Gmsh
{
    public:
	  Gmsh(Mesh &M);
          Gmsh(const Gmsh& GR);
          void save (const char* filename);
	  //GmshReader(const Mesh &M);
          int read(const char* file);
          void print_info();
          bool getPhysicalIDByName(const char* pname,GmshPhysicalParameter &id);
          bool getPhysicalIDByNumber(int id, GmshPhysicalParameter &p);
          void getElementMarkedByPhysicalID(GmshPhysicalParameter &id, std::vector<int> &elem);
          void getFacesMarkedByPhysicalID(GmshPhysicalParameter &id, std::vector<int> &elem);
          void getNodesMarkedByPhysicalID(GmshPhysicalParameter &id, std::vector<int>& nodeSet);
          std::vector<int>& getElementConnectivity(int iel);
          std::vector<int>& getFaceConnectivity(int iel);
          double getX(int node_id);
          double getY(int node_id);
          double getZ(int node_id);
          void   setX(int node, double x);
          void   setY(int node, double y);
          void   setZ(int node, double z);
          int    getNumberofNodes();
          int    getNumberofElements();
          int    getNumberofFaces();
          int    GetNumberofPhysicalParameters();
          GmshElement* getElement(int i);
          GmshPhysicalParameter* getPhysicalParameter(int i);
          ~Gmsh();

    private:
      Mesh *_mesh;

      int readFile(const char* filename);
      // private methods
      void read_header();
      void read_nodes_bin();
      void read_elements_bin();
      void read_physical_names();

      // private attributes
      fstream                       gmsh_file;
      float                         gmsh_version;
      int                           gmsh_format;
      int                           n_nodes;
      int                           n_elements;
      vector<int>                   nodeIDs;
      vector<double>                coords;
      vector<GmshElement>           elements;
      vector<GmshPhysicalParameter> PhysicalNames;
      int                           face_type;
      int                           elem_type;
      int                           n_faces;
      int                           n_total_elements;
      int                           n_dim;
      int                           n_noel;
      int                           n_nodes_per_face;
};


inline Gmsh::Gmsh(Mesh &M)
{
  this->_mesh = &M;
  this->n_faces = 0;
  this->face_type = 0;
  this->elem_type = 0;
  this->n_dim = 0;
  this->n_noel = 0;
  this->n_nodes = 0;
  this->n_elements = 0;
  this->n_total_elements = 0;
  this->n_nodes_per_face = 0;
  this->gmsh_format = 0;
  this->gmsh_version = 0.0;

}


#endif /* GMSH_H_ */
