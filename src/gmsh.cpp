/*
 * gmsh.cpp
 *
 *  Created on: Mar 5, 2013
 *      Author: camata
 */

#include "gmsh.h"

#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <string.h>
using namespace std;

GmshElement::GmshElement()
{
    this->element_id = -1;
    this->element_type = 0;
    this->elementary = -1;
    this->physicalID = -1;
    this->num_tags = 0;
    this->num_nodes = 0;
}

int GmshElement::getPhysicalID()
{
    return this->physicalID;
}

void GmshElement::read(fstream* efile, int elemtype, int nnodes, int ntags)
{
    int i,j;
    this->element_type = elemtype;
    this->num_nodes    = nnodes;
    this->num_tags     = ntags;

    this->nodeIDs.resize(this->num_nodes);
    this->tags.resize(this->num_tags);

    // reading element id from gmsh file
    efile->read( (char*)&this->element_id, sizeof(int));
    for(i=0; i < this->num_tags; i++)
        efile->read((char*)&this->tags[i], sizeof(int));
    for(i=0; i < this->num_nodes; i++)
        efile->read((char*)&this->nodeIDs[i], sizeof(int));

    if(this->num_tags >= 2 )
    {
        this->physicalID = this->tags[0];
        this->elementary = this->tags[1];
    }

}

int GmshElement::getNumNodes()
{
    return this->num_nodes;
}

int GmshElement::getNumTags()
{
    return this->num_tags;
}

int GmshElement::node(int ith)
{
    return this->nodeIDs[ith];
}

int GmshElement::tag(int ith)
{
    return this->tags[ith];
}

int GmshElement::type(){
    return this->element_type;
}

GmshElement::~GmshElement()
{
    this->nodeIDs.clear();
    this->tags.clear();
}


GmshElement& GmshElement::operator=(const GmshElement& gelem)
{
    this->element_id   = gelem.element_id;
    this->element_type = gelem.element_type;
    this->elementary   = gelem.elementary;
    //this->nodeIDs     = gelem.nodeIDs;
    this->num_nodes    = gelem.num_nodes;
    this->num_tags     = gelem.num_tags;
    this->physicalID   = gelem.physicalID;
    //this->tags         = gelem.tags;
    return *this;
}


std::vector<int>& GmshElement::getConnectivity()
{
    return this->nodeIDs;
}



int Gmsh::read(const char* filename)
{
    std::vector<int> nodeset;
    GmshPhysicalParameter     id;
    int nsed, nfluidh, nfluidl, nfluid;
    
    nsed = nfluidh = nfluidl = nfluid = 0;
    
    this->readFile(filename);
    this->print_info();
    
    
    GmshPhysicalParameter sed, fluid, fluidh, fluidl;
    
    getPhysicalIDByName("SEDIMENT",sed);
    getPhysicalIDByName("FLUID",fluid);
    getPhysicalIDByName("FLUID_HEAVY",fluidh);
    getPhysicalIDByName("FLUID_LIGHT",fluidl);
          
    _mesh->set_n_nodes(this->n_nodes);
    for(int i=0; i < this->n_nodes; i++)
    {
        
        _mesh->set_nodal_coordinate(i, this->getX(i), this->getY(i), this->getZ(i));
    }

    _mesh->set_n_element(this->n_elements);
    int iel =0;
    for(int i=this->n_faces; i < this->n_total_elements; i++)
    {
        std::vector<int> conn = this->elements[i].getConnectivity();
	_mesh->set_element_conn(iel, conn[0]-1, conn[1]-1, conn[2]-1, conn[3]-1);
        int id = this->elements[i].getPhysicalID();
        
        if( id == sed.physical_number )
        {
            _mesh->set_element_tag(iel, SEDIMENT);
            nsed++;
        }
        
        if( id == fluid.physical_number )
        {
            _mesh->set_element_tag(iel, FLUID);
            nfluid++;
        }
        
        if( id == fluidh.physical_number )
        {
            _mesh->set_element_tag(iel, FLUID_HEAVY);
            nfluidh++;
        }
        
        if( id == fluidl.physical_number )
        {
            _mesh->set_element_tag(iel, FLUID_LIGHT);
            nfluidl++;
        }
       
        
        iel++;
     }
    
    cout << "# Elements on sediment region: "    << nsed << endl;
    cout << "# Elements on fluid heavy region: " << nfluidh << endl;
    cout << "# Elements on fluid light region: " << nfluidl << endl;
    cout << "# Elements on fluid region: "       << nfluid << endl;
     cout<< "-----------------------------------------------------------" << endl;
     cout<< " Boundary Nodes: " << endl;
    
     // set nodal tag
    if (getPhysicalIDByName("NOSLIP",id))
    {
	getNodesMarkedByPhysicalID(id,nodeset);
        
	for(int i = 0; i < nodeset.size(); i++)
	{
            int no_id = nodeset[i];
            _mesh->set_nodal_tag(no_id, NOSLIP);
	}
        cout<< " NOSLIP wall has " << nodeset.size() << " nodes." << endl;
        
    }

    if (getPhysicalIDByName("SLIPX",id))
    {
		getNodesMarkedByPhysicalID(id,nodeset);
		for(int i = 0; i < nodeset.size(); i++)
		{
			int no_id = nodeset[i];
			_mesh->set_nodal_tag(no_id, SLIPX);
		}
                cout<< " SLIPX wall has " << nodeset.size() << " nodes." << endl;
    }

    if (getPhysicalIDByName("SLIPY",id))
    {
		getNodesMarkedByPhysicalID(id,nodeset);
		for(int i = 0; i < nodeset.size(); i++)
		{
			int no_id = nodeset[i];
			_mesh->set_nodal_tag(no_id, SLIPY);
		
                }
                cout<< " SLIPY wall has " << nodeset.size() << " nodes." << endl;
    }

    if (getPhysicalIDByName("SLIPZ",id))
    {
		getNodesMarkedByPhysicalID(id,nodeset);
		for(int i = 0; i < nodeset.size(); i++)
		{
			int no_id = nodeset[i];
			_mesh->set_nodal_tag(no_id, SLIPZ);
		}
                cout<< " SLIPZ wall has " << nodeset.size() << " nodes." << endl;
     }

     if (getPhysicalIDByName("PNULL",id))
     {
		getNodesMarkedByPhysicalID(id,nodeset);
		for(int i = 0; i < nodeset.size(); i++)
		{
			int no_id = nodeset[i];
			_mesh->set_nodal_tag(no_id, PNULL);
		}
                cout<< " PNULL wall has " << nodeset.size() << " nodes." << endl;
    }
 
     if (getPhysicalIDByName("MOVING",id))
    {
                getNodesMarkedByPhysicalID(id,nodeset);
		for(int i = 0; i < nodeset.size(); i++)
		{
			int no_id = nodeset[i];
			_mesh->set_nodal_tag(no_id, MOVING);
		}
                cout<< " MOVING wall has " << nodeset.size() << " nodes." << endl;
    }
     
    if (getPhysicalIDByName("SOLID",id))
    {
                getNodesMarkedByPhysicalID(id,nodeset);
		for(int i = 0; i < nodeset.size(); i++)
		{
			int no_id = nodeset[i];
			_mesh->set_nodal_tag(no_id, SOLID);
		}
                cout<< " SOLID wall has " << nodeset.size() << " nodes." << endl;
    }
     
    if (getPhysicalIDByName("OUTFLOWL",id))
    {
                getNodesMarkedByPhysicalID(id,nodeset);
		for(int i = 0; i < nodeset.size(); i++)
		{
			int no_id = nodeset[i];
			_mesh->set_nodal_tag(no_id, OUTFLOWL);
		}
                cout<< " OUTFLOWL wall has " << nodeset.size() << " nodes." << endl;
    }
     
    if (getPhysicalIDByName("OUTFLOWR",id))
    {
                getNodesMarkedByPhysicalID(id,nodeset);
		for(int i = 0; i < nodeset.size(); i++)
		{
			int no_id = nodeset[i];
			_mesh->set_nodal_tag(no_id, OUTFLOWR);
		}
                cout<< " OUTFLOWR wall has " << nodeset.size() << " nodes." << endl;
    }
     
    if (getPhysicalIDByName("OUTFLOWB",id))
    {
                getNodesMarkedByPhysicalID(id,nodeset);
		for(int i = 0; i < nodeset.size(); i++)
		{
			int no_id = nodeset[i];
			_mesh->set_nodal_tag(no_id, OUTFLOWB);
		}
                cout<< " OUTFLOWB wall has " << nodeset.size() << " nodes." << endl;
    }
     
    if (getPhysicalIDByName("INLET1",id))
    {
                getNodesMarkedByPhysicalID(id,nodeset);
		for(int i = 0; i < nodeset.size(); i++)
		{
			int no_id = nodeset[i];
			_mesh->set_nodal_tag(no_id, INLET1);
		}
                cout<< " INLET1 wall has " << nodeset.size() << " nodes." << endl;
    }
     
    if (getPhysicalIDByName("INLET2",id))
    {
                getNodesMarkedByPhysicalID(id,nodeset);
		for(int i = 0; i < nodeset.size(); i++)
		{
			int no_id = nodeset[i];
			_mesh->set_nodal_tag(no_id, INLET2);
		}
                cout<< " INLET2 wall has " << nodeset.size() << " nodes." << endl;
    }
     
    return 0;
}




int Gmsh::readFile(const char* filename)
{
    string beginSection;
    this->gmsh_file.open(filename, ios::in | ios::binary);
    if(!this->gmsh_file.is_open())
    {
       cerr<<"ERROR: GMSH file, "<< string(filename)
           <<", cannot be opened. Does it exist? Have you read permission?\n";
       exit(-1);
    }

    read_header();
    this->gmsh_file >> beginSection;
    // Skip newline character
    char newLineChar[2];
    this->gmsh_file.read( newLineChar, 1 );

    if (beginSection == string("$PhysicalNames"))
    {
       read_physical_names();

       this->gmsh_file >> beginSection;

    }
    else if( beginSection != "$Nodes" ) {
       cerr << "Can't find $Nodes tag\n";
       exit(1);
    }

    read_nodes_bin();
    read_elements_bin();

    this->gmsh_file.close();
    return 0;

}

void Gmsh::print_info()
{
           //01234567890123456789012345678901234567890123456789
    cout << "-----------------------------------------------------------" << endl
         << "  GmshReader: Information Board                   " << endl
         << "-----------------------------------------------------------" << endl
         << "n-dim      = " << this->n_dim        << endl
         << "n-elements = " << this->n_elements   << endl
         << "n-faces    = " << this->n_faces      << endl
         << "n-nodes    = " << this->n_nodes      << endl
         << "-----------------------------------------------------------" << endl;
     cout<< " Physical Parameters: " << endl;
     for(int i = 0; i < this->PhysicalNames.size(); i++)
     {
       cout <<  " dim = " <<  this->PhysicalNames[i].physical_dim
            <<  "  number =  "   << this->PhysicalNames[i].physical_number
            <<  "  " <<  this->PhysicalNames[i].physical_name << endl;
     }
     cout << "----------------------------------------------------------" << endl;
}





void Gmsh::getElementMarkedByPhysicalID(GmshPhysicalParameter &id, std::vector<int> &elem)
{
   // int start;
    int num_elem = 0;
    if( id.physical_number < 0 || id.physical_number > this->PhysicalNames.size() )
    {
        cout << "Invalid Physical parameter " << endl;
        return;
    }
    
    elem.clear();

    int start= this->n_faces;
    int end  = this->n_total_elements;
    
    if (id.physical_dim == 2 && this->n_faces != 0 ) return;
    else{
        start = 0;
        end   = this->n_total_elements;
    }

    for(int i = start; i < end; i++)
    {
	if(this->elements[i].getPhysicalID() == id.physical_number)
        {
            elem.push_back(i);
        }
    }

}

void Gmsh::getFacesMarkedByPhysicalID(GmshPhysicalParameter &id, std::vector<int> &elem)
{

    int start;
    if(  id.physical_number < 0 || id.physical_number > this->PhysicalNames.size() )
    {
        cout << "Invalid Physical parameter " << endl;
        return;
    }

    if (id.physical_dim == 3 ) return;
    
    int n_faces = this->n_faces;
    if(this->n_dim == 2 && this->n_faces == 0) n_faces = this->n_elements;
    elem.clear();

    for(int i = 0; i < n_faces; i++)
    {
        if(this->elements[i].getPhysicalID() == id.physical_number)
        {
            elem.push_back(i);
        }
    }

}
void Gmsh::getNodesMarkedByPhysicalID(GmshPhysicalParameter &id, std::vector<int>& nodeSet)
{
    vector<int>::iterator it;
    
    nodeSet.clear();

    if(id.physical_number < 0 || id.physical_number > this->PhysicalNames.size())
    {
        cout << "Invalid Physical parameter " << endl;
        return;
    }

	if( id.physical_dim == 2 )
	{
            int n_faces = this->n_faces;
            if(this->n_dim == 2 && this->n_faces == 0) n_faces = this->n_elements;
   
		for(int i = 0; i < n_faces; i++)
		   if(this->elements[i].getPhysicalID() == id.physical_number)
                        for(int j=0; j < this->elements[i].getNumNodes(); j++)
                            nodeSet.push_back(this->elements[i].node(j));
	}
    else
	{

		for(int i = this->n_faces; i < this->n_total_elements; i++)
			if(this->elements[i].getPhysicalID() == id.physical_number)
				for(int j=0; j < this->elements[i].getNumNodes(); j++)
					nodeSet.push_back(this->elements[i].node(j));
    }

    sort(nodeSet.begin(),nodeSet.end());
    it = unique(nodeSet.begin(),nodeSet.end());
    nodeSet.resize(it-nodeSet.begin());

}


std::vector<int>& Gmsh::getElementConnectivity(int iel)
{
	return this->elements[iel+this->n_faces].getConnectivity();
}

std::vector<int>& Gmsh::getFaceConnectivity(int iel)
{
	return this->elements[iel].getConnectivity();
}


double Gmsh::getX(int node_id)
{
    return this->coords[node_id*3];
}


void Gmsh::setX(int node_id, double x)
{
  this->coords[node_id*3] = x;   
}

double Gmsh::getY(int node_id)
{
    return this->coords[node_id*3+1];
}

void Gmsh::setY(int node_id, double y)
{
  this->coords[node_id*3+1] = y;   
}

double Gmsh::getZ(int node_id)
{
    return this->coords[node_id*3+2];
}

void Gmsh::setZ(int node_id, double z)
{
  this->coords[node_id*3+2] = z;   
}


void Gmsh::read_physical_names()
{
    int num_physical;
    char sectionTag[80];
    char newLineChar[2]="\n";
    this->gmsh_file >> num_physical ;
    this->gmsh_file.read( newLineChar, 1 );
    this->PhysicalNames.resize(num_physical);
    for(int i = 0; i < num_physical; i++) {
         this->gmsh_file >> this->PhysicalNames[i].physical_dim
                         >> this->PhysicalNames[i].physical_number
                         >> this->PhysicalNames[i].physical_name ;
         this->gmsh_file.read( newLineChar, 1 );
    }

    this->gmsh_file.getline( sectionTag, sizeof("$EndPhysicalNames") );
    this->gmsh_file.read( newLineChar, 1 );

}


void Gmsh::read_header()
{
    int errorCode=1;
    string beginTag;
    int doubleFloatLen;

    this->gmsh_file >> beginTag >> this->gmsh_version >> this->gmsh_format >> doubleFloatLen;

     if( beginTag != "$MeshFormat" ) {
        cerr << "Can't find $MeshFormat tag\n";
        exit(errorCode);
     }

    if( this->gmsh_version < 2 || this->gmsh_version >= 3 ) {
        cerr << "Currently only GMSH format v2.x is supported\n";
        exit(errorCode);
    }

    if( doubleFloatLen != sizeof(double) ) {
        cerr << "Double float size in GMSH file differs from system double\n";
        exit(errorCode);
    }

    // For performance reasons only binary GMSH files are permitted with fldecomp
    if( this->gmsh_version == 0 ) {
       cerr << "** GMSH ASCII files are verboten:\n"
            << "** please use 'gmsh -bin' to generate binary format.\n";
       exit(errorCode);
    }

    // Skip newline character
    char newLineChar[2];
    this->gmsh_file.read( newLineChar, 1 );

     // Read in integer '1' written in binary format: check for endianness
     int oneInt;
     this->gmsh_file.read( (char *)&oneInt, sizeof(int) );
     this->gmsh_file.read( newLineChar, 1 );


     /*if( oneInt != 1 ){
        cerr << "** Oh dear: internal and file binary formats for integers differ (endianness?)\n";
        exit(errorCode);
     }*/


     char sectionTag[80];
     this->gmsh_file.getline( sectionTag, sizeof("$EndMeshFormat") );

     // Check for end of formatting section
     if( string(sectionTag) != "$EndMeshFormat" ) {
         cerr << "Can't find $EndMeshFormat tag\n";
         exit(errorCode);
     }

}

void Gmsh::read_nodes_bin()
{
    int errorCode=1;
    char sectionTag[80];
    char newLineChar[2]="\n";

    // Start of nodes section
    //

    this->gmsh_file >> this->n_nodes;
    this->gmsh_file.read( newLineChar, 1 );

    // Read in binary node data
    this->nodeIDs.resize(this->n_nodes);
    this->coords.resize(this->n_nodes*3);

    for(int n=0; n<this->n_nodes; n++)
    {
      this->gmsh_file.read( (char *)&this->nodeIDs[n], sizeof(int) );
      this->gmsh_file.read( (char *)&this->coords[n*3], sizeof(double) );
      this->gmsh_file.read( (char *)&this->coords[n*3+1], sizeof(double) );
      this->gmsh_file.read( (char *)&this->coords[n*3+2], sizeof(double) );
    }

    this->gmsh_file.read( newLineChar, 1 );


    // End of nodes section
    this->gmsh_file.getline( sectionTag, sizeof("$EndNodes") );

    if( string(sectionTag) != "$EndNodes" )
    {
        cerr << "Can't find $EndNodes tag\n";
        exit(errorCode);
    }

}


void Gmsh::read_elements_bin()
{
    int errorCode=1;
    char sectionTag[80];
    char newLineChar[2]="\n";
    int elemType;
    int numBlockElems;
    int numTags;

    // Elements section
    this->gmsh_file.getline( sectionTag, sizeof("$Elements") );


    if( string(sectionTag) != "$Elements" ){
        cerr << "Can't find $Elements tag\n";
        exit(errorCode);
    }

    this->gmsh_file >> this->n_total_elements;

    this->gmsh_file.read( newLineChar, 1 );

    int eIx=0;

    // Element array
    this->elements.resize(this->n_total_elements);

    int numEdges=0, numTriangles=0, numQuads=0, numTets=0, numHexes=0;

    while( eIx < this->n_total_elements )
    {
        this->gmsh_file.read( (char *)&elemType, sizeof(int) );
        this->gmsh_file.read( (char *)&numBlockElems, sizeof(int) );
        this->gmsh_file.read( (char *)&numTags, sizeof(int) );

        // Do we know about this type of element?
        if( elemNumNodes[elemType] == -1 )
        {
            cerr << "Element type not supported by fldecomp\n";
            exit(errorCode);
        }

        // Read in all the elements of this type
        for(int e=0; e<numBlockElems; e++)
        {
            this->elements[eIx].read(&this->gmsh_file, elemType, elemNumNodes[elemType], numTags);

            // Here, we count up the number of different types of elements.
            // This allows us to decide what are faces, and what are internal
            // elements.
          switch( elemType )
            {
            case 1:
              numEdges++;
              break;
            case 2:
              numTriangles++;
              break;
            case 3:
              numQuads++;
              break;
            case 4:
              numTets++;
              break;
            case 5:
              numHexes++;
              break;
            case 15:
              break;
            default:
              cerr << "Unsupported element type in GMSH mesh\n";
              cerr << "type: "<< elemType << "\n";
              exit(errorCode);
              break;
            }

          eIx++;
        }
    }

  this->gmsh_file.read( newLineChar, 1 );

  // End of Elements section
  this->gmsh_file.getline( sectionTag, sizeof("$EndElements") );


  if( string(sectionTag) != "$EndElements" ) {
      cerr << "Can't find $EndElements tag\n";
      exit(errorCode);
    }

  // Make some calculations based on the different types of elements
  // collected.

  if(numTets>0) {
    this->n_faces   = numTriangles;
    this->face_type = 2;
    this->elem_type = 4;
    this->n_dim     = 3;

  } else if(numTriangles>0) {
    this->n_faces   = numEdges;
    this->face_type = 1;
    this->elem_type = 2;
    this->n_dim  = 2;

  } else if(numHexes>0) {
    this->n_faces   = numQuads;
    this->face_type = 3;
    this->elem_type = 5;
    this->n_dim = 3;

  } else if(numQuads>0) {
    this->n_faces = numEdges;
    this->face_type = 1;
    this->elem_type = 3;
    this->n_dim = 2;

  } else {
       cerr << "Unsupported mixture of face/element types\n";
       exit(errorCode);

  }

  // Set some handy variables to be used elsewhere
  this->n_elements        = this->n_total_elements - this->n_faces;
  this->n_noel           = elemNumNodes[this->elem_type];
  this->n_nodes_per_face = elemNumNodes[this->face_type];

}



bool Gmsh::getPhysicalIDByName (const char* PhysicalName, GmshPhysicalParameter& id){


    for(int i = 0; i < this->PhysicalNames.size(); i++) {

       if (  this->PhysicalNames[i].physical_name == ("\"" + string(PhysicalName) + "\"") ){

            id = this->PhysicalNames[i];
			return true;

       }

     }

     return false;

}


bool Gmsh::getPhysicalIDByNumber (int id, GmshPhysicalParameter& p){


    for(int i = 0; i < this->PhysicalNames.size(); i++) {

       if (  this->PhysicalNames[i].physical_number == id ){

            p = this->PhysicalNames[i];
      return true;

       }

     }

     return false;

}

int  Gmsh::getNumberofNodes(){
     return this->n_nodes;
}

int  Gmsh::getNumberofElements(){
     return this->n_elements;
}

int  Gmsh::getNumberofFaces(){
     return this->n_faces;
}

Gmsh::~Gmsh()
{
    nodeIDs.clear();
    coords.clear();
    elements.clear();
    PhysicalNames.clear();
}


void Gmsh::save(const char *filename)
{
    string buffer;
    FILE *fout;
    fout = fopen(filename, "w");
    if(!fout)
    {
       cerr<<"ERROR: GMSH file, "<< string(filename)
           <<", cannot be opened. Does it exist? Have you read permission?\n";
       exit(-1);
    }
    
    printf("Saving mesh in gmsh format..\n");
    fprintf(fout,"$MeshFormat\n");
    fprintf(fout,"2.2 0 %d\n",sizeof(double));
    //int one = 1;
    //fwrite(&one, sizeof(int), 1, fout);
    fprintf(fout,"$EndMeshFormat\n");
    fprintf(fout,"$PhysicalNames\n");
    fprintf(fout, "%d\n", this->PhysicalNames.size());
    for(int i = 0; i < this->PhysicalNames.size(); i++)
    {
        fprintf(fout,"%d %d %s\n", this->PhysicalNames[i].physical_dim,
                                   this->PhysicalNames[i].physical_number,
                                   this->PhysicalNames[i].physical_name.c_str());
        
    }
    fprintf(fout,"$EndPhysicalNames\n");
    fprintf(fout,"$Nodes\n");
    int nnodes = this->n_nodes; //this->_mesh->get_n_nodes();
    //fwrite(&nnodes,sizeof(int),1,fout);
    fprintf(fout,"%d\n", nnodes);
    for(int i =0; i < nnodes ;++i)
    {
        int id = i+1;
        double xyz[3];
        //xyz[0] = getX(i);
        //xyz[1] = getY(i);
        //xyz[2] = getZ(i);
        this->_mesh->get_nodal_coordinate(i,xyz);
        //fwrite(&id, sizeof(int)   , 1, fout);
        //fwrite(xyz, sizeof(double), 3, fout);
        fprintf(fout,"%d %f %f %f\n",id, xyz[0], xyz[1], xyz[2]);
    }
    fprintf(fout,"$EndNodes\n");
    fprintf(fout,"$Elements\n");
    
    int nelem = this->n_total_elements;
    //fwrite(&nelem,sizeof(int),1,fout);
    fprintf(fout,"%d\n", nelem);
    for(int i = 0; i < nelem; i++)
    {
        int header[3];
        header[0] = i+1;
        header[1] = this->elements[i].type();
        header[2] = this->elements[i].getNumTags();
        //fwrite(header, sizeof(int), 3, fout);
        fprintf(fout,"%d %d %d ", header[0], header[1], header[2]);
        for(int j = 0; j < this->elements[i].getNumTags(); j++)
            fprintf(fout, "%d ", this->elements[i].tag(j));
        vector<int> conn = this->elements[i].getConnectivity();
        for(int j = 0; j < this->elements[i].getNumNodes(); j++)
            fprintf(fout, "%d ", conn[j]);
        
        fprintf(fout, "\n");
    }
    fprintf(fout,"$EndElements\n");
    fclose(fout);
}


int Gmsh::GetNumberofPhysicalParameters() 
{
    return this->PhysicalNames.size();
}

GmshPhysicalParameter* Gmsh::getPhysicalParameter(int i)
{
    return &this->PhysicalNames[i];
}

GmshElement* Gmsh::getElement(int i)
{
    return &this->elements[i];
}


