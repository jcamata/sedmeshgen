
#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <iomanip>

#include "mesh.h"
#include "edgecfd_writer.h"


using namespace std;


EdgeCFDWriter::EdgeCFDWriter(Mesh &M)
{
    mesh = &M;
}


void EdgeCFDWriter::write_mesh(const char *fname)
{
     ofstream    MshFile;
     double      xyz[3];
     long        conn[4];
     
     cout << " Writting EdgeCFD .msh file... " << endl;
     MshFile.open(fname, fstream::out);
     MshFile.flags(ios::right | ios::scientific);
     MshFile.precision(10);
     int nnodes = mesh->get_n_nodes();
     for(int node = 0; node < nnodes; node++){

         mesh->get_nodal_coordinate(node, xyz);
         MshFile << setw(20) << xyz[0]
                 << setw(20) << xyz[1]
                 << setw(20) << xyz[2] << endl;
     }
     int nel = mesh->get_n_elements();
     for(int iel=0; iel < nel; iel++)
     {
         mesh->get_element_conn(iel, conn);
         for(int i=0; i < 4 ; i++)
            MshFile << setw(10) << conn[i]+1;
         MshFile << setw(10) << 1 << endl;  
     }

   MshFile.close(); 
    
}

void EdgeCFDWriter::write_input(const char *fname)
{
 
    int   conn[4];
    double xyz[3];
    std::vector<int>& noslip  = mesh->getNodesOnBoundary(NOSLIP);
    std::vector<int>& mov     = mesh->getNodesOnBoundary(MOVING);
    std::vector<int>& slipx   = mesh->getNodesOnBoundary(SLIPX);
    std::vector<int>& slipy   = mesh->getNodesOnBoundary(SLIPY);
    std::vector<int>& slipz   = mesh->getNodesOnBoundary(SLIPZ);
    std::vector<int>& pnull   = mesh->getNodesOnBoundary(PNULL);
    std::vector<int>& solid   = mesh->getNodesOnBoundary(SOLID);
    std::vector<int>& outflowr   = mesh->getNodesOnBoundary(OUTFLOWR);
    std::vector<int>& outflowl   = mesh->getNodesOnBoundary(OUTFLOWL);
    std::vector<int>& outflowb   = mesh->getNodesOnBoundary(OUTFLOWB);
    std::vector<int>& inlet1     = mesh->getNodesOnBoundary(INLET1);
    std::vector<int>& inlet2     = mesh->getNodesOnBoundary(INLET2);
    
    //std::vector<int>& top    = mesh->getNodesOnBoundary(TOP);
    double vx,vy,vz;
    double s = 1.0;
//#define MAPA2    
#ifdef MAPA2
    vx  = 0.0990;
    vy  = 0.0;
    vz  = -0.0139;
#else
    
    vx  =  0.0706137*s;
    vy  = -0.0706144*s;
    vz  = -0.0052264*s;
#endif
    
    //std::set<int>   top;
    //std::vector<int> inteception;


    
    int ncar =  (noslip.size()+mov.size())*3 + 2*pnull.size() + 
                 outflowr.size()*3 + outflowl.size()*3 + outflowb.size()*3  +
                 inlet1.size()*4   + inlet2.size()*4;

   
    vector<int> sediment;
    vector<int> water;
 
    
    int lp = 0;
    int nic = 0;
    
   ofstream    InFile;
   InFile.open(fname, fstream::out);
   InFile.flags(ios::right | ios::scientific);
   InFile.precision(6);
  
   int nel    = mesh->get_n_elements();
   int nnodes = mesh->get_n_nodes(); 
   //LINE 1
   InFile << setw(40) << "DADOS PARA O EDGECFD" << endl;
   
   //LINE 2
   InFile << setw(10) << "4" //nnoel
          << setw(10) << nel //nel
          << setw(10) << "2" //numat
          << setw(10) << "3" //ndim
          << setw(10) << "4" //ngl
          << setw(10) << nnodes //nnos
          << setw(10) << "8" << endl;//npar
   
   //LINE 3
   InFile << setw(10) << "20" //maxbct
          << setw(10) << "2" //nbct
          << setw(10) << "1" //imet
          << setw(10) << "3" //isolver
          << setw(10) << "35" << endl; //kmax1
   
   //LINE 4        
   InFile << setw(10) << "10" //iprint
          << setw(10) << "7" //maxit
          << setw(10) << "-200" //maxslv
          << setw(10) << "0" //nansys
          << setw(10) << "1000" << endl; //ntan
   
   //LINE 5      
   InFile << setw(10) << "1" //nfun
          << setw(10) << ncar //ncar
          << setw(10) << nic << endl; //nic
   
   
   //LINE 6        
   InFile << setw(10) << "1.000e-4" //alpha
          << setw(10) << "0.100" //sigma0
          << setw(10) << "0.500" //sigma1
          << setw(10) << "0.900" //gamma
          << setw(10) << "0.100" //eta0
          << setw(10) << "1.000e-3" << endl; //etamax
   
   
   //LINE 7      
   InFile << setw(10) << "1.000e-3" //rtol
          << setw(10) << "1.000e-3" << endl; //dutol
   
   //LINE 8      
   InFile << setw(10) << "-0.050" //dt
          << setw(10) << "3.000" << endl; //tmax
   
   
   //LINE 9      
   InFile << setw(10) << "0.500" //cflmin
          << setw(10) << "1.000" //cflmax
          << setw(10) << "0.200" << endl; //tol_v
   
   //LINE 10      
   InFile << setw(10) << "2" //nptf
          << setw(10) << "CONSTANTE" << endl; //titfun(ifun)
   
   //TIME FUNCTION     
   InFile << setw(10) << "0.000e+0" //tf
          << setw(10) << "1.000" << endl; //val
        
   InFile << setw(10) << "1.000e+3" //tf
          << setw(10) << "1.000" << endl; //val
    
   {
        // NOSLIP BOUNDARY CONDITION
        vector<int>::iterator iter;
        
        for( iter = mov.begin(); iter != mov.end(); iter++)
       {
               int node = *(iter)+1;
               InFile <<setw(10)<< node << setw(7)  << -1 << setw(7)   << 1  << setw(7) << -1 << setw(15) << 0.0 << endl; lp++;
               InFile <<setw(10)<< node << setw(7)  << -2 << setw(7)   << 1  << setw(7) << -1 << setw(15) << 0.0 << endl; lp++;
               InFile <<setw(10)<< node << setw(7)  << -3 << setw(7)   << 1  << setw(7) << -1 << setw(15) << 0.0 << endl; lp++;
       }
        
        for( iter = noslip.begin(); iter != noslip.end(); iter++)
       {
               int node = *(iter)+1;
               InFile <<setw(10)<< node << setw(7)  << -1 << setw(7)   << 1  << setw(7) << 0 << setw(15) << 0.0 << endl; lp++;
               InFile <<setw(10)<< node << setw(7)  << -2 << setw(7)   << 1  << setw(7) << 0 << setw(15) << 0.0 << endl; lp++;
               InFile <<setw(10)<< node << setw(7)  << -3 << setw(7)   << 1  << setw(7) << 0 << setw(15) << 0.0 << endl; lp++;
       }
        
       for( iter = outflowb.begin(); iter != outflowb.end(); iter++)
       {
               int node = *(iter)+1;
               InFile <<setw(10)<< node << setw(7)  << -1 << setw(7)   << 1  << setw(7) << 0 << setw(15) <<   0.0     << endl; lp++;
               InFile <<setw(10)<< node << setw(7)  << -2 << setw(7)   << 1  << setw(7) << 0 << setw(15) <<  -0.0021*s  << endl; lp++;
               InFile <<setw(10)<< node << setw(7)  << -3 << setw(7)   << 1  << setw(7) << 0 << setw(15) <<   0.0      << endl; lp++;
       }
        
       for( iter = outflowl.begin(); iter !=outflowl.end(); iter++)
       {
               int node = *(iter)+1;
               InFile <<setw(10)<< node << setw(7)  << -1 << setw(7)   << 1  << setw(7) << 0 << setw(15) <<  -0.00210*s << endl; lp++;
               InFile <<setw(10)<< node << setw(7)  << -2 << setw(7)   << 1  << setw(7) << 0 << setw(15) <<  0.0      << endl; lp++;
               InFile <<setw(10)<< node << setw(7)  << -3 << setw(7)   << 1  << setw(7) << 0 << setw(15) <<  0.0       << endl; lp++;
       }
                    
        for( iter = outflowr.begin(); iter != outflowr.end(); iter++)
       {
               int node = *(iter)+1;
               InFile <<setw(10)<< node << setw(7)  << -1 << setw(7)   << 1  << setw(7) << 0 << setw(15) <<  0.0021*s    << endl; lp++;
               InFile <<setw(10)<< node << setw(7)  << -2 << setw(7)   << 1  << setw(7) << 0 << setw(15) <<  0.0       << endl; lp++;
               InFile <<setw(10)<< node << setw(7)  << -3 << setw(7)   << 1  << setw(7) << 0 << setw(15) <<  0.0       << endl; lp++;
       }
       
       
  

       for( iter = pnull.begin(); iter != pnull.end(); iter++)
       {
               int node = *(iter)+1;
               InFile <<setw(10)<< node << setw(7)  << -3 << setw(7)   << 1  << setw(7) << 0 << setw(15) << 0.0 << endl; lp++;
               InFile <<setw(10)<< node << setw(7)  << -4 << setw(7)   << 1  << setw(7) << 0 << setw(15) << 0.0 << endl; lp++;
       }
        
       //InFile         <<setw(10)<< 9    << setw(7)  << -4    << setw(7) << 1  << setw(7) << 0 << setw(15) << 0.0 << endl; lp++;
       //InFile         <<setw(10)<< 10   << setw(7)  << -4    << setw(7) << 1  << setw(7) << 0 << setw(15) << 0.0 << endl; lp++; 

       for( iter = inlet1.begin(); iter != inlet1.end(); iter++)
       {
               int node = *(iter)+1;
               InFile <<setw(10)<< node << setw(7)  << -1 << setw(7)   << 1  << setw(7) << 0 << setw(15) <<   vx   << endl; lp++;
               InFile <<setw(10)<< node << setw(7)  << -2 << setw(7)   << 1  << setw(7) << 0 << setw(15) <<   vy  << endl; lp++;
               InFile <<setw(10)<< node << setw(7)  << -3 << setw(7)   << 1  << setw(7) << 0 << setw(15) <<   vz   << endl; lp++;
               InFile <<setw(10)<< node << setw(7)  << -20 << setw(7)   << 1  << setw(7) << 0 << setw(15) <<   1.0      << endl; lp++;
       }
        
       // velocity: (0.08,-0.06,-0.01) canal 
       // velocity: 0.0990, 0, -0.0139
        
        
       for( iter = inlet2.begin(); iter != inlet2.end(); iter++)
       {
               int node = *(iter)+1;
               InFile <<setw(10)<< node << setw(7)  << -1 << setw(7)   << 1  << setw(7) << 0 << setw(15) <<   vx   << endl; lp++;
               InFile <<setw(10)<< node << setw(7)  << -2 << setw(7)   << 1  << setw(7) << 0 << setw(15) <<   vy  << endl; lp++;
               InFile <<setw(10)<< node << setw(7)  << -3 << setw(7)   << 1  << setw(7) << 0 << setw(15) <<   vz   << endl; lp++;
               InFile <<setw(10)<< node << setw(7)  << -20 << setw(7)   << 1  << setw(7) << 0 << setw(15) <<  0.0   << endl; lp++;
       }
        

        
   }


    
    
    InFile.close();
    
}


void EdgeCFDWriter::write(const char *rootname)
{
    char meshfile[80];
    char inputfile[80];
    
    sprintf(meshfile, "%s_ecfd.msh", rootname);
    sprintf(inputfile, "%s_ecfd.in", rootname);
    
    write_mesh(meshfile);
    write_input(inputfile);
    
}
