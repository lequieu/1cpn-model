#include <stdlib.h>
#include <math.h>
//#include "myrandom.h"
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include "trajectory_iterator.h"

#define TWID 500

class Bond{
    public:
      int stea;
      int steb;
      Bond(){};
      Bond(int a, int b):stea(a),steb(b){};
};

using namespace LAMMPS_NS;
int main(int argc, char**argv){
    
    if (argc != 3) {
        std::cout<<"Usage: "<<argv[0]<<" <dump file> <new file prefix (will write xyz and psf)>"<<std::endl;
        exit(1);
    }

    bool writexyz = true ;
    bool writedump = false;
    bool centercom = true;
    if (writexyz && writedump){
        printf("Error! cannot have both xyz and dump flags set to true!\n");
        exit(1);
    }

    char dumpfile[TWID], newdumpfile[TWID], newpsffile[TWID];
    sprintf(dumpfile,"%s",argv[1]);
    if (writedump)
        sprintf(newdumpfile,"%s.dump",argv[2]);
    if (writexyz)
        sprintf(newdumpfile,"%s.xyz",argv[2]);
        sprintf(newpsffile,"%s.psf",argv[2]);

    int natomtypes = 5;
    double r[natomtypes],a[natomtypes],d[natomtypes], nrot[natomtypes],c[natomtypes];
    int n[natomtypes];
    //r = radius of `cylinder` of extra sites
    //n = number of points that makeup the cylinder
    //a = radius of sites that compose cylinder, currently not used
    //d = stem height
    double axis[natomtypes][3];

    double ls = 55.0 / 1.0; //length scale
    //nucleosome
    r[0] = 45.37; //38.5; //0.7 * ls; //25;
    a[0] = 1100; //20 * ls;
    n[0] = 50; //16;
    d[0] = 32.370058; // from 1kx5 //31.0 * 1.5 ;  //stem height only for nucl
    c[0] = 80.0; //stem length

    //stem_interp_angle is gamma from init_1cpn for 11bp unwrap
    //double stem_interp_angle = 0.628318530718; 
    double stem_interp_angle = 90 * M_PI / 180.; //where does the DNA start to leave the surface of the nucleosome
    nrot[0] = 2;
    axis[0][0] = 1; axis[0][1] = 0; axis[0][2] = 0; //must be unit vector!

    //dna
    r[1] = 0.06 * ls;//5;
    a[1] = 5 * ls;
    axis[1][0] = 1; axis[1][1] = 0; axis[1][2] = 0; //must be unit vector!
    n[1] = 1;
    d[1] = 0 * ls;
    c[1] = 0. * ls;
    nrot[1] = 0;

    //ghost
    r[2] = 0.06 * ls;//5;
    a[2] = 5 * ls;
    axis[2][0] = 1; axis[2][1] = 0; axis[2][2] = 0; //must be unit vector!
    n[2] = 0;
    d[2] = 0 * ls;
    c[2] = 0. * ls;
    nrot[2] = 0;

    //linker histone haven't decided how to viz yet
    for (size_t lh_iter=3; lh_iter<natomtypes; lh_iter++) {
        r[lh_iter] = 0.00;
        a[lh_iter] = 0;
        n[lh_iter] = 0;
        d[lh_iter] = 0;
        c[lh_iter] = 0.;
        nrot[lh_iter] = 0;
    }

    //MASSES are hardcoded!
    std::vector<double> masses(5+1);
    masses[1] = 196666.0000;
    masses[2] = 1950.000000;  
    masses[3] = 19500.00000;  
    masses[4] = 1466.000000;
    masses[5] = 550.0000000;



    //myrandom_set_seed(12); 

    FILE *ofptr, *opsfptr;
    int timestep, natom;
    double *coloring;

    //Set up all vectors needed for the trajectory parser class
    TrajectoryIterator parser;
    std::vector<Bond> bonds;
    std::vector<int> atom_types;
    std::vector<float> box_dim;
    std::vector<std::vector<double>> atoms;
    std::vector<std::vector<double>> quats;
    std::vector<std::vector<double>> vects_f;
    std::vector<std::vector<double>> vects_u;
       
    //Load the trajectory into the parser
    parser.load_dump(argv[1]); 
    //Get the number of frames
    timestep = parser.get_numFrames();
    std::cout<<timestep<<std::endl;
    //Get the vector for the types of atoms
    atom_types = parser.get_type(); 
    //Get the number of atoms
    natom = parser.get_numAtoms();
    box_dim = parser.get_boxDim();

    //Open the output files (old code ideals once this is integrated into the parser, think about cleaning up and converting to c++)
    ofptr = fopen(newdumpfile,"w");
    opsfptr = fopen(newpsffile,"w");
    if (ofptr == NULL) {printf("Error! %s does not exist!\n",newdumpfile); exit(1);}
    if (opsfptr == NULL) {printf("Error! %s does not exist!\n",newpsffile); exit(1);}

    bool firstframe = true;
    //Loop through the dump file using the parser
    for(size_t i=0; i<timestep; i++) {

        //Write new dump file
        int nnewatom;
        int type,prevtype;
        std::vector<int> natomoftype(natomtypes);

        //The actual functions from the parser
        atoms = parser.get_coord();
        quats = parser.get_quat();
        vects_f = parser.get_vect(quats,'f');
        vects_u = parser.get_vect(quats,'u');
        //Make sure to move to the next frame
        //Since everything after this is processing the data we can put the next frame option here
        parser.next_frame();

        for (size_t k=0;k<natom;k++){
            type = atom_types[k]-1 ;
            if (type > natomtypes){
                printf("ERROR: more than two atom types in dump! FIX ME!!!\n");
                exit(1);
            }
            natomoftype[type]++;
        }

        nnewatom = 0;
        for (size_t k=0;k<natomtypes;k++){
            nnewatom += natomoftype[k]*(n[k] + 1);
        }

        //calculate COM
        double com[3], totalmass;
        com[0] =0.0;
        com[1] =0.0;
        com[2] =0.0;
        totalmass = 0;
        for (size_t iatom=0; iatom < natom; iatom++){
          type = atom_types[iatom]-1;
          com[0] += atoms[iatom][0]*masses[type];
          com[1] += atoms[iatom][1]*masses[type];
          com[2] += atoms[iatom][2]*masses[type];
          totalmass += masses[type];
        }
        com[0] /= totalmass;
        com[1] /= totalmass;
        com[2] /= totalmass;


        //Write colors for visualization
        if (firstframe){
            coloring = (double*) calloc(nnewatom,sizeof(double));
        //    double ransum,ran;
        //    ransum = 0;
            for (size_t j=0;j<nnewatom; j++){
        //        ran = myrandom_gauss(1.0);
        //        ransum += ran;
        //        coloring[j] = ransum;
                coloring[j] = 0.0;
            }
        }

        //Write psf
        if (firstframe){
            fprintf(opsfptr,"*\n*\n*\n*\n*\n\n");
            fprintf(opsfptr,"\t%d !NATOMS\n",nnewatom);
            bonds.reserve(2*nnewatom);            
        }
 
        //For writing xyz 
        if (writexyz){
            fprintf(ofptr,"%d\n\n",nnewatom);
        }

        //For writing to dump
        if (writedump){
            fprintf(ofptr,"ITEM: TIMESTEP\n%d\n",timestep);
            fprintf(ofptr,"ITEM: NUMBER OF ATOMS\n%d\n",nnewatom);
            fprintf(ofptr,"ITEM: BOX BOUNDS\n%f %f\n%f %f\n%f %f\n",box_dim[0],box_dim[1],box_dim[2],box_dim[3],box_dim[4],box_dim[5]);
            fprintf(ofptr,"ITEM: ATOMS id type x y z vx vy vz\n");
        }
        double angle;
        double x,y,z;
        int rotaxis, initaxis;
        int iatom=0;
        int ibond=0;
        char name[TWID];
        vector vin, vout;
        for(size_t k=0;k<natom;k++) {
            type = atom_types[k]-1;

            if (type==0) sprintf(name,"N");
            else if (type==1) sprintf(name,"D");
            else if (type==2) sprintf(name,"G");
            else if (type==3) sprintf(name,"GH");
            else if (type==4) sprintf(name,"CTD");

            //print central atom
            x = atoms[k][0];
            y = atoms[k][1];
            z = atoms[k][2];

            if (centercom){
              x -= com[0];
              y -= com[1];
              z -= com[2];
            }

            if (writedump){
              fprintf(ofptr,"%d %d %f %f %f %f %f %f\n",iatom + 1,type + 1,x,y,z,coloring[iatom], 0.0,0.0); //dump
            }
            if (writexyz){
              fprintf(ofptr,"%s %f %f %f\n",name,x,y,z); //xyz
            }

            if (firstframe) {//write psf

                fprintf(opsfptr,"%8d %s %-4s %s  %-5s  %3d %13.6e  %7.3f            0\n",iatom+1, "MOL1", "RESI", "RESN", name,1,0.0, 100.0);

                if (k >= 1) {
                    //write dna-dna bonds
                    if ((atom_types[k] == 2) && (atom_types[k-1] == 2)){
                        bonds.push_back(Bond(iatom,iatom-1-n[1]));
                        ibond++;
                    }
                }
            }
            iatom ++;

            double rotaxis[3];
            double initaxis[3];
            if (type == 0){ //nucleosome
                for(size_t j = 0; j < 3; j++) {
                    rotaxis[j] = vects_f[k][j];
                    initaxis[j] = vects_u[k][j];
                }
            }
            else if (type == 1){ //dna
                for(size_t j = 0; j < 3; j++) {
                    rotaxis[j] = vects_u[k][j];
                    initaxis[j] = vects_f[k][j];
                }
            }

            vin[0] = initaxis[0];
            vin[1] = initaxis[1];
            vin[2] = initaxis[2];

            double dstem[3];
            double cstem[3];
            
            int rotcount = 0;
            int nsite4stem;
            int siteperrot = n[type] / nrot[type];
            double factor; //stem length scaling factor
            //double anglepersite = (2*M_PI * nrot[type]) / n[type]; //this was when stem was explicitly drawn
            double anglepersite;
            
            if (n[type] > 1 )
              anglepersite = (2*M_PI * nrot[type] - 2*stem_interp_angle) / (n[type]-1);
            else 
              anglepersite = (2*M_PI * nrot[type] - 2*stem_interp_angle) / (n[type]);
            //nsite4stem = ceil(stem_interp_angle / anglepersite);

            if (type==0) sprintf(name,"N1");
            else if (type==1) sprintf(name,"D1");
            else if (type==2) sprintf(name,"G1");
            else if (type ==3) sprintf(name,"GH1");
            else if (type ==4) sprintf(name,"CTD1");

            for (size_t j=0; j<n[type]; j++){
                angle = anglepersite * j + stem_interp_angle; 
                  
                if ((j % siteperrot) == 0){
                   rotcount ++; 
                }

                quaternion rq;
                rq[0] = cos(0.5*angle);
                rq[1] = rotaxis[0] * sin(0.5*angle);
                rq[2] = rotaxis[1] * sin(0.5*angle);
                rq[3] = rotaxis[2] * sin(0.5*angle);
                quat_vec_rot(vout,vin,rq);          
                
                dstem[0] =  rotaxis[0] * d[type] * (0.5 -(double) j/n[type]);
                dstem[1] =  rotaxis[1] * d[type] * (0.5 -(double) j/n[type]);
                dstem[2] =  rotaxis[2] * d[type] * (0.5 -(double) j/n[type]);
                
                factor = 0.0;

                cstem[0] = factor * vin[0] * (c[type] - r[type]);
                cstem[1] = factor * vin[1] * (c[type] - r[type]);
                cstem[2] = factor * vin[2] * (c[type] - r[type]);

                x = atoms[k][0] + vout[0] * r[type] + dstem[0] + cstem[0];
                y = atoms[k][1] + vout[1] * r[type] + dstem[1] + cstem[1];
                z = atoms[k][2] + vout[2] * r[type] + dstem[2] + cstem[2];

                if (centercom){
                  x -= com[0];
                  y -= com[1];
                  z -= com[2];
                }

                if (writedump){
                    fprintf(ofptr,"%d %d %f %f %f %f %f %f\n",iatom+1,type+1+natomtypes,x,y,z,coloring[iatom], 0.0,0.0); //dump
                }
                if (writexyz){
                    fprintf(ofptr,"%s %f %f %f\n",name ,x,y,z); //xyz
                }

                if (firstframe){//write psf
                    fprintf(opsfptr,"%8d %s %-4s %s  %-5s  %3d %13.6e  %7.3f            0\n",iatom+1, "MOL1", "RESI", "RESN", name,1,0.0, 100.0);

                    if (j == 0){ //if first site
                        if ((type == 0) && (prevtype == 1)) { 
                            //for nucl, bond it to the last dna site
                            int bondsite = iatom-2-n[1]; //no ghost
                            //int bondsite = iatom-2-n[1]-1-n[2]; //ghost
                            if (bondsite >= 0){
                                bonds.push_back(Bond(bondsite,iatom));
                                ibond++;
                            }
                        }
                        if ((type==1) && (n[type] = 1) && (k>0) && prevtype==1){
                            //for DNA:  bond between the previous new site, this will help viz helix
                            bonds.push_back(Bond(iatom-1-n[1],iatom));
                            ibond++;
                        }
                    }
                    else{ //if not first site

                        if (type == 0){ //if type == nucl, then bond all the extra sites together
                            bonds.push_back(Bond(iatom-1,iatom));
                            ibond++;
                        }
                    }

                    if (j==(n[type]-1)){// also if last site
                        
                        if ((k != natom-3) && (type == 0) && (atom_types[k+2] == 2)){
                            //int bondsite = iatom + 1; //bond to next dna, no ghost
                            int bondsite = iatom + 1 + 1 + n[2]; //bond to next dna, with ghost
                            if (bondsite < nnewatom){
                                bonds.push_back(Bond(iatom,bondsite));
                                ibond++;
                            }
                        }
                    }
                }
                iatom ++;
            }
            prevtype = type;
        }

        //write bonds to psf
        if (firstframe){
            fprintf(opsfptr,"\n%8d !NBONDS\n",bonds.size());
            int count = 0;
            while (count < bonds.size()){
                for (int ii = 0; ii<4;ii++){
                    fprintf(opsfptr,"%8d%8d", bonds[count].stea+1, bonds[count].steb+1);
                    count++;
                    if (count >= bonds.size()) break;
                }
                fprintf(opsfptr,"\n");
            }
        }
        if (firstframe) firstframe = false;
    }  
}
