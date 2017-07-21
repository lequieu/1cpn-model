#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "string.h"
#include "math_vector.h"
//#include "myrandom.h"

#define TWID 500
//typedef double vector[3];                // 0:x  1:y  2:z
//typedef double quaternion[4];                // quaternion
//typedef double form[6];                        // 0:xx 1:yy 2:zz 3:zy 4:zx 5:yx
//inline void quat_vec_rot(vector &dest, vector &src, quaternion &q) {
//  // dest = q*src*conj(q)
//  quaternion aa={q[0]*q[0], q[1]*q[1], q[2]*q[2], q[3]*q[3]};
//  form ab={q[0]*q[1], q[0]*q[2], q[0]*q[3], q[1]*q[2], q[1]*q[3], q[2]*q[3]};
//  dest[0] = (aa[0]+aa[1]-aa[2]-aa[3])*src[0]+
//            ((ab[3]-ab[2])*src[1]+(ab[1]+ab[4])*src[2])*2.0;
//  dest[1] = (aa[0]-aa[1]+aa[2]-aa[3])*src[1]+
//            ((ab[2]+ab[3])*src[0]+(ab[5]-ab[0])*src[2])*2.0;
//  dest[2] = (aa[0]-aa[1]-aa[2]+aa[3])*src[2]+
//            ((ab[4]-ab[1])*src[0]+(ab[0]+ab[5])*src[1])*2.0;
//}
using namespace LAMMPS_NS;
int main(int argc, char**argv){
    
   if (argc != 4){
     printf("Usage: %s <dump file> <new xyz file> <length of u vector>\n",argv[0]);
     exit(1);
   }


   char dumpfile[TWID], newdumpfile[TWID];
   sprintf(dumpfile,"%s",argv[1]);
   sprintf(newdumpfile,"%s",argv[2]);
   double length;
   length = atof(argv[3]);
   
   double fvu0[3][3];
   fvu0[0][0] = 1; fvu0[0][1] = 0; fvu0[0][2] = 0; 
   fvu0[1][0] = 0; fvu0[1][1] = 1; fvu0[1][2] = 0; 
   fvu0[2][0] = 0; fvu0[2][1] = 0; fvu0[2][2] = 1; 

   //myrandom_set_seed(12354); 
   //myrandom_set_seed(12); 
  
   FILE *fptr, *ofptr;
   char subline[TWID],junk[TWID];
   int timestep, natom;
   double xboxmin, xboxmax, yboxmin, yboxmax, zboxmin, zboxmax;
   double **atoms, **quats;
   int *atom_types;
   int i,j;

   bool writexyz = true ;
   bool writedump = false;
   if (writexyz && writedump){
        printf("Error! cannot have both xyz and dump flags set to true!\n");
        exit(1);
   }

   bool firstframe = true;
   double *coloring;

   fptr = fopen(dumpfile,"r");
   ofptr = fopen(newdumpfile,"w");

   if (fptr == NULL) {printf("Error! %s does not exist!\n",dumpfile); exit(1);}
   if (ofptr == NULL) {printf("Error! %s does not exist!\n",newdumpfile); exit(1);}
   while (fscanf(fptr,"%s",subline) != EOF){
     //fscanf(fptr,"%s",subline);
     if (!strcmp(subline,"ITEM:")){
       fscanf(fptr,"%s",subline);
       //fgets(subline,TWID,fptr);
       if (!strcmp(subline,"TIMESTEP")){
         fgets(junk,TWID,fptr);
         fscanf(fptr,"%d",&timestep);
       }
       else if (!strcmp(subline,"NUMBER")){
         fgets(junk,TWID,fptr);
         fscanf(fptr,"%d",&natom);
       }
       else if (!strcmp(subline,"BOX")){
         fgets(junk,TWID,fptr);
         fscanf(fptr,"%lf",&xboxmin);
         fscanf(fptr,"%lf",&xboxmax);
         fgets(junk,TWID,fptr);
         fscanf(fptr,"%lf",&yboxmin);
         fscanf(fptr,"%lf",&yboxmax);
         fgets(junk,TWID,fptr);
         fscanf(fptr,"%lf",&zboxmin);
         fscanf(fptr,"%lf",&zboxmax);
       }
       else if (!strcmp(subline,"ATOMS")){
         fgets(junk,TWID,fptr);
         //printf("Timestep: %d, Natom: %d, Box: %f %f %f %f %f %f\n",timestep,natom,xboxmin, xboxmax, yboxmin, yboxmax,zboxmin, zboxmax);
         
         //allocate memory
         atoms = (double**)calloc(natom,sizeof(double*));
         quats = (double**)calloc(natom,sizeof(double*));
         for(i=0;i<natom;i++){
            atoms[i] = (double*)calloc(3,sizeof(double));
            quats[i] = (double*)calloc(4,sizeof(double*));
         }
         atom_types = (int*) calloc(natom,sizeof(int));
        
         //read atom positions
         int idx; 
         for(i=0;i<natom;i++){
            fscanf(fptr,"%d",&idx);
            idx -= 1;
            fscanf(fptr,"%d",&atom_types[idx]);
            atom_types[idx] -= 1;
            fscanf(fptr,"%lf",&atoms[idx][0]);
            fscanf(fptr,"%lf",&atoms[idx][1]);
            fscanf(fptr,"%lf",&atoms[idx][2]);

            fscanf(fptr,"%lf",&quats[idx][0]);
            fscanf(fptr,"%lf",&quats[idx][1]);
            fscanf(fptr,"%lf",&quats[idx][2]);
            fscanf(fptr,"%lf",&quats[idx][3]);
            fgets(junk,TWID,fptr);
         }

         //write new dump file
         int nnewatom;
         int type;

         nnewatom = natom * (3+1);
         
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

         
         //for writing xyz 
         if (writexyz){
           fprintf(ofptr,"%d\n\n",nnewatom);
         }

         // for writing to dump
         if (writedump){
           fprintf(ofptr,"ITEM: TIMESTEP\n%d\n",timestep);
           fprintf(ofptr,"ITEM: NUMBER OF ATOMS\n%d\n",nnewatom);
           fprintf(ofptr,"ITEM: BOX BOUNDS\n%f %f\n%f %f\n%f %f\n",xboxmin, xboxmax, yboxmin, yboxmax, zboxmin, zboxmax);
           fprintf(ofptr,"ITEM: ATOMS id type x y z vx vy vz\n");
         }
         double angle;
         double x,y,z;
         double fvu[3][3];
         int rotaxis, initaxis;
         int iatom=0;
         vector vin, vout;
         for(i=0;i<natom;i++){
            type = atom_types[i];

            //print central atom
            x = atoms[i][0];
            y = atoms[i][1];
            z = atoms[i][2];
            if (writedump){
              fprintf(ofptr,"%d %d %f %f %f %f %f %f\n",iatom + 1,type + 1,x,y,z,coloring[iatom], 0.0,0.0); //dump
            }
            if (writexyz){
              fprintf(ofptr,"%d %f %f %f\n",type + 1,x,y,z); //xyz
            }
            iatom ++;

            quaternion rquat;
            double norm;
            norm = sqrt(quats[i][0]*quats[i][0] +  quats[i][1]*quats[i][1] + quats[i][2]*quats[i][2] + quats[i][3]*quats[i][3]);
            if (fabs(norm - 1.0) > 1e-3){
                printf("Norm of quat is too big (%f)! You've screwed something up with quat math!\n",norm);
                exit(1);
            }
            
            //rquat[0] = cos(0.5*quats[i][0]);
            //rquat[1] = quats[i][1]*sin(0.5*quats[i][0]);
            //rquat[2] = quats[i][2]*sin(0.5*quats[i][0]);
            //rquat[3] = quats[i][3]*sin(0.5*quats[i][0]);
            rquat[0] = quats[i][0];
            rquat[1] = quats[i][1];
            rquat[2] = quats[i][2];
            rquat[3] = quats[i][3];

            
            //get f,v,u
            quat_vec_rot(fvu[0],fvu0[0],rquat);          
            quat_vec_rot(fvu[1],fvu0[1],rquat);          
            quat_vec_rot(fvu[2],fvu0[2],rquat);          

            for (j=0; j<3; j++){
                
                x = atoms[i][0] + fvu[j][0] * length;
                y = atoms[i][1] + fvu[j][1] * length;
                z = atoms[i][2] + fvu[j][2] * length;

                if (writedump){
                  //fprintf(ofptr,"%d %d %f %f %f %f %f %f\n",iatom+1,type+j+3,x,y,z,coloring[iatom], 0.0,0.0); //dump
                  if (j==0) fprintf(ofptr,"%d %s %f %f %f %f %f %f\n",iatom+1,"F",x,y,z,coloring[iatom], 0.0,0.0); //dump
                  if (j==1) fprintf(ofptr,"%d %s %f %f %f %f %f %f\n",iatom+1,"V",x,y,z,coloring[iatom], 0.0,0.0); //dump
                  if (j==2) fprintf(ofptr,"%d %s %f %f %f %f %f %f\n",iatom+1,"U",x,y,z,coloring[iatom], 0.0,0.0); //dump
                }
                if (writexyz){
                  //fprintf(ofptr,"%d %f %f %f\n",type + 3 + j,x,y,z); //xyz
                  if (j==0) fprintf(ofptr,"%s %f %f %f\n","F",x,y,z); //xyz
                  if (j==1) fprintf(ofptr,"%s %f %f %f\n","V",x,y,z); //xyz
                  if (j==2) fprintf(ofptr,"%s %f %f %f\n","U",x,y,z); //xyz
                }
                iatom ++;
            }
         }


         //free memory
         free(atom_types);
         for(i=0;i<natom;i++){
           free(atoms[i]);
         }
       }
        
     }
   }
   
}
