
#include "book.h"
#include "errhand.h"

#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>
#define NTHREAD   32 


#define CUDA_CALL(x) do { if((x) != cudaSuccess) { \
printf("Error at %s:%d\n",__FILE__,__LINE__); \
return EXIT_FAILURE;}} while(0)
#define CURAND_CALL(x) do { if((x)!=CURAND_STATUS_SUCCESS) { \
printf("Error at %s:%d\n",__FILE__,__LINE__);\
return EXIT_FAILURE;}} while(0)

/*  Declaration of variables that will be preserved accross kernel calls */

static int *Map_Cell_dev, *Offset_dev;
static int *Nl_Cell_dev, *Nat_Cell_dev, *Npatches_dev;
static int *Ntor_angle_dev, *List_Cell_dev;
static int *overlap_error, *overlap_error_dev;
static int *Ntras_ac_dev, *Nrot_ac_dev, *Ntras_ac, *Nrot_ac;
static int *Nat_Cell_new, *List_Cell_new;

static double *W_RU_dev;
static float *patchdev, *hdev, *DR_Cell_dev, *Sphere_Celldev;
static float *energydev, *range_dev, *XOP_dev, *VL0_dev, *sigma_LJ_dev;
static float *ref_tor_dev, *tor_angle_dev, *Vpot_Matrix_dev, *sigma_jon_dev;
static float *Delta_E_dev, *Delta_E, *energy;
static float *Sphere_Cell_new;


/* This kernel initializes state per thread */

__global__ void setup_random(int *Nl_Cell, 
		     unsigned long long seed,
		     unsigned long long offset, 
		     curandState_t *state)
{

  int Indx_Cell = threadIdx.x + blockIdx.x * blockDim.x;
  curand_init(seed, Indx_Cell , offset, &state[Indx_Cell]);
}

///////////////////////////////////////////////////////////////
__device__ void vlj( float dist, float sigma, float VL0, float& Vrad) {

  float R126, R1212, distr;

  distr = dist/sigma;
  R126=distr*distr*distr*distr*distr*distr;
  R1212=R126*R126;
  Vrad=-4.0*((R1212 - R126)/(R1212*R126) + VL0);

}


//////////////////////////////////////////////////////////////////////
__device__ void get_torsion(float PX, float PY, float PZ, float DX, float DY, float DZ,
                            float& vx, float& vy, float& vz, float& modv) {

   vx=PZ*DY-PY*DZ;
   vy=PX*DZ-PZ*DX;
   vz=PY*DX-PX*DY;
   modv=__fsqrt_rd(vx*vx+vy*vy+vz*vz);

}
//////////////////////////////////////////////////////////////////////
__device__ void get_torsion_angle(float v1x, float v1y, float v1z, float mod_v1,
                                  float v2x, float v2y, float v2z, float mod_v2,
                                  float DX, float DY, float DZ, float dist, float& torsion) {

      float cangle, v12x, v12y, v12z, v12b2;

      cangle=(v1x*v2x+v1y*v2y+v1z*v2z)/mod_v1/mod_v2;
      if(fabsf(cangle) > 1.0) {
         if( fabsf(cangle) - 1.0 < 0.0001 ) {
            cangle = copysignf(1.0,cangle); }
         else {
            printf("\n error in cangle %f",cangle);
         }
      }
      v12x=v1y*v2z-v1z*v2y;
      v12y=v1z*v2x-v1x*v2z;
      v12z=v1x*v2y-v1y*v2x;
      v12b2=(v12x*DX+v12y*DY+v12z*DZ)/dist/mod_v1/mod_v2;
      if( fabsf(v12b2) > 1.0) {
          if( fabsf(v12b2) - 1.0 < 0.0001 ) {
                v12b2 = copysignf(1.0,v12b2); }
          else {
                printf("\n error in v12b2 %f",v12b2);
          }
      }
      torsion= atan2(v12b2,cangle);

}

/////////////////////////////////////////////////////////////////////////////////
 __device__ void get_dist(float XJ,float YJ,float ZJ,float XI,float YI,float ZI,
                          float H11, float H12, float H13,
                          float H21, float H22, float H23,
                          float H31, float H32, float H33,
                          float& DX, float& DY, float& DZ, float& dist )  {

  float DXP, DYP, DZP;

  DX = XJ - XI ;
  DY = YJ - YI ;
  DZ = ZJ - ZI ;

  DXP = DX - nearbyintf(DX);
  DYP = DY - nearbyintf(DY);
  DZP = DZ - nearbyintf(DZ);

  DX = H11*DXP + H12*DYP + H13* DZP;
  DY = H21*DXP + H22*DYP + H23* DZP;
  DZ = H31*DXP + H32*DYP + H33* DZP;
  dist = __fsqrt_rd(DX*DX + DY*DY + DZ*DZ);

}

/////////////////////////////////////////////////////////
__device__ void rot_matrix(float Q0,  float Q1,  float Q2, float Q3, float R[9]) {

 float Q02, Q12, Q22, Q32, Q0Q1, Q0Q2, Q0Q3, Q1Q2, Q1Q3, Q2Q3;

  Q02 = Q0*Q0;
  Q12 = Q1*Q1;
  Q22 = Q2*Q2;
  Q32 = Q3*Q3;

  Q0Q1 = Q0*Q1;
  Q0Q2 = Q0*Q2;
  Q0Q3 = Q0*Q3;
  Q1Q2 = Q1*Q2;
  Q1Q3 = Q1*Q3;
  Q2Q3 = Q2*Q3;

  R[0] = Q02+Q12-Q22-Q32;
  R[1] = 2.00*(Q1Q2-Q0Q3);
  R[2] = 2.00*(Q1Q3+Q0Q2);
  R[3] = 2.00*(Q1Q2+Q0Q3);
  R[4] = Q02-Q12+Q22-Q32;
  R[5] = 2.00*(Q2Q3-Q0Q1);
  R[6] = 2.00*(Q1Q3-Q0Q2);
  R[7] = 2.00*(Q2Q3+Q0Q1);
  R[8] = Q02-Q12-Q22+Q32;

 }


//////////////////////////////////////////////////////////////////////
__device__ void patch_real ( float PX, float PY, float PZ,
                             float R[9], float& PX2, float& PY2, float& PZ2 ) {

 PX2 = R[0]*PX+R[1]*PY+R[2]*PZ;
 PY2 = R[3]*PX+R[4]*PY+R[5]*PZ;
 PZ2 = R[6]*PX+R[7]*PY+R[8]*PZ;

}

//////////////////////////////////////////////////////////////////////
__device__ void get_cosine(float CSOMG, float& OMG) {


       if( fabsf(CSOMG) > 1.0 ) {
          if ( fabsf(CSOMG) - 1.0 < 0.0001) {
              CSOMG = copysignf(1.0,CSOMG); }
          else {
              printf("\n error in csomg %f",CSOMG);
          }
       }
       OMG = acos(CSOMG);
}


/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////

__device__ void genranor(curandState_t *RND_state, int Indx_Rnd, float Rnd_or[4]) {

float xran1, xran2, xran3, xran4, ransq12, ransq34, ranh;

ransq12=2.;
do {
xran1 = 1. - 2.*curand_uniform(&RND_state[Indx_Rnd]); 
xran2 = 1. - 2.*curand_uniform(&RND_state[Indx_Rnd]); 
ransq12 = xran1*xran1 + xran2*xran2;
} while ( ransq12 >= 1. );


ransq34=2.;
do {
xran3 = 1. - 2.*curand_uniform(&RND_state[Indx_Rnd]); 
xran4 = 1. - 2.*curand_uniform(&RND_state[Indx_Rnd]); 
ransq34 = xran3*xran3 + xran4*xran4;
} while ( ransq34 >= 1. );

ranh = sqrtf( (1. - ransq12) / ransq34);

Rnd_or[0] = xran1;
Rnd_or[1] = xran2;
Rnd_or[2] = xran3*ranh;
Rnd_or[3] = xran4*ranh;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void enerGPU(int bool_tor, int Ndim, int Nat_Cell_Max, int Ntot_Cell, int Npatches_max, int Ntor_max, int Ntypes,
                        int *Npatches, int *Ntor_angle, float *tor_angle, float *patch, float *ref_tor, 
                        float *range, float *sigma_jon, float sigma_tor, float *sigma_LJ, float *VL0, float *XOP, 
                        float *Vpot_Matrix, 
                        float *h, float *Sphere_Cell, int *List_Cell, int *Map_Cell, int *Nl_Cell, int *Nat_Cell, 
                        float *DR_Cell, float *energy)
{

  int imol, jmol, kpatch, k;
  int Indx_Cell2, Indx_Cell3, l, m, indx_atm_aux;
  int Nat_Cell_ijk, Nat_Cell_neigh;
  int Nl0, Nl1, Nl2, ICell, indx_neigh_aux2;
  int imol2, jmol2, indx_atm_aux2, itypei, itypej;
  int Indx2, Indx3, Indx4, indx_neigh_auxp;
  int j3, j4, IndxT1, IndxT2, Indx;
  int jcell_x, jcell_y, jcell_z;

  float xi, yi, zi, energ, DX, DY, DZ, dist;
  float H11, H12, H13, H21, H22, H23, H31, H32, H33, DXP, DYP, DZP;
  float Q0i, Q1i, Q2i, Q3i;
  float Q102, Q112, Q122, Q132, Q10Q11, Q10Q12;
  float Q10Q13, Q11Q12, Q11Q13, Q12Q13, R111, R112, R113, R121, R122, R123, R131, R132, R133;
  float Q0j, Q1j, Q2j, Q3j, Q202, Q212, Q222, Q232, Q20Q21, Q20Q22;
  float Q20Q23, Q21Q22, Q21Q23, Q22Q23, R211, R212, R213, R221, R222, R223, R231, R232, R233;
  float CSM1M, CSM2M, PX, PY, PZ, PX1, PY1, PZ1, PX2, PY2, PZ2, CSOMG1, CSOMG2;
  float dx_neigh, dy_neigh, dz_neigh, factor_energy;
  float PX1_ref[20], PX2_ref[20];
  float PY1_ref[20], PY2_ref[20];
  float PZ1_ref[20], PZ2_ref[20];
  float COS1[20], COS2[20];
  float Vangtor_max, v1x, v1y, v1z, mod_v1, FX, FXI;
  float Vangtor, Vang, Vrad, v2x, v2y, v2z, mod_v2, torsion, OMG1, OMG2;
  float v12x, v12y, v12z, v12b2, Utor_max, Utor, cangle, R126, R1212;
  float dtor, sigma_jon_ij, s_LJ_aux, distr, XOP_aux, VL0_aux;
  float range_aux;
  float a =1000.0;
  float pi= acosf(-1.0), dospi=2.*pi;

  int Indx_Cell = threadIdx.x + blockIdx.x * blockDim.x;

  Nl0 = Nl_Cell[0];
  Nl1 = Nl_Cell[1];
  Nl2 = Nl_Cell[2];

  jcell_z = Indx_Cell/(Nl0*Nl1)  ;
  jcell_y = ((Indx_Cell - Nl0*Nl1*jcell_z)/Nl0)  ;
  if (jcell_y == 0 && jcell_z == 0) {
     jcell_x = Indx_Cell ;
 }
  else {
     jcell_x = ( Indx_Cell % (jcell_y*Nl0+jcell_z*Nl0*Nl1)) ;
  }

  if ( jcell_x < Nl0  && jcell_y < Nl1 && jcell_z < Nl2 ) {

  Nat_Cell_ijk = Nat_Cell[Indx_Cell];
  energ = -0.0000006;
  energy[Indx_Cell] = energ; // energ;

  if ( Nat_Cell_ijk > 0 ) {

    H11=h[0];
    H12=h[1];
    H13=h[2];
    H21=h[3];
    H22=h[4];
    H23=h[5];
    H31=h[6];
    H32=h[7];
    H33=h[8];

    indx_atm_aux = Indx_Cell * Nat_Cell_Max * 7;
    indx_atm_aux2 = Indx_Cell * Nat_Cell_Max * 2;

    for ( m=0; m<Nat_Cell_ijk; m++) {   // loop over all atoms in the chosen cell

        imol = indx_atm_aux + m*7;
        imol2 = indx_atm_aux2 + m*2;

        xi = Sphere_Cell[imol];
        yi = Sphere_Cell[imol+1];
        zi = Sphere_Cell[imol+2];

        Q0i = Sphere_Cell[imol+3];
        Q1i = Sphere_Cell[imol+4];
        Q2i = Sphere_Cell[imol+5];
        Q3i = Sphere_Cell[imol+6];

        itypei = List_Cell [ imol2 +1 ];

        Q102 = Q0i*Q0i;
        Q112 = Q1i*Q1i;
        Q122 = Q2i*Q2i;
        Q132 = Q3i*Q3i;

        Q10Q11 = Q0i*Q1i;
        Q10Q12 = Q0i*Q2i;
        Q10Q13 = Q0i*Q3i;
        Q11Q12 = Q1i*Q2i;
        Q11Q13 = Q1i*Q3i;
        Q12Q13 = Q2i*Q3i;

        R111 = Q102+Q112-Q122-Q132;
        R112 = 2.00*(Q11Q12-Q10Q13);
        R113 = 2.00*(Q11Q13+Q10Q12);
        R121 = 2.00*(Q11Q12+Q10Q13);
        R122 = Q102-Q112+Q122-Q132;
        R123 = 2.00*(Q12Q13-Q10Q11);
        R131 = 2.00*(Q11Q13-Q10Q12);
        R132 = 2.00*(Q12Q13+Q10Q11);
        R133 = Q102-Q112-Q122+Q132;

        // loop over the atoms in the same cell

        for ( l=0; l<Nat_Cell_ijk; l++) {   // loop over all atoms in the chosen cell
           if ( l != m) {

	    jmol = indx_atm_aux + l*7;
	    jmol2 = indx_atm_aux2 + l*2;

	    DX = Sphere_Cell[jmol] - xi ;   
	    DY = Sphere_Cell[jmol+1]- yi ;
	    DZ = Sphere_Cell[jmol+2] - zi ;
 
            itypej = List_Cell [ jmol2 +1 ];

            Indx = itypei * Ntypes + itypej;
            s_LJ_aux = sigma_LJ[ Indx ];
            range_aux = range[ Indx ];
            VL0_aux = VL0[ Indx ];
            XOP_aux = XOP[ Indx ];

	    DXP = DX - nearbyintf(DX);
	    DYP = DY - nearbyintf(DY);
	    DZP = DZ - nearbyintf(DZ);

	    DX = H11*DXP + H12*DYP + H13* DZP;
	    DY = H21*DXP + H22*DYP + H23* DZP;
	    DZ = H31*DXP + H32*DYP + H33* DZP;

	    dist = sqrtf(DX*DX + DY*DY + DZ*DZ );

	    if ( dist <= 0.5000000) {
	       energy[Indx_Cell] = 9.e+99;
	       return;
	    }
	    else if (dist < 0.98000000*s_LJ_aux ) {
                distr=dist/s_LJ_aux;
                R126=distr*distr*distr*distr*distr*distr;
                R1212=R126*R126;
                Vrad=-4.0*((R1212 - R126)/(R1212*R126) + VL0_aux);
	        energ=energ + Vrad ;
            }
	    else if (dist < range_aux ) {

	       Q0j = Sphere_Cell[jmol+3];
	       Q1j = Sphere_Cell[jmol+4];
	       Q2j = Sphere_Cell[jmol+5];
	       Q3j = Sphere_Cell[jmol+6];

	       Q202 = Q0j*Q0j;
	       Q212 = Q1j*Q1j;
	       Q222 = Q2j*Q2j;
	       Q232 = Q3j*Q3j;

	       Q20Q21 = Q0j*Q1j;
	       Q20Q22 = Q0j*Q2j;
	       Q20Q23 = Q0j*Q3j;
	       Q21Q22 = Q1j*Q2j;
	       Q21Q23 = Q1j*Q3j;
	       Q22Q23 = Q2j*Q3j;

	       R211 = Q202+Q212-Q222-Q232;
	       R212 = 2.00*(Q21Q22-Q20Q23);
	       R213 = 2.00*(Q21Q23+Q20Q22);
	       R221 = 2.00*(Q21Q22+Q20Q23);
	       R222 = Q202-Q212+Q222-Q232;
	       R223 = 2.00*(Q22Q23-Q20Q21);
	       R231 = 2.00*(Q21Q23-Q20Q22);
	       R232 = 2.00*(Q22Q23+Q20Q21);
	       R233 = Q202-Q212-Q222+Q232;

	       CSM1M=-100.;
	       CSM2M=-100.;


	       for ( k=0; k < Npatches[itypei]; k++) {

		  kpatch = Npatches_max*Ndim*itypei + k*Ndim;

		  PX = patch[kpatch];
		  PY = patch[kpatch+1];
		  PZ = patch[kpatch+2];

		  PX1 = R111*PX+R112*PY+R113*PZ;
		  PY1 = R121*PX+R122*PY+R123*PZ;
		  PZ1 = R131*PX+R132*PY+R133*PZ;

		  CSOMG1=(PX1*DX+PY1*DY+PZ1*DZ)/dist;
                  COS1[k] =  CSOMG1;
		  if(CSOMG1 > CSM1M) {
                         CSM1M = CSOMG1;
                  }

		  PX = ref_tor[kpatch];
		  PY = ref_tor[kpatch+1];
		  PZ = ref_tor[kpatch+2];

		  PX1_ref[k] = R111*PX+R112*PY+R113*PZ;
		  PY1_ref[k] = R121*PX+R122*PY+R123*PZ;
		  PZ1_ref[k] = R131*PX+R132*PY+R133*PZ;

	       } //for k < Npatches
 
	       for ( k=0; k < Npatches[itypej]; k++) {

		  kpatch = Npatches_max*Ndim*itypej + k*Ndim;

		  PX = patch[kpatch];
		  PY = patch[kpatch+1];
		  PZ = patch[kpatch+2];

		  PX2 = R211*PX+R212*PY+R213*PZ;
		  PY2 = R221*PX+R222*PY+R223*PZ;
		  PZ2 = R231*PX+R232*PY+R233*PZ;

		  CSOMG2=-(PX2*DX+PY2*DY+PZ2*DZ)/dist;
                  COS2[k] =  CSOMG2;
		  if(CSOMG2 > CSM2M) {
                         CSM2M = CSOMG2;
                  }

		  PX = ref_tor[kpatch];
		  PY = ref_tor[kpatch+1];
		  PZ = ref_tor[kpatch+2];

		  PX2_ref[k] = R211*PX+R212*PY+R213*PZ;
		  PY2_ref[k] = R221*PX+R222*PY+R223*PZ;
		  PZ2_ref[k] = R231*PX+R232*PY+R233*PZ;

	       } //for k < Npatches

               Vangtor_max = -1.0;
                
               for ( j3 = 0; j3 < Npatches[itypei]; j3++) {

                  v1x=PZ1_ref[j3]*DY-PY1_ref[j3]*DZ;
                  v1y=PX1_ref[j3]*DZ-PZ1_ref[j3]*DX;
                  v1z=PY1_ref[j3]*DX-PX1_ref[j3]*DY;

                  mod_v1=sqrtf(v1x*v1x+v1y*v1y+v1z*v1z);
                  CSOMG1 = COS1[j3];

                  if( fabsf(CSOMG1) > 1.0 ) {
                     if ( fabsf(CSOMG1) - 1.0 < 0.0001) {
                         CSOMG1 = copysignf(1.0,CSOMG1); }
                     else { 
                         printf("\n error in csomg1 %f",CSOMG1);
                     } 
                  }
                  OMG1 = acosf(CSOMG1);

                  for ( j4=0; j4 <  Npatches[itypej]; j4++) {

                      Indx2 = itypei * Npatches_max + j3 ;  // elijo fila
                      Indx3 = Indx2 * Npatches_max * Ntypes;  // pongo el indice al inicio de la fila
                      Indx4 = Indx3 + itypej * Npatches_max;  // me muevo en la fila hasta el atomo de tipo Itypej
                      factor_energy = Vpot_Matrix[ Indx4+ j4 ];

                      if (factor_energy > 0.000001 ) { 

                          CSOMG2 = COS2[j4];
                          if( fabsf(CSOMG2) > 1.0 ) {
                            if ( fabsf(CSOMG2) - 1.0 < 0.0001) {
                                CSOMG2 = copysignf(1.0,CSOMG2); }
                            else {
                                printf("\n error in csomg2 %f",CSOMG2);
                            } 
                          }
                          OMG2 = acosf(CSOMG2);

                          IndxT1 =  Npatches_max*itypei + j3;
                          IndxT2 =  Npatches_max*itypej + j4;

                          sigma_jon_ij = 2.0*sigma_jon[IndxT1]*sigma_jon[IndxT2];

                          Vang = __expf(-(OMG1*OMG1+OMG2*OMG2)/sigma_jon_ij);

                          if (bool_tor == 1 && Npatches[itypei] > 1 && Npatches[itypej] > 1 && Ntor_angle[IndxT1] > 0 && Ntor_angle[IndxT2] > 0) {

                              v2x=PZ2_ref[j4]*DY-PY2_ref[j4]*DZ;
                              v2y=PX2_ref[j4]*DZ-PZ2_ref[j4]*DX;
                              v2z=PY2_ref[j4]*DX-PX2_ref[j4]*DY;

                              mod_v2=sqrtf(v2x*v2x+v2y*v2y+v2z*v2z);

                              if( (mod_v1 != 0.0) &&  (mod_v2 != 0.0 ) ) {

                                    cangle=(v1x*v2x+v1y*v2y+v1z*v2z)/mod_v1/mod_v2;

                                    if(fabsf(cangle) > 1.0) {
                                       if( fabsf(cangle) - 1.0 < 0.0001 ) {
                                          cangle = copysignf(1.0,cangle); }
                                       else {
                                          printf("\n error in cangle %f",cangle);
                                       }     
                                    }  
                                    v12x=v1y*v2z-v1z*v2y;
                                    v12y=v1z*v2x-v1x*v2z;
                                    v12z=v1x*v2y-v1y*v2x;
                                    v12b2=(v12x*DX+v12y*DY+v12z*DZ)/dist/mod_v1/mod_v2;

                                    if( fabsf(v12b2) > 1.0) {
                                        if( fabsf(v12b2) - 1.0 < 0.0001 ) {   
                                              v12b2 = copysignf(1.0,v12b2); }
                                        else {
                                              printf("\n error in v12b2 %f",v12b2);
                                        }     
                                    }    

                                    torsion= atan2(v12b2,cangle);

                                    IndxT1 =  Npatches_max*itypei + j3;

                                    Utor_max = -100.0;

                                    for(k=0; k < Ntor_angle[ IndxT1 ]; k++) {

                                          dtor = torsion-tor_angle[IndxT1*Ntor_max+k] ;
                                          if (dtor > pi) {
                                              dtor = dtor - dospi; }
                                          else if (dtor < -pi) {
                                              dtor = dtor + dospi ;  }
                                          else if (abs(dtor) > dospi) {
                                                   printf("\n ERROR");
                                          }
                                          Utor = __expf (-dtor*dtor/sigma_tor);
                                          if (Utor  > Utor_max) Utor_max = Utor;

                                    }      

                                    Vangtor = factor_energy*Vang*Utor_max;
                              }
                              else {

                              } /* if( (mod_v1 != 0.0) &&  (mod_v2 != 0.0 )  */
                          }
                          else { /* if (bool_tor == 0)  No torsion */

                               Vangtor = factor_energy*Vang;
                                
                          }  /* if (bool_tor) */

                          if( Vangtor  > Vangtor_max)  { 
                                 Vangtor_max = Vangtor;   
                          }

                      } /* if (factor_energy > 0.00000001 ) */
                  } /* for ( j4=0; j4 <  Npatches[itypej]; j4++) */
               } /* for ( j3 = 0; j3 < Npatches[itypei]; j3++) */

     

               if(Vangtor_max == -1.0) {    
	            energ=energ + 0.0  ;
               }      
               else {
                   Vangtor=Vangtor_max;
                   FX=(1.0 + tanhf(a*(-XOP_aux + dist)))/2.0;
                   FXI=(1.0 - tanhf(a*(-XOP_aux + dist)))/2.0;
                   distr=dist/s_LJ_aux;
                   R126=distr*distr*distr*distr*distr*distr;
                   R1212=R126*R126;
                   Vrad=-4.0*((R1212 - R126)/(R1212*R126) + VL0_aux);
	           energ=energ + Vrad * (FXI + Vangtor * FX);
               }
	    }  // if dist < range
         }
       }


       for ( ICell= 0; ICell< 26; ICell++)  {      // loop over the neighbouring cells

	   Indx_Cell2 = Indx_Cell*26 + ICell ;
	   Indx_Cell2 = Map_Cell [Indx_Cell2]; 

	   Nat_Cell_neigh = Nat_Cell[Indx_Cell2];

	   Indx_Cell3 = Indx_Cell*26*3 + ICell*3 ;
	   dx_neigh = DR_Cell[Indx_Cell3];
	   dy_neigh = DR_Cell[Indx_Cell3+1];
	   dz_neigh = DR_Cell[Indx_Cell3+2];

	   if ( Nat_Cell_neigh > 0 ) {

	       indx_neigh_aux2 = Indx_Cell2 * Nat_Cell_Max * 7; 
	       indx_neigh_auxp = Indx_Cell2 * Nat_Cell_Max * 2; 

	       for (l=0; l< Nat_Cell_neigh; l++) {

	  	   jmol = indx_neigh_aux2 + l*7;
	  	   jmol2 = indx_neigh_auxp + l*2;

		   DX = Sphere_Cell[jmol] - xi + dx_neigh;  
		   DY = Sphere_Cell[jmol+1]- yi + dy_neigh;
		   DZ = Sphere_Cell[jmol+2] - zi + dz_neigh;

                   itypej = List_Cell [ jmol2 +1 ];

                   Indx = itypei * Ntypes + itypej;
                   s_LJ_aux = sigma_LJ[ Indx ];
                   range_aux = range[ Indx ];
                   VL0_aux = VL0[ Indx ];
                   XOP_aux = XOP[ Indx ];

		   DXP = DX - nearbyintf(DX);
		   DYP = DY - nearbyintf(DY);
		   DZP = DZ - nearbyintf(DZ);

		   DX = H11*DXP + H12*DYP + H13* DZP;
		   DY = H21*DXP + H22*DYP + H23* DZP;
		   DZ = H31*DXP + H32*DYP + H33* DZP;

		   dist = sqrtf(DX*DX + DY*DY + DZ*DZ);
			    
		   if ( dist < 0.500000) {
		      energy[Indx_Cell] = 9.e+99;
		      return;
		   }
	           else if (dist < 0.98000000*s_LJ_aux ) {
                       distr=dist/s_LJ_aux;
                       R126=distr*distr*distr*distr*distr*distr;
                       R1212=R126*R126;
                       Vrad=-4.0*((R1212 - R126)/(R1212*R126) + VL0_aux);
	               energ=energ + Vrad ;
                   }
		   else if (dist < range_aux) {

		      Q0j = Sphere_Cell[jmol+3];
		      Q1j = Sphere_Cell[jmol+4];
		      Q2j = Sphere_Cell[jmol+5];
		      Q3j = Sphere_Cell[jmol+6];

		      Q202 = Q0j*Q0j;
		      Q212 = Q1j*Q1j;
		      Q222 = Q2j*Q2j;
		      Q232 = Q3j*Q3j;

		      Q20Q21 = Q0j*Q1j;
		      Q20Q22 = Q0j*Q2j;
		      Q20Q23 = Q0j*Q3j;
		      Q21Q22 = Q1j*Q2j;
		      Q21Q23 = Q1j*Q3j;
		      Q22Q23 = Q2j*Q3j;

		      R211 = Q202+Q212-Q222-Q232;
		      R212 = 2.00*(Q21Q22-Q20Q23);
		      R213 = 2.00*(Q21Q23+Q20Q22);
		      R221 = 2.00*(Q21Q22+Q20Q23);
		      R222 = Q202-Q212+Q222-Q232;
		      R223 = 2.00*(Q22Q23-Q20Q21);
		      R231 = 2.00*(Q21Q23-Q20Q22);
		      R232 = 2.00*(Q22Q23+Q20Q21);
		      R233 = Q202-Q212-Q222+Q232;

		      CSM1M=-100.;
		      CSM2M=-100.;


		      for ( k=0; k < Npatches[itypei]; k++) {

		         kpatch = Npatches_max*Ndim*itypei + k*Ndim;
			 PX = patch[kpatch];
			 PY = patch[kpatch+1];
			 PZ = patch[kpatch+2];


			 PX1 = R111*PX+R112*PY+R113*PZ;
			 PY1 = R121*PX+R122*PY+R123*PZ;
			 PZ1 = R131*PX+R132*PY+R133*PZ;


			 CSOMG1=(PX1*DX+PY1*DY+PZ1*DZ)/dist;
                         COS1[k] = CSOMG1;
			 if(CSOMG1 > CSM1M) {
                               CSM1M = CSOMG1;
                         }

		         PX = ref_tor[kpatch];
		         PY = ref_tor[kpatch+1];
		         PZ = ref_tor[kpatch+2];

		         PX1_ref[k] = R111*PX+R112*PY+R113*PZ;
		         PY1_ref[k] = R121*PX+R122*PY+R123*PZ;
		         PZ1_ref[k] = R131*PX+R132*PY+R133*PZ;

		      } //for k < Npatches

		      for ( k=0; k < Npatches[itypej]; k++) {

		         kpatch = Npatches_max*Ndim*itypej + k*Ndim;
			 PX = patch[kpatch];
			 PY = patch[kpatch+1];
			 PZ = patch[kpatch+2];

			 PX2 = R211*PX+R212*PY+R213*PZ;
			 PY2 = R221*PX+R222*PY+R223*PZ;
			 PZ2 = R231*PX+R232*PY+R233*PZ;

			 CSOMG2=-(PX2*DX+PY2*DY+PZ2*DZ)/dist;
                         COS2[k] = CSOMG2;
			 if(CSOMG2 > CSM2M) {
                               CSM2M = CSOMG2;
                         }

		         PX = ref_tor[kpatch];
		         PY = ref_tor[kpatch+1];
		         PZ = ref_tor[kpatch+2];

		         PX2_ref[k] = R211*PX+R212*PY+R213*PZ;
		         PY2_ref[k] = R221*PX+R222*PY+R223*PZ;
		         PZ2_ref[k] = R231*PX+R232*PY+R233*PZ;

		      } //for k < Npatches

                
                      Vangtor_max = -1.0;
                
                      for ( j3 = 0; j3 < Npatches[itypei]; j3++) {

                         v1x=PZ1_ref[j3]*DY-PY1_ref[j3]*DZ;
                         v1y=PX1_ref[j3]*DZ-PZ1_ref[j3]*DX;
                         v1z=PY1_ref[j3]*DX-PX1_ref[j3]*DY;

                         mod_v1=sqrtf(v1x*v1x+v1y*v1y+v1z*v1z);
                         CSOMG1 = COS1[j3];

                         if( fabsf(CSOMG1) > 1.0 ) {
                            if ( fabsf(CSOMG1) - 1.0 < 0.0001) {
                                CSOMG1 = copysignf(1.0,CSOMG1); }
                            else { 
                                printf("\n error in csomg1 %f",CSOMG1);
                            } 
                         }
                         OMG1 = acosf(CSOMG1);

                         for ( j4=0; j4 <  Npatches[itypej]; j4++) {

                             Indx2 = itypei * Npatches_max + j3;  // elijo fila
                             Indx3 = Indx2 * Npatches_max * Ntypes;  // pongo el indice al inicio de la fila
                             Indx4 = Indx3 + itypej * Npatches_max;  // me muevo en la fila hasta el atomo de tipo Itypej
                             factor_energy = Vpot_Matrix[ Indx4+ j4];

                             if (factor_energy > 0.000001 ) { 

                                 CSOMG2 = COS2[j4];
                                 if( fabsf(CSOMG2) > 1.0 ) {
                                   if ( fabsf(CSOMG2) - 1.0 < 0.0001) {
                                       CSOMG2 = copysignf(1.0,CSOMG2); }
                                   else {
                                       printf("\n error in csomg2 %f",CSOMG2);
                                   } 
                                 }
                                 OMG2 = acosf(CSOMG2);
                                 IndxT1 =  Npatches_max*itypei + j3;
                                 IndxT2 =  Npatches_max*itypej + j4;

                                 sigma_jon_ij = 2.0 * sigma_jon[IndxT1] * sigma_jon[IndxT2];
                                 Vang = __expf(-(OMG1*OMG1+OMG2*OMG2)/sigma_jon_ij);

                                 if (bool_tor == 1 && Npatches[itypei] > 1 && Npatches[itypej] > 1 && Ntor_angle[IndxT1]>0 && Ntor_angle[IndxT2] > 0 ) {

                                   v2x=PZ2_ref[j4]*DY-PY2_ref[j4]*DZ;
                                   v2y=PX2_ref[j4]*DZ-PZ2_ref[j4]*DX;
                                   v2z=PY2_ref[j4]*DX-PX2_ref[j4]*DY;

                                   mod_v2=sqrtf(v2x*v2x+v2y*v2y+v2z*v2z);

                                   if( (mod_v1) != 0.0 &&  (mod_v2 != 0.0 ) ) {

                                         cangle=(v1x*v2x+v1y*v2y+v1z*v2z)/mod_v1/mod_v2;

                                         if(fabsf(cangle) > 1.0) {
                                            if( fabsf(cangle) - 1.0 < 0.0001 ) {
                                               cangle = copysignf(1.0,cangle); }
                                            else {
                                               printf("\n error in cangle %f",cangle);
                                            }     
                                         }  
                                         v12x=v1y*v2z-v1z*v2y;
                                         v12y=v1z*v2x-v1x*v2z;
                                         v12z=v1x*v2y-v1y*v2x;
                                         v12b2=(v12x*DX+v12y*DY+v12z*DZ)/dist/mod_v1/mod_v2;

                                         if( fabsf(v12b2) > 1.0) {
                                             if( fabsf(v12b2) - 1.0 < 0.0001 ) {   
                                                   v12b2 = copysignf(1.0,v12b2); }
                                             else {
                                                   printf("\n error in v12b2 %f",v12b2);
                                             }     
                                         }    

                                         torsion= atan2(v12b2,cangle);

                                         IndxT1 =  Npatches_max*itypei + j3;

                                         Utor_max = -100.0;

                                         for(k=0; k < Ntor_angle[ IndxT1 ]; k++) {

                                               dtor = torsion-tor_angle[IndxT1*Ntor_max+k] ;
                                               if (dtor > pi) {
                                                   dtor = dtor - dospi; }
                                               else if (dtor < -pi) {
                                                   dtor = dtor + dospi ; }
                                               else if (abs(dtor) > dospi) {
                                                   printf("\n ERROR");
                                               }
                                               Utor = __expf (-dtor*dtor/sigma_tor);
                                               if (Utor  > Utor_max) Utor_max = Utor;
                                         }      

                                         Vangtor = factor_energy*Vang*Utor_max;
                                   }
                                   else {
                                     //    Vangtor = Vang;
                                   }

                                 }
                                 else { /* if (bool_tor==0 ) No torsion */
                                      Vangtor = factor_energy*Vang;
                                 }

                                 if( Vangtor  > Vangtor_max)  { 
                                         Vangtor_max = Vangtor;   
                                 }

                               } /* if (factor_energy > 0.00000001 ) */
                           } /* for ( j4=0; j4 <  Npatches[itypej]; j4++) */
                        } /* for ( j3 = 0; j3 < Npatches[itypei]; j3++) */


                        if(Vangtor_max == -1.0) {    
	                   energ=energ + 0.0 ;
                        }      
                        else {

                           Vangtor=Vangtor_max;
                           FX=(1.0 + tanhf(a*(-XOP_aux + dist)))/2.0;
                           FXI=(1.0 - tanhf(a*(-XOP_aux + dist)))/2.0;
                           distr=dist/s_LJ_aux;
                           R126=distr*distr*distr*distr*distr*distr;
                           R1212=R126*R126;
                           Vrad=-4.0*((R1212 - R126)/(R1212*R126) + VL0_aux);
	                   energ=energ + Vrad * (FXI + Vangtor * FX);
                        }
      
		   }  // if dist < range
                   /* } // if Napcthes > 0 */
		   /* }  // if dist < range_long */
		}  //  for l < Nat_neigh_cell
		} // if neigh cell is not empty
	   }   // loop over neighbour cells

    }  // for m < Nat_cell
  } // if cell is not empty
  energy[Indx_Cell] = energ; // energ; 
  }  // if cell within bounds 
}  // fin kernel


///////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////


extern "C" void gpucount_(int *nogpus) {
   int ndev;
   int mpc, nthr,warps, nreg;
   size_t mem, shared_mem;
   cudaDeviceProp properties;
   // get the number of GPUs
   cudaError_t cudaRescode = cudaGetDeviceCount( &ndev ) ;
   if (cudaRescode != 0){
      *nogpus = 0;
   } else {
      *nogpus = ndev;
   }
   if(*nogpus == 1) {
   cudaGetDeviceProperties(&properties,0);
   printf("\n\n\n\n %%%%%%%%%%%%%%%%%%%%%%\n GPU properties \n ");
   mem=properties.totalGlobalMem;
   printf("\n GPU global memory %zu",mem);
   shared_mem=properties.sharedMemPerBlock;
   printf("\n Shared memory  per block %zu",shared_mem);
   mpc=properties.multiProcessorCount;
   printf("\n Number of symmetric multiprocessors %d",mpc);
   nthr=properties.maxThreadsPerBlock;
   printf("\n Maximum number of threads per block %d",nthr);
   nreg=properties.regsPerBlock;
   printf("\n Maximum number of registers per block %d",nreg);
   warps=properties.warpSize;
   printf("\n Warp size in threads %d",warps);
   printf("\n %%%%%%%%%%%%%%%%%%%%%%\n \n ");
}
}

///////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

extern "C" void choosegpu_(int *gpuid) {

   int gpuid_dev, mpc, nthr,warps, nreg;
   size_t mem, shared_mem;
   gpuid_dev = *gpuid;
   cudaDeviceProp properties;
   cudaError_t cudaRescode = cudaSetDevice( gpuid_dev ) ;
   if (cudaRescode != 0){
      printf("\n Error in choosing the GPU\n");
     exit(1);
   }

   cudaGetDeviceProperties(&properties,gpuid_dev);
   printf("\n\n\n\n %%%%%%%%%%%%%%%%%%%%%%\n GPU properties \n ");
   mem=properties.totalGlobalMem;
   printf("\n GPU global memory %zu",mem);
   shared_mem=properties.sharedMemPerBlock;
   printf("\n Shared memory  per block %zu",shared_mem);
   mpc=properties.multiProcessorCount;
   printf("\n Number of symmetric multiprocessors %d",mpc);
   nthr=properties.maxThreadsPerBlock;
   printf("\n Maximum number of threads per block %d",nthr);
   nreg=properties.regsPerBlock;
   printf("\n Maximum number of registers per block %d",nreg);
   warps=properties.warpSize;
   printf("\n Warp size in threads %d",warps);
   printf("\n %%%%%%%%%%%%%%%%%%%%%%\n \n ");
}







////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

extern "C" void computeenergygpu_(int *bool_tor, int *init, int *Ndim, int *Nat_Cell_Max, int *Ntot_Cell, int *Ntypes, 
                          int *Npatches_max, int *Ntor_max, int *Npatches, float *patch, float *ref_tor_vec, 
                          float *sigma_jon, float *sigma_tor_jon, float *sigma_LJ, float *range, float *VL0, float *XOP, 
                          int *Ntor_angle, float *tor_angle, float *Vpot_Matrix, 
                          float *h, float *Sphere_Cell, int *List_Cell, int *Map_Cell, 
			  int *Nl_Cell, int *Nat_Cell, float *DR_Cell, float *En_tot) {

  int Ndim_dev, Nat_Cell_Max_dev, Npatches_max_dev;
  int Ntot_Cell_dev, Ntor_max_dev, Ntypes_dev; 
  int i , bool_tor_dev;

  float sigma_tor_dev;


  Ndim_dev = *Ndim;
  Nat_Cell_Max_dev = *Nat_Cell_Max;
  Ntot_Cell_dev = *Ntot_Cell;
  Npatches_max_dev = *Npatches_max;
  Ntor_max_dev = *Ntor_max;
  Ntypes_dev = *Ntypes;
  bool_tor_dev = *bool_tor;

  sigma_tor_dev = *sigma_tor_jon;

/* Allocate memory on the GPU and transfer positions */

  int Ndim1= Ntypes_dev*Npatches_max_dev*Ntypes_dev*Npatches_max_dev;
  int Ndim2= Ntypes_dev*Npatches_max_dev;
  int Ndim3= Ndim2*Ntor_max_dev;


  if(*init==0) {

     energy = (float*)malloc(Ntot_Cell_dev*sizeof(float));

     MANEJA_ERROR(cudaMalloc((void**)&energydev, Ntot_Cell_dev*sizeof(float)));
     MANEJA_ERROR(cudaMalloc((void**)&DR_Cell_dev, Ntot_Cell_dev*26*3*sizeof(float)));
     MANEJA_ERROR(cudaMalloc((void**)&hdev, Ndim_dev*Ndim_dev*sizeof(float)));
     MANEJA_ERROR(cudaMalloc((void**)&Map_Cell_dev, Ntot_Cell_dev*26*sizeof(int)));
     MANEJA_ERROR(cudaMalloc((void**)&Sphere_Celldev, Ntot_Cell_dev*Nat_Cell_Max_dev*(Ndim_dev+4)*sizeof(float)));
     MANEJA_ERROR(cudaMalloc((void**)&Nat_Cell_dev, Ntot_Cell_dev*sizeof(int)));
     MANEJA_ERROR(cudaMalloc((void**)&Nl_Cell_dev, Ndim_dev*sizeof(int)));
     MANEJA_ERROR(cudaMalloc((void**)&List_Cell_dev, Ntot_Cell_dev*Nat_Cell_Max_dev*2*sizeof(int)));

     /* Todas estas son definiciones del modelo y permanecen constantes en toda la simulación */

     MANEJA_ERROR(cudaMalloc((void**)&Npatches_dev, Ntypes_dev*sizeof(int)));
     MANEJA_ERROR(cudaMalloc((void**)&patchdev, Ntypes_dev*Npatches_max_dev*Ndim_dev*sizeof(float)));
     MANEJA_ERROR(cudaMalloc((void**)&ref_tor_dev, Ntypes_dev*Npatches_max_dev*Ndim_dev*sizeof(float)));
     MANEJA_ERROR(cudaMalloc((void**)&tor_angle_dev, Ndim3*sizeof(float)));
     MANEJA_ERROR(cudaMalloc((void**)&Ntor_angle_dev, Ndim2*sizeof(int)));
     MANEJA_ERROR(cudaMalloc((void**)&Vpot_Matrix_dev, Ndim1*sizeof(float)));
     MANEJA_ERROR(cudaMalloc((void**)&sigma_jon_dev, Ntypes_dev*Npatches_max_dev*sizeof(float)));

     MANEJA_ERROR(cudaMalloc((void**)&sigma_LJ_dev, Ntypes_dev*Ntypes_dev*sizeof(float)));
     MANEJA_ERROR(cudaMalloc((void**)&XOP_dev, Ntypes_dev*Ntypes_dev*sizeof(float)));
     MANEJA_ERROR(cudaMalloc((void**)&range_dev, Ntypes_dev*Ntypes_dev*sizeof(float)));
     MANEJA_ERROR(cudaMalloc((void**)&VL0_dev, Ntypes_dev*Ntypes_dev*sizeof(float)));

     MANEJA_ERROR(cudaMemcpy( DR_Cell_dev, DR_Cell, Ntot_Cell_dev*26*3*sizeof(float), cudaMemcpyHostToDevice));  
     MANEJA_ERROR(cudaMemcpy( Map_Cell_dev, Map_Cell, Ntot_Cell_dev*26*sizeof(int), cudaMemcpyHostToDevice));  
     MANEJA_ERROR(cudaMemcpy( Nl_Cell_dev, Nl_Cell, Ndim_dev*sizeof(int), cudaMemcpyHostToDevice));  

     MANEJA_ERROR(cudaMemcpy( Npatches_dev, Npatches, Ntypes_dev*sizeof(int), cudaMemcpyHostToDevice));  
     MANEJA_ERROR(cudaMemcpy( patchdev, patch, Ntypes_dev*Npatches_max_dev*Ndim_dev*sizeof(float), cudaMemcpyHostToDevice));  
     MANEJA_ERROR(cudaMemcpy( ref_tor_dev, ref_tor_vec, Ntypes_dev*Npatches_max_dev*Ndim_dev*sizeof(float), cudaMemcpyHostToDevice));  
     MANEJA_ERROR(cudaMemcpy( tor_angle_dev, tor_angle, Ndim3*sizeof(float), cudaMemcpyHostToDevice));  
     MANEJA_ERROR(cudaMemcpy( Ntor_angle_dev, Ntor_angle, Ndim2*sizeof(int), cudaMemcpyHostToDevice));  
     MANEJA_ERROR(cudaMemcpy( Vpot_Matrix_dev, Vpot_Matrix, Ndim1*sizeof(float), cudaMemcpyHostToDevice));  
     MANEJA_ERROR(cudaMemcpy( sigma_jon_dev, sigma_jon, Ntypes_dev*Npatches_max_dev*sizeof(float), cudaMemcpyHostToDevice));  
     MANEJA_ERROR(cudaMemcpy( sigma_LJ_dev, sigma_LJ, Ntypes_dev*Ntypes_dev*sizeof(float), cudaMemcpyHostToDevice));  
     MANEJA_ERROR(cudaMemcpy( range_dev, range, Ntypes_dev*Ntypes_dev*sizeof(float), cudaMemcpyHostToDevice));  
     MANEJA_ERROR(cudaMemcpy( XOP_dev, XOP, Ntypes_dev*Ntypes_dev*sizeof(float), cudaMemcpyHostToDevice));  
     MANEJA_ERROR(cudaMemcpy( VL0_dev, VL0, Ntypes_dev*Ntypes_dev*sizeof(float), cudaMemcpyHostToDevice));  

     *init=1;
  }
  else if(*init==3) {

     MANEJA_ERROR(cudaFree(DR_Cell_dev));
     MANEJA_ERROR(cudaFree(Map_Cell_dev));
     MANEJA_ERROR(cudaFree(Sphere_Celldev));
     MANEJA_ERROR(cudaFree(Nat_Cell_dev));
     MANEJA_ERROR(cudaFree(Nl_Cell_dev));
     MANEJA_ERROR(cudaFree(energydev));
     free(energy);

     energy = (float*)malloc(Ntot_Cell_dev*sizeof(float));

     MANEJA_ERROR(cudaMalloc((void**)&energydev, Ntot_Cell_dev*sizeof(float)));
     MANEJA_ERROR(cudaMalloc((void**)&DR_Cell_dev, Ntot_Cell_dev*26*3*sizeof(float)));
     MANEJA_ERROR(cudaMalloc((void**)&Map_Cell_dev, Ntot_Cell_dev*26*sizeof(int)));
     MANEJA_ERROR(cudaMalloc((void**)&Sphere_Celldev, Ntot_Cell_dev*Nat_Cell_Max_dev*(Ndim_dev+4)*sizeof(float)));
     MANEJA_ERROR(cudaMalloc((void**)&Nat_Cell_dev, Ntot_Cell_dev*sizeof(int)));
     MANEJA_ERROR(cudaMalloc((void**)&Nl_Cell_dev, Ndim_dev*sizeof(int)));
     MANEJA_ERROR(cudaMalloc((void**)&List_Cell_dev, Ntot_Cell_dev*Nat_Cell_Max_dev*2*sizeof(int)));
     
     MANEJA_ERROR(cudaMemcpy( DR_Cell_dev, DR_Cell, Ntot_Cell_dev*26*3*sizeof(float), cudaMemcpyHostToDevice));  
     MANEJA_ERROR(cudaMemcpy( Map_Cell_dev, Map_Cell, Ntot_Cell_dev*26*sizeof(int), cudaMemcpyHostToDevice));  
     MANEJA_ERROR(cudaMemcpy( Nl_Cell_dev, Nl_Cell, Ndim_dev*sizeof(int), cudaMemcpyHostToDevice));  

     *init=1;
   }

  MANEJA_ERROR(cudaMemcpy( List_Cell_dev, List_Cell, Ntot_Cell_dev*Nat_Cell_Max_dev*2*sizeof(int), cudaMemcpyHostToDevice));  
  MANEJA_ERROR(cudaMemcpy( Sphere_Celldev, Sphere_Cell, Ntot_Cell_dev*Nat_Cell_Max_dev*(Ndim_dev+4)*sizeof(float), cudaMemcpyHostToDevice));
  MANEJA_ERROR(cudaMemcpy( Nat_Cell_dev, Nat_Cell, Ntot_Cell_dev*sizeof(int), cudaMemcpyHostToDevice));  
  MANEJA_ERROR(cudaMemcpy( hdev, h, Ndim_dev*Ndim_dev*sizeof(float), cudaMemcpyHostToDevice));  

  MANEJA_ERROR(cudaMemset(energydev, 0, Ntot_Cell_dev*sizeof(float)));

  int nblock = Ntot_Cell_dev/NTHREAD;
  if( nblock== 0 ) { nblock = 1 ; }

  // Call GPU kernel for energy calculation

  int Ntot_jobs = nblock*NTHREAD;

  if( Ntot_jobs < Ntot_Cell_dev ) {
       nblock = nblock+1;
  }

  enerGPU<<<nblock,NTHREAD>>>( bool_tor_dev, Ndim_dev, Nat_Cell_Max_dev, Ntot_Cell_dev, Npatches_max_dev, Ntor_max_dev,
                                 Ntypes_dev, Npatches_dev, Ntor_angle_dev, tor_angle_dev,
                                 patchdev, ref_tor_dev, range_dev, sigma_jon_dev, sigma_tor_dev, sigma_LJ_dev, 
                                 VL0_dev, XOP_dev,
                                 Vpot_Matrix_dev, hdev, Sphere_Celldev, List_Cell_dev,
                                 Map_Cell_dev, Nl_Cell_dev, Nat_Cell_dev, DR_Cell_dev, energydev);
  //printf("\n after kernel");

  // Recall energy per particle from the GPU
  MANEJA_ERROR(cudaMemcpy( energy, energydev, Ntot_Cell_dev*sizeof(float), cudaMemcpyDeviceToHost ) );
  MANEJA_ERROR(cudaMemcpy( Sphere_Cell, Sphere_Celldev, Ntot_Cell_dev*Nat_Cell_Max_dev*(Ndim_dev+4)*sizeof(float), cudaMemcpyDeviceToHost));


  //   Compute total energy
   float etot = 0 ;
   for (i=0; i < Ntot_Cell_dev; i++) {
      etot = etot +energy[i];
   }
   etot = etot/2.;
   *En_tot = etot;

  if(*init==2) {

     MANEJA_ERROR(cudaFree(Map_Cell_dev));
     MANEJA_ERROR(cudaFree(Offset_dev));
     MANEJA_ERROR(cudaFree(Nl_Cell_dev));
     MANEJA_ERROR(cudaFree(Nat_Cell_dev));
     MANEJA_ERROR(cudaFree(Npatches_dev));
     MANEJA_ERROR(cudaFree(Ntor_angle_dev));
     MANEJA_ERROR(cudaFree(List_Cell_dev));

     MANEJA_ERROR(cudaFree(W_RU_dev));
     MANEJA_ERROR(cudaFree(patchdev));
     MANEJA_ERROR(cudaFree(hdev));
     MANEJA_ERROR(cudaFree(DR_Cell_dev));
     MANEJA_ERROR(cudaFree(Sphere_Celldev));
     MANEJA_ERROR(cudaFree(ref_tor_dev));
     MANEJA_ERROR(cudaFree(tor_angle_dev));
     MANEJA_ERROR(cudaFree(Vpot_Matrix_dev));

  }
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

__global__ void subsweep(curandState_t *RND_state, int *overlap_error, int bool_tor, int Ndim, int Nmove, int Nat_Cell_Max, 
		 int Ntot_Cell, int Nsubset_Cell, int Npatches_max, int Ntor_max, int Ntypes, 
                 int *Npatches, float *patch, float *ref_tor, float *sigma_jon, 
                 float sigma_tor, float *sigma_LJ, float *rangeP, float *VL0, float *XOP, 
                 int *Ntor_angle, float *tor_angle, float *Vpot_Matrix, 
		 float *h, float beta, float hmax, float omax, float *Sphere_Cell, int *List_Cell, int *Map_Cell, 
		 int *Nl_Cell, int *Nat_Cell, float *DR_Cell, int *Offset_CB, double *W_RU, 
		 float *Delta_energy, int *Ntras_ac, int *Nrot_ac)
{

  int Nl0, Nl1, Nl2, Nl02, Nl12, Nl22;
  __shared__ int Indx_Cell, Nat_Cell_ijk, itypei, Npatchesi;
  int Indx_Rnd,indx_atm_aux, indx_atm_aux2;
  int j, k, l, s, overlap, kpatch, imol; 
  int ICELL, indx_test, jmol, ntrasa, nrota;
  int Indx_Cell2, Indx_Cell3, Nat_Cell_neigh,  indx_neigh_aux2, itras;
  int indx_neigh_auxp, imol2, jmol2, itypej;
  int Indx2, Indx3, Indx4, Indxk;
  int j3, j4, IndxT1, IndxT2;
  int i, indx_thread, indx_rnd_x, indx_rnd_y, indx_rnd_z;
  int Npatchesj, kpatch_aux, jcell_x, jcell_y, jcell_z;

  float DX, DY, DZ;
  float DXn, DYn, DZn, distn, dist;
  float H11, H12, H13, H21, H22, H23, H31, H32, H33;
  float Q0j, Q1j, Q2j, Q3j;
  float PX, PY, PZ, PX1, PY1, PZ1, PX2, PY2, PZ2;
  float dx_neigh, dy_neigh, dz_neigh, PX1n, PY1n, PZ1n;
  float norm, factor_energy;
  float XJ, YJ, ZJ, prob;
  float PX1_ref, PX2_ref[20];
  float PY1_ref, PY2_ref[20];
  float PZ1_ref, PZ2_ref[20];
  float PX1N_ref, PY1N_ref, PZ1N_ref;
  float Rnd_or[4];
  float COS2[20], R2[9];
  float COS2N[20];
  float Vangtor_max, v1x, v1y, v1z, mod_v1, FX, FXI;
  float Vangtor, Vang, Vrad, v2x, v2y, v2z, mod_v2, torsion ;
  float Utor_max, Utor;
  float Vangtor_n, VangN, v2xn, v2yn, v2zn, mod_v2n, torsionn;
  float Utor_max_n, sigma_jon_ij;
  float Vangtor_max_n, v1xn, v1yn, v1zn, mod_v1n;
  float dtor, VradN, sigma_LJ_aux, XOP_aux, VL0_aux, rangeP_aux;
  float a =1000.0;
  float pi= acosf(-1.0), dospi=2.*pi;
  float DE, ener_new, ener_old, de1;
  float OMG1, OMG2, OMG1N, OMG2N, CSOMG1N, CSOMG1;
  float Q0_old, Q1_old, Q2_old, Q3_old, Qmove_0, Qmove_1, Qmove_2, Qmove_3;

  __shared__ float R1[9], R1N[9] ;
  __shared__ float Dx_old, Dy_old, Dz_old, Dmove_x, Dmove_y, Dmove_z;
  __shared__ float entot_new[27], entot_old[27];
  __shared__ int overlap_aux[27];

   // ahora 27 hilos ( o 32 para hacer que sea medio warp) se ocupan de una celda 
   // primero tengo que elegir la celda sobre la que trabajo
   // puedo considerar que cada bloque se ocupa de una celda

   int indx_rnd = blockIdx.x ;  // cada bloque se ocupa de una celda
   Nl0 = Nl_Cell[0];
   Nl1 = Nl_Cell[1];
   Nl2 = Nl_Cell[2];
   Nl02 = Nl0/2;
   Nl12 = Nl1/2;
   Nl22 = Nl2/2;


   indx_rnd_z = indx_rnd/(Nl02*Nl12)  ;
   indx_rnd_y = ((indx_rnd - Nl02*Nl12*indx_rnd_z)/Nl02)  ;
   if (indx_rnd_y == 0 && indx_rnd_z == 0) {
      indx_rnd_x = indx_rnd ;
   }
   else {
      indx_rnd_x = ( indx_rnd % (indx_rnd_y*Nl02+indx_rnd_z*Nl02*Nl12)) ;
   }

   jcell_x = 2*indx_rnd_x+Offset_CB[0];
   jcell_y = 2*indx_rnd_y+Offset_CB[1];
   jcell_z = 2*indx_rnd_z+Offset_CB[2];    // celda elegida: (i,j,k)
   indx_thread = threadIdx.x;

   if ( jcell_x < Nl0  && jcell_y < Nl1 && jcell_z < Nl2 ) {


         Indx_Cell=((jcell_x + Nl0)%Nl0) + ((jcell_y+Nl1)%Nl1)*Nl0 + ((jcell_z+Nl2)%Nl2)* Nl0*Nl1;
         Indx_Rnd = ((indx_rnd_x + Nl02)%Nl02) + ((indx_rnd_y+Nl12)%Nl12)*Nl02 + ((indx_rnd_z+Nl22)%Nl22)* Nl02*Nl12;
         Nat_Cell_ijk = Nat_Cell[Indx_Cell];
         DE = 0.0;
         ntrasa = 0;
         nrota = 0;
     __syncthreads();

     // set logical variable overlap equal to zero, using an integer, 0-> no overlap, 1-> overlap

     if ( Nat_Cell_ijk > 0 ) {

         H11=h[0];
         H12=h[1];
         H13=h[2];
         H21=h[3];
         H22=h[4];
         H23=h[5];
         H31=h[6];
         H32=h[7];
         H33=h[8];

         indx_atm_aux = Indx_Cell * Nat_Cell_Max * 7;
         indx_atm_aux2 = Indx_Cell * Nat_Cell_Max * 2;

         for (s=0; s< Nmove; s++) {   // do Nmove trials

///////////////////////  only done in the master thread   //////////////

           if(indx_thread == 0 ) { 


               if( Indx_Rnd >= Nsubset_Cell ) { printf("\n error in rnd 1");}
               indx_test= floorf(Nat_Cell_ijk*curand_uniform(&RND_state[Indx_Rnd]));  
               if(indx_test >= Nat_Cell_ijk) {
                     indx_test = Nat_Cell_ijk - 1;
               }
               if(indx_test < 0) {indx_test=0;}

               imol = indx_atm_aux+ indx_test*7;
               imol2 = indx_atm_aux2+ indx_test*2;

               itypei = List_Cell [ imol2 +1];
               Npatchesi = Npatches[itypei];

               Dx_old = Sphere_Cell[imol];
               Dy_old = Sphere_Cell[imol+1];
               Dz_old = Sphere_Cell[imol+2];
               Q0_old = Sphere_Cell[imol+3];
               Q1_old = Sphere_Cell[imol+4];
               Q2_old = Sphere_Cell[imol+5];
               Q3_old = Sphere_Cell[imol+6]; 

               if( curand_uniform(&RND_state[Indx_Rnd]) < 0.5 ) {

                   itras = 1;
                   Dmove_x = Dx_old + (curand_uniform(&RND_state[Indx_Rnd])-0.5)*hmax;
                   Dmove_y = Dy_old + (curand_uniform(&RND_state[Indx_Rnd])-0.5)*hmax;
                   Dmove_z = Dz_old + (curand_uniform(&RND_state[Indx_Rnd])-0.5)*hmax;

                   Qmove_0 = Q0_old;
                   Qmove_1 = Q1_old;
                   Qmove_2 = Q2_old;
                   Qmove_3 = Q3_old;

               } // end of traslational move
               else {

                   itras = 0;
                   Dmove_x = Dx_old;
                   Dmove_y = Dy_old;
                   Dmove_z = Dz_old;

                   genranor( RND_state, Indx_Rnd, Rnd_or);

                   Qmove_0 = Q0_old + Rnd_or[0]*omax ;
                   Qmove_1 = Q1_old + Rnd_or[1]*omax ;
                   Qmove_2 = Q2_old + Rnd_or[2]*omax ;
                   Qmove_3 = Q3_old + Rnd_or[3]*omax ;
                   norm=sqrtf(Qmove_0*Qmove_0+Qmove_1*Qmove_1+Qmove_2*Qmove_2+Qmove_3*Qmove_3);
                   Qmove_0 = Qmove_0 / norm;
                   Qmove_1 = Qmove_1 / norm;
                   Qmove_2 = Qmove_2 / norm;
                   Qmove_3 = Qmove_3 / norm;

               } // end of rotational move

               rot_matrix(Q0_old,Q1_old,Q2_old,Q3_old,R1);   // Dmove, DX_old, R1 and R1N need to be shared
               rot_matrix(Qmove_0,Qmove_1,Qmove_2,Qmove_3,R1N); 

          }    /////   end of tasks of master thread



          entot_old[indx_thread] = 0.0;
          entot_new[indx_thread] = 0.0;//}
          overlap_aux[indx_thread] = 0;

          __syncthreads();

          if ( fabsf(Dmove_x)  >= 0.5* W_RU[0] ) {
          }
          else if ( fabsf( Dmove_y ) >= 0.5* W_RU[1] ) {
          }
          else if ( fabsf( Dmove_z ) >= 0.5* W_RU[2] ) {
          }
          else {

            
           // calculate the energy of the old and new configuration


////////////////////// master thread /////////////////////////////////////////

          if( indx_thread == 0 ) {

              ener_old =0.0;
              ener_new =0.0;
              overlap = 0;

              for (j =0; j < Nat_Cell_ijk; j++) {

	         if ( j != indx_test ) {  

	             jmol= indx_atm_aux + j*7;
                     jmol2 =  indx_atm_aux2+j*2;

	             XJ = Sphere_Cell[jmol];
	             YJ = Sphere_Cell[jmol+1];
	             ZJ = Sphere_Cell[jmol+2];


                     itypej = List_Cell[jmol2+1 ];

                     Indxk = Ntypes * itypei + itypej;
                     sigma_LJ_aux = sigma_LJ[Indxk];
                     XOP_aux = XOP[Indxk];
                     VL0_aux = VL0[Indxk];
                     rangeP_aux = rangeP[Indxk];

                     get_dist(XJ,YJ,ZJ,Dx_old,Dy_old,Dz_old,H11,H12,H13,H21,H22,H23,H31,H32,H33,DX,DY,DZ,dist);
                     get_dist(XJ,YJ,ZJ,Dmove_x,Dmove_y,Dmove_z,H11,H12,H13,H21,H22,H23,H31,H32,H33,DXn,DYn,DZn,distn);

	             if ( dist < 0.4999998) {
		         printf("\n ERROR 1: overlap in old configuration, i: %d %f %f %f, j: %d %f %f %f, dist: %f %f %f: %f",indx_test,Dx_old,Dy_old,Dz_old,j,XJ,YJ,ZJ,DX,DY,DZ,dist);
		         overlap_error[Indx_Cell] =1;
                         printf("\n Dmove_x %f 0.5* W_RU[0] %f",Dmove_x, 0.5* W_RU[0]);
                         printf("\n Dmove_y %f 0.5* W_RU[1] %f",Dmove_y, 0.5* W_RU[1]);
                         printf("\n Dmove_z %f 0.5* W_RU[2] %f",Dmove_z, 0.5* W_RU[2]);
		         return;
	             }
	             if ( distn < 0.5000000) {
		         overlap=1;
                         ener_new=ener_new+99999999999.;
	             }
                     else if (dist < 0.98000000*sigma_LJ_aux && distn < 0.98000000*sigma_LJ_aux) {
                          vlj(dist,sigma_LJ_aux,VL0_aux,Vrad);
                          ener_old=ener_old + Vrad ;

                          vlj(distn,sigma_LJ_aux,VL0_aux,Vrad);
                          ener_new=ener_new + Vrad ;
                     }
                     else if (dist < rangeP_aux || distn < rangeP_aux) {

                        Q0j = Sphere_Cell[jmol+3];
                        Q1j = Sphere_Cell[jmol+4];
                        Q2j = Sphere_Cell[jmol+5];
                        Q3j = Sphere_Cell[jmol+6];
                        rot_matrix(Q0j,Q1j,Q2j,Q3j,R2);
                        Npatchesj = Npatches[itypej];
                        kpatch_aux= Npatches_max*Ndim*itypej;

                        for ( k=0; k < Npatchesj; k++) {

                           kpatch = kpatch_aux + k*Ndim;
                           PX = patch[kpatch];
                           PY = patch[kpatch+1];
                           PZ = patch[kpatch+2];
                           patch_real (PX, PY, PZ, R2, PX2, PY2, PZ2 );
                           COS2[k]=-(PX2*DX+PY2*DY+PZ2*DZ)/dist;
                           COS2N[k]=-(PX2*DXn+PY2*DYn+PZ2*DZn)/distn;

                           PX = ref_tor[kpatch];
                           PY = ref_tor[kpatch+1];
                           PZ = ref_tor[kpatch+2];
                           patch_real (PX, PY, PZ, R2, PX2_ref[k], PY2_ref[k],PZ2_ref[k] );

                        } //for k < Npatches

                        Vangtor_max = -1.0;
                        Vangtor_max_n = -1.0;

                        kpatch_aux=Npatches_max*Ndim*itypei;

                        for ( j3 = 0; j3 < Npatchesi; j3++) {

                           kpatch = kpatch_aux + j3*Ndim;
                           PX = patch[kpatch];
                           PY = patch[kpatch+1];
                           PZ = patch[kpatch+2];
                           patch_real (PX, PY, PZ, R1, PX2, PY2, PZ2 );
                           CSOMG1=(PX2*DX+PY2*DY+PZ2*DZ)/dist;
                           patch_real (PX, PY, PZ, R1N, PX2, PY2, PZ2 );
                           CSOMG1N=(PX2*DXn+PY2*DYn+PZ2*DZn)/distn;

                           PX = ref_tor[kpatch];
                           PY = ref_tor[kpatch+1];
                           PZ = ref_tor[kpatch+2];
                           patch_real (PX, PY, PZ, R1, PX1_ref, PY1_ref, PZ1_ref );
                           patch_real (PX, PY, PZ, R1N, PX1N_ref, PY1N_ref, PZ1N_ref );

                           // old configuration //
                           get_torsion(PX1_ref, PY1_ref, PZ1_ref, DX, DY, DZ, v1x, v1y, v1z, mod_v1);
                           get_cosine(CSOMG1, OMG1);

                           // new configuration //
                          get_torsion(PX1N_ref, PY1N_ref, PZ1N_ref, DXn, DYn, DZn, v1xn, v1yn, v1zn, mod_v1n);
                           get_cosine(CSOMG1N, OMG1N);

                          Indx2 = itypei * Npatches_max + j3 ;  // elijo fila
                          Indx3 = Indx2 * Npatches_max * Ntypes;  // pongo el indice al inicio de la fila
                          Indx4 = Indx3 + itypej * Npatches_max;  // me muevo en la fila hasta el atomo de tipo Itypej
                          IndxT1 =  Npatches_max*itypei + j3;

                          for ( j4=0; j4 <  Npatchesj; j4++) {

                              factor_energy = Vpot_Matrix[ Indx4+ j4 ];

                              if (factor_energy > 0.000001 ) {
                              // old configuration //
                              get_cosine(COS2[j4], OMG2);

                              IndxT2 =  Npatches_max*itypej + j4;
                              sigma_jon_ij = 2.0* sigma_jon[IndxT1]*sigma_jon[IndxT2];

                              Vang = __expf(-(OMG1*OMG1+OMG2*OMG2)/sigma_jon_ij);

                              // new configuration //
                              get_cosine(COS2N[j4], OMG2N);
                              VangN = __expf(-(OMG1N*OMG1N+OMG2N*OMG2N)/sigma_jon_ij);

                              if (bool_tor == 1 && Npatchesi > 1 && Npatchesj > 1 && Ntor_angle[IndxT1]> 0 && Ntor_angle[IndxT2]> 0 ) {

                                 get_torsion(PX2_ref[j4], PY2_ref[j4], PZ2_ref[j4],
                                          DX, DY, DZ, v2x, v2y, v2z, mod_v2);
                                 get_torsion(PX2_ref[j4], PY2_ref[j4], PZ2_ref[j4],
                                          DXn, DYn, DZn, v2xn, v2yn, v2zn, mod_v2n);

                                 if( mod_v1 != 0.0 &&  mod_v2 != 0.0 ) {

                                    // old configuration //
                                    get_torsion_angle(v1x, v1y, v1z, mod_v1, v2x, v2y, v2z,
                                           mod_v2, DX, DY, DZ, dist, torsion);

                                    Utor_max = -100.0;
                                    for(k=0; k < Ntor_angle[ IndxT1 ]; k++) {
                                          dtor = torsion-tor_angle[IndxT1*Ntor_max+k] ;
                                          if (dtor > pi) {
                                              dtor = dtor - dospi; }
                                          else if (dtor < -pi) {
                                              dtor = dtor + dospi ;  }
                                          else if (abs(dtor) > dospi) {
                                                   printf("\n ERROR");
                                          }
                                          Utor = __expf (-dtor*dtor/sigma_tor);
                                          if (Utor  > Utor_max) Utor_max = Utor;
                                     }
                                     Vangtor = factor_energy*Vang*Utor_max;
                                 }
                                 else {
                                    //Vangtor = Vang;
                                 }
                                 if(   mod_v1n != 0.0 &&  mod_v2n != 0.0  ) {
                                    // new configuration //
                                    get_torsion_angle(v1xn, v1yn, v1zn, mod_v1n, v2xn, v2yn, v2zn,
                                           mod_v2n, DXn, DYn, DZn, distn, torsionn);

                                    IndxT1 =  Npatches_max*itypei + j3;
                                    Utor_max_n = -100.0;
                                    for(k=0; k < Ntor_angle[ IndxT1 ]; k++) {
                                          dtor = torsionn-tor_angle[IndxT1*Ntor_max+k] ;
                                          if (dtor > pi) {
                                              dtor = dtor - dospi; }
                                          else if (dtor < -pi) {
                                              dtor = dtor + dospi ;  }
                                          else if (abs(dtor) > dospi) {
                                                   printf("\n ERROR");
                                          }
                                          Utor = __expf (-dtor*dtor/sigma_tor);
                                          if (Utor  > Utor_max_n) Utor_max_n = Utor;
                                    }
                                    Vangtor_n = factor_energy*VangN*Utor_max_n;
                                }
                                else {
                                    //Vangtor_n = VangN;
                                }
                              }
                              else { // if (bool_tor==1 ) //
                               Vangtor = factor_energy*Vang;
                               Vangtor_n = factor_energy*VangN;
                              }  // if (bool_tor) //

                              if( Vangtor  > Vangtor_max)  Vangtor_max = Vangtor;
                              if( Vangtor_n  > Vangtor_max_n)  Vangtor_max_n = Vangtor_n;
                             } // if (factor_energy > 0.00000001 ) //
                          } // for ( j4=0; j4 <  Npatches[itypej]; j4++) //
                       } // for ( j3 = 0; j3 < Npatches[itypei]; j3++) //

                       // old configuration //
                       if (dist < 0.98000000*sigma_LJ_aux ) {
                            vlj(dist, sigma_LJ_aux, VL0_aux, Vrad);
                            ener_old=ener_old + Vrad;
                       }
                       else if (dist < rangeP_aux ) {
                          if(Vangtor_max == -1.0) {
                             ener_old=ener_old + 0.0  ;
                          }
                          else {
                             Vangtor=Vangtor_max;
                             FX=(1.0 + tanh(a*(-XOP_aux + dist)))/2.0;
                             FXI=(1.0 - tanh(a*(-XOP_aux + dist)))/2.0;
                             vlj(dist, sigma_LJ_aux, VL0_aux, Vrad);
                             ener_old=ener_old + Vrad * (FXI + Vangtor * FX);
                          }
                       }
                       // new configuration //
                       if (distn < 0.98000000*sigma_LJ_aux ) {
                             vlj(distn, sigma_LJ_aux, VL0_aux, Vrad);
                             ener_new=ener_new + Vrad;
                       }
                       else if (distn < rangeP_aux ) {
                          if(Vangtor_max_n == -1.0) {
                             ener_new=ener_new + 0.0 ;
                          }
                          else {
                             Vangtor=Vangtor_max_n;
                             FX=(1.0 + tanh(a*(-XOP_aux + distn)))/2.0;
                             FXI=(1.0 - tanh(a*(-XOP_aux + distn)))/2.0;
                             vlj(distn, sigma_LJ_aux, VL0_aux, VradN);
                             ener_new=ener_new + VradN * (FXI + Vangtor * FX);
                          }
                       }


	        } // if dist < range_KF
             } // if j =! itest
          } // for j < Nat_Cell_ij
          //if(indx_thread > 31) { printf("\n error indx_thread larger than 31");}
          entot_old[indx_thread] = ener_old;
          entot_new[indx_thread] = ener_new;
          overlap_aux[indx_thread] = overlap;

///////////////////////  only done in the master thread   //////////////
           }
////////////////////// remaining threads /////////////////////////////////////////
           else { 

                 ICELL = indx_thread-1;

                 ener_old=0.0;
                 ener_new=0.0;
                 overlap= 0;

                 if(ICELL < 26 ) {

	             Indx_Cell2 = Indx_Cell*26 + ICELL ;
	             Indx_Cell2 = Map_Cell [Indx_Cell2];

	             Nat_Cell_neigh = Nat_Cell[Indx_Cell2];

	             Indx_Cell3 = Indx_Cell*26*3 + ICELL*3 ;
	             dx_neigh = DR_Cell[Indx_Cell3];
	             dy_neigh = DR_Cell[Indx_Cell3+1];
	             dz_neigh = DR_Cell[Indx_Cell3+2];

	             if ( Nat_Cell_neigh > 0 ) {

	                 indx_neigh_aux2 = Indx_Cell2 * Nat_Cell_Max * 7;
	                 indx_neigh_auxp = Indx_Cell2 * Nat_Cell_Max * 2;

	                 for (l=0; l< Nat_Cell_neigh; l++) {

		               jmol = indx_neigh_aux2 + l*7;
		               jmol2 = indx_neigh_auxp + l*2;

		               XJ = Sphere_Cell[jmol] + dx_neigh;
		               YJ = Sphere_Cell[jmol+1] + dy_neigh;
		               ZJ = Sphere_Cell[jmol+2] + dz_neigh ;

                               itypej = List_Cell [ jmol2 +1 ];

                               Indxk = Ntypes * itypei + itypej;
                               sigma_LJ_aux = sigma_LJ[Indxk];
                               XOP_aux = XOP[Indxk];
                               VL0_aux = VL0[Indxk];
                               rangeP_aux = rangeP[Indxk];

                               get_dist(XJ,YJ,ZJ,Dx_old,Dy_old,Dz_old,H11,H12,H13,H21,H22,H23,H31,H32,H33,DX,DY,DZ,dist);
                               get_dist(XJ,YJ,ZJ,Dmove_x,Dmove_y,Dmove_z,H11,H12,H13,H21,H22,H23,H31,H32,H33,DXn,DYn,DZn,distn);

		               if ( dist < 0.4999998) {
		                  printf("\n ERROR 2: overlap in old configuration");
		                  printf("\n DX_old %f %f %f ",Dx_old,Dy_old,Dz_old);
		                  printf("\n Dmove_x %f %f %f ",Dmove_x,Dmove_y,Dmove_z);
		                  printf("\n dist distn %f %f ",dist,distn);
                                  printf("\n dx_neigh %f %f %f",dx_neigh,dy_neigh,dz_neigh);
                                  printf("\n s %d",s);
                                  printf("\n itras %d",itras);
                                  printf("\n jcell %d %d",ICELL,Indx_Cell2);
                                  printf("\n icell %d ",Indx_Cell);

                         printf("\n Dmove_x %f 0.5* W_RU[0] %f",Dmove_x, 0.5* W_RU[0]);
                         printf("\n Dmove_y %f 0.5* W_RU[1] %f",Dmove_y, 0.5* W_RU[1]);
                         printf("\n Dmove_z %f 0.5* W_RU[2] %f",Dmove_z, 0.5* W_RU[2]);
		                  overlap_error[Indx_Cell]=1;
		                  return; 
		               } 
		               if ( distn < 0.5000000) {
		                  overlap=1;
                                  ener_new=ener_new+99999999999.;
		               }
		               else if (dist < 0.98000000*sigma_LJ_aux && distn < 0.98000000*sigma_LJ_aux)  {
                                   vlj(dist,sigma_LJ_aux,VL0_aux,Vrad);
                                   ener_old=ener_old + Vrad ;

                                   vlj(distn,sigma_LJ_aux,VL0_aux,Vrad);
                                   ener_new=ener_new + Vrad ;
                               }
                               else if (dist < rangeP_aux || distn < rangeP_aux) {

                                    Q0j = Sphere_Cell[jmol+3];
                                    Q1j = Sphere_Cell[jmol+4];
                                    Q2j = Sphere_Cell[jmol+5];
                                    Q3j = Sphere_Cell[jmol+6];
                                    rot_matrix(Q0j,Q1j,Q2j,Q3j,R2);

                                    Npatchesj= Npatches[itypej];
                                    kpatch_aux = Npatches_max*Ndim*itypej;

                                    for ( k=0; k < Npatchesj; k++) {

                                       kpatch = kpatch_aux + k*Ndim;
                                       PX = patch[kpatch];
                                       PY = patch[kpatch+1];
                                       PZ = patch[kpatch+2];
                                       patch_real (PX, PY, PZ, R2, PX2, PY2, PZ2 );
                                       COS2[k]= -(PX2*DX+PY2*DY+PZ2*DZ)/dist;
                                       COS2N[k]= -(PX2*DXn+PY2*DYn+PZ2*DZn)/distn;

                                       PX = ref_tor[kpatch];
                                       PY = ref_tor[kpatch+1];
                                       PZ = ref_tor[kpatch+2];
                                       patch_real (PX, PY, PZ, R2, PX2_ref[k], PY2_ref[k], PZ2_ref[k] ); 

                                    } //for k < Npatches

                                    Vangtor_max = -1.0;
                                    Vangtor_max_n = -1.0;

                                    kpatch_aux = Npatches_max*Ndim*itypei;

                                    for ( j3 = 0; j3 < Npatchesi; j3++) {

                                       kpatch = kpatch_aux + j3*Ndim;
                                       PX = patch[kpatch];
                                       PY = patch[kpatch+1];
                                       PZ = patch[kpatch+2];
                                       patch_real (PX, PY, PZ, R1, PX1, PY1, PZ1 );
                                       patch_real (PX, PY, PZ, R1N, PX1n, PY1n, PZ1n );

                                       CSOMG1=(PX1*DX+PY1*DY+PZ1*DZ)/dist;
                                       CSOMG1N=(PX1n*DXn+PY1n*DYn+PZ1n*DZn)/distn;

                                       PX = ref_tor[kpatch];
                                       PY = ref_tor[kpatch+1];
                                       PZ = ref_tor[kpatch+2];
                                       patch_real (PX, PY, PZ, R1, PX1_ref, PY1_ref, PZ1_ref );
                                       patch_real (PX, PY, PZ, R1N, PX1N_ref, PY1N_ref, PZ1N_ref );


                                      // old configuration //
                                       get_torsion(PX1_ref, PY1_ref, PZ1_ref, DX, DY, DZ, v1x, v1y, v1z, mod_v1);
                                       get_cosine(CSOMG1,OMG1);

                                       // new configuration //
                                       get_torsion(PX1N_ref, PY1N_ref, PZ1N_ref, DXn, DYn, DZn, v1xn, v1yn, v1zn, mod_v1n);
                                       get_cosine(CSOMG1N,OMG1N);

                                       Indx2 = itypei * Npatches_max + j3 ;  // elijo fila
                                       Indx3 = Indx2 * Npatches_max * Ntypes;  // pongo el indice al inicio de la fila
                                       Indx4 = Indx3 + itypej * Npatches_max;  // me muevo en la fila hasta el atomo de tipo Itypej
                                       IndxT1 =  Npatches_max*itypei + j3;

                                       for ( j4=0; j4 <  Npatchesj; j4++) {

                                           factor_energy = Vpot_Matrix[ Indx4+ j4 ];

                                           if (factor_energy > 0.000001 ) {

                                               // old configuration //
                                               get_cosine(COS2[j4],OMG2);
                                               IndxT2 =  Npatches_max*itypej + j4;
                                               sigma_jon_ij = 2.0* sigma_jon[IndxT1]*sigma_jon[IndxT2];
                                               Vang = __expf(-(OMG1*OMG1+OMG2*OMG2)/sigma_jon_ij);

                                               // new configuration //
                                               get_cosine(COS2N[j4],OMG2N);
                                               VangN = __expf(-(OMG1N*OMG1N+OMG2N*OMG2N)/sigma_jon_ij);

                                              if (bool_tor == 1 && Npatchesi > 1 && Npatchesj > 1 && Ntor_angle[IndxT1] > 0 && Ntor_angle[IndxT2] > 0 ) {

                                                 get_torsion(PX2_ref[j4], PY2_ref[j4], PZ2_ref[j4], DX, DY, DZ,
                                                          v2x, v2y, v2z, mod_v2);

                                                 get_torsion(PX2_ref[j4], PY2_ref[j4], PZ2_ref[j4], DXn, DYn, DZn,
                                                          v2xn, v2yn, v2zn, mod_v2n);

                                                 if( mod_v1 != 0.0 &&  mod_v2 != 0.0 ){
                                                     // old configuration //
                                                     get_torsion_angle(v1x, v1y, v1z, mod_v1, v2x, v2y, v2z,
                                                            mod_v2, DX, DY, DZ, dist, torsion);
                                                     Utor_max = -100.0;
                                                     for(k=0; k < Ntor_angle[ IndxT1 ]; k++) {
                                                           dtor = torsion-tor_angle[IndxT1*Ntor_max+k] ;
                                                           if (dtor > pi) {
                                                               dtor = dtor - dospi; }
                                                           else if (dtor < -pi) {
                                                               dtor = dtor + dospi ;  }
                                                           else if (abs(dtor) > dospi) {
                                                               printf("\n ERROR");
                                                           }

                                                           Utor = __expf (-dtor*dtor/sigma_tor);
                                                           if (Utor  > Utor_max) Utor_max = Utor;
                                                     }
                                                     Vangtor = factor_energy*Vang*Utor_max;
                                                 }
                                                 else {
                                                      //Vangtor = Vang;
                                                 }
                                                 ///////////////
                                                 if(  mod_v1n != 0.0 &&  mod_v2n != 0.0  ) {
                                                     // new configuration //
                                                     get_torsion_angle(v1xn, v1yn, v1zn, mod_v1n, v2xn, v2yn, v2zn,
                                                        mod_v2n, DXn, DYn, DZn, distn, torsionn);
                                                     Utor_max_n = -100.0;
                                                     for(k=0; k < Ntor_angle[ IndxT1 ]; k++) {
                                                           dtor = torsionn-tor_angle[IndxT1*Ntor_max+k] ;
                                                           if (dtor > pi) {
                                                               dtor = dtor - dospi; }
                                                           else if (dtor < -pi) {
                                                               dtor = dtor + dospi ;  }
                                                           else if (abs(dtor) > dospi) {
                                                               printf("\n ERROR");
                                                           }
                                                           Utor = __expf (-dtor*dtor/sigma_tor);
                                                           if (Utor  > Utor_max_n) Utor_max_n = Utor;
                                                     }
                                                     Vangtor_n = factor_energy*VangN*Utor_max_n;
                                                 }
                                                 else {
                                                     //Vangtor_n = VangN;
                                                 }
                                                 ///////////////
                                              }
                                             else { // if (bool_tor ==0 ) No torsion //

                                                    Vangtor = factor_energy*Vang;
                                                    Vangtor_n = factor_energy*VangN;

                                              }

                                              if( Vangtor  > Vangtor_max)  Vangtor_max = Vangtor;
                                              if( Vangtor_n  > Vangtor_max_n)  Vangtor_max_n = Vangtor_n;

                                          } // if (factor_energy > 0.00000001 ) //
                                      } // for ( j4=0; j4 <  Npatches[itypej]; j4++) //
                                   } // for ( j3 = 0; j3 < Npatches[itypei]; j3++) //

                                   // old configuration //
                                   if( dist < 0.98000000*sigma_LJ_aux) {
                                       vlj(dist, sigma_LJ_aux, VL0_aux, Vrad);
                                       ener_old=ener_old + Vrad ;
                                   }
                                   else if( dist < rangeP_aux) {
                                       if(Vangtor_max == -1.0) {
                                          ener_old=ener_old + 0.0  ;
                                       }
                                       else {
                                          Vangtor=Vangtor_max;
                                          FX=(1.0 + tanh(a*(-XOP_aux + dist)))/2.0;
                                          FXI=(1.0 - tanh(a*(-XOP_aux + dist)))/2.0;
                                          vlj(dist, sigma_LJ_aux, VL0_aux, Vrad);
                                          ener_old=ener_old + Vrad * (FXI + Vangtor * FX);
                                       }
                                   }

                                   // new configuration //
                                   if( distn < 0.98000000*sigma_LJ_aux) {
                                       vlj(distn, sigma_LJ_aux, VL0_aux, Vrad);
                                       ener_new=ener_new + Vrad ;
                                   }
                                   else if( distn < rangeP_aux) {
                                       if(Vangtor_max_n == -1.0) {
                                         ener_new=ener_new + 0.0  ;
                                       }
                                       else {
                                         Vangtor=Vangtor_max_n;
                                         FX=(1.0 + tanh(a*(-XOP_aux + distn)))/2.0;
                                         FXI=(1.0 - tanh(a*(-XOP_aux + distn)))/2.0;
                                         vlj(distn, sigma_LJ_aux, VL0_aux, Vrad);
                                         ener_new=ener_new + Vrad * (FXI + Vangtor * FX);
                                       }
                                   }

		              }  // if dist < range_KF
		           }  //  for l < Nat_neigh_cell
	               } // if neigh cell is not empty 

                       entot_new[indx_thread]=ener_new;
                       entot_old[indx_thread]=ener_old;
                       overlap_aux[indx_thread]=overlap;

	        } // if icell < 26 
         }  //  endif index_threads ==0 
////////////////////// end of job of remaining threads /////////////////////////////////////////
         __syncthreads();

////////////////////// only in the master thread /////////////////////////////////////////

         if (indx_thread == 0) {

               // if overlap eq true, reject move
               overlap=0;
               for ( i=0; i < 27; i++) {
                   if( overlap_aux[i] == 1) {
                       overlap=1; }
               }
               if ( overlap == 1 ) {
	             //printf("\n overlap 2");
	             //continue;   // cycles to the next iteration in the loop
               }
               // otherwise, calculate acceptance probability and acept or reject accordingly
               else {
    
                  // collect results for each thread //

                  ener_old=0.0;
                  ener_new=0.0;
                  for ( i=0; i < 27; i++) {
                      ener_old=ener_old+entot_old[i];
                      ener_new=ener_new+entot_new[i];
                  }
                  //
                  de1 = ener_new - ener_old ;
                  prob = __expf ( - beta* de1);

                  if ( curand_uniform(&RND_state[Indx_Rnd]) < prob ) {   // accept
	
	              Sphere_Cell[imol] = Dmove_x ;
	              Sphere_Cell[imol+1] = Dmove_y ;
	              Sphere_Cell[imol+2] = Dmove_z ;
                      Sphere_Cell[imol+3] = Qmove_0 ;
                      Sphere_Cell[imol+4] = Qmove_1 ;
                      Sphere_Cell[imol+5] = Qmove_2 ;
                      Sphere_Cell[imol+6] = Qmove_3 ;

	              DE = DE + de1 ;
	              if ( itras ) {
	                 ntrasa = ntrasa + 1 ;}
	              else {
	                 nrota = nrota +1 ; }
	 
                  }
               }
         }
////////////////////// end of job of the master thread //////////////////////////////////////
         __syncthreads();

      }
         __syncthreads();


}  // endfor Nmove trials

         __syncthreads();
} //endif, chqueo que la celda elegida no esta vacia
         __syncthreads();

if (indx_thread == 0 ) {

  Delta_energy[Indx_Cell] = DE;
  Ntras_ac[Indx_Cell] = ntrasa;
  Nrot_ac[Indx_Cell] = nrota;  }

}  // endif, chequeo que la celda elegida esta en los limites

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////



extern "C" void subsweepgpu_(int *bool_tor, int *seed, int *No_tras_ac, int *No_rot_ac, int *Ndim, int *Nmove, 
                             int *Nat_Cell_Max, int *Ntot_Cell, int *Nsubset_Cell, 
                             int *Npatches_max, int *Ntor_max, int *Npart_types,
                             int *Npatches, float *patch, float *ref_tor_vec, 
                             float *sigma_jon, float *sigma_tor, float *sigma_LJ, float *rangeP, 
                             float *VL0, float *XOP,
                             int *Ntor_angle, float *tor_angle, float *Vpot_Matrix,
                             float *h, float *beta, float *hmax, float *omax, float  *Sphere_Cell, int *List_Cell,
                             int *Map_Cell, int *Nl_Cell, int *Nat_Cell, float *DR_Cell, 
                             double *W_RU, int *Offset_CB, float *En_tot){

  int i, Ndim_dev, Nat_Cell_Max_dev, Ntot_Cell_dev;
  int Nsubset_Cell_dev, Nmove_dev;
  int Npatches_max_dev, Ntor_max_dev, Ntypes_dev;
  int bool_tor_dev;
  unsigned long long seed_dev;

  static curandState_t *RND_States;
  static int init_subsweep=0;

  float beta_dev, hmax_dev, omax_dev;
  float sigma_tor_dev;


  Ndim_dev = *Ndim;
  Nat_Cell_Max_dev = *Nat_Cell_Max;
  Ntot_Cell_dev = *Ntot_Cell;
  Npatches_max_dev = *Npatches_max;
  Ntor_max_dev = *Ntor_max;
  Ntypes_dev = *Npart_types;
  bool_tor_dev = *bool_tor;

  Nsubset_Cell_dev = *Nsubset_Cell;
  Nmove_dev = *Nmove;
  beta_dev = *beta;
  hmax_dev = *hmax;
  omax_dev = *omax;
  seed_dev = *seed; 
  sigma_tor_dev = *sigma_tor;

  if(init_subsweep==0) {


     /* Allocate memory for one state per cell */
     int nblock = 10*Nsubset_Cell_dev/NTHREAD; // i add an extra factor 10, to account for possible volume changes
     if( nblock == 0) { nblock=1; }
     int Ntot_jobs = nblock* NTHREAD;
     if( Ntot_jobs < Nsubset_Cell_dev ) {
       nblock = nblock+1;
     }
     int TOTAL_THREADS = NTHREAD * nblock;

     MANEJA_ERROR(cudaMalloc((void **)&RND_States, TOTAL_THREADS * sizeof(curandState_t)));
  
     unsigned long  long offset = 0;

     setup_random<<<nblock,NTHREAD>>>(Nl_Cell_dev, seed_dev, offset, RND_States);

     /* Allocate memory for variables */


     MANEJA_ERROR(cudaMalloc((void**)&W_RU_dev, Ndim_dev*sizeof(double)));
     MANEJA_ERROR(cudaMalloc((void**)&Offset_dev, Ndim_dev*sizeof(int)));
     MANEJA_ERROR(cudaMalloc((void**)&Ntras_ac_dev, Ntot_Cell_dev*sizeof(int)));
     MANEJA_ERROR(cudaMalloc((void**)&Nrot_ac_dev, Ntot_Cell_dev*sizeof(int)));
     MANEJA_ERROR(cudaMalloc((void**)&Delta_E_dev, Ntot_Cell_dev*sizeof(float)));
     MANEJA_ERROR(cudaMalloc((void**)&overlap_error_dev,Ntot_Cell_dev*sizeof(int)));

     MANEJA_ERROR(cudaMemcpy( W_RU_dev, W_RU, Ndim_dev*sizeof(double), cudaMemcpyHostToDevice));

     Delta_E = (float*)malloc(Ntot_Cell_dev*sizeof(float));
     Ntras_ac = (int*)malloc(Ntot_Cell_dev*sizeof(int));
     Nrot_ac = (int*)malloc(Ntot_Cell_dev*sizeof(int));
     overlap_error = (int*)malloc(Ntot_Cell_dev*sizeof(int));

     init_subsweep=1;

  }

  MANEJA_ERROR(cudaMemcpy( Offset_dev, Offset_CB, Ndim_dev*sizeof(int), cudaMemcpyHostToDevice));

  MANEJA_ERROR(cudaMemcpy( Sphere_Celldev, Sphere_Cell, Ntot_Cell_dev*Nat_Cell_Max_dev*(Ndim_dev+4)*sizeof(float), cudaMemcpyHostToDevice));
  MANEJA_ERROR(cudaMemcpy( List_Cell_dev, List_Cell, Ntot_Cell_dev*Nat_Cell_Max_dev*2*sizeof(int), cudaMemcpyHostToDevice));
  MANEJA_ERROR(cudaMemcpy( Nat_Cell_dev, Nat_Cell, Ntot_Cell_dev*sizeof(int), cudaMemcpyHostToDevice));

  MANEJA_ERROR(cudaMemset(overlap_error_dev, 0, Ntot_Cell_dev*sizeof(int)));
  MANEJA_ERROR(cudaMemset(Delta_E_dev, 0.0, Ntot_Cell_dev*sizeof(float)));
  MANEJA_ERROR(cudaMemset(Ntras_ac_dev, 0, Ntot_Cell_dev*sizeof(int)));
  MANEJA_ERROR(cudaMemset(Nrot_ac_dev, 0, Ntot_Cell_dev*sizeof(int)));
  MANEJA_ERROR(cudaMemcpy( hdev, h, Ndim_dev*Ndim_dev*sizeof(float), cudaMemcpyHostToDevice));  

  int nblock = Nsubset_Cell_dev;
  if( nblock == 0) { nblock=1; }
  int Ntot_jobs = nblock * NTHREAD;
  if( Ntot_jobs < Nsubset_Cell_dev ) {
       nblock = nblock+1;
  }

  subsweep<<<nblock,NTHREAD>>>( RND_States, overlap_error_dev, bool_tor_dev, Ndim_dev, Nmove_dev, Nat_Cell_Max_dev, 
			    Ntot_Cell_dev, Nsubset_Cell_dev, Npatches_max_dev, Ntor_max_dev, Ntypes_dev,
                            Npatches_dev, patchdev, ref_tor_dev, sigma_jon_dev, sigma_tor_dev, sigma_LJ_dev, range_dev, 
                            VL0_dev, XOP_dev, Ntor_angle_dev, tor_angle_dev, Vpot_Matrix_dev, 
			    hdev, beta_dev, hmax_dev, omax_dev, Sphere_Celldev, List_Cell_dev, Map_Cell_dev, 
			    Nl_Cell_dev, Nat_Cell_dev, DR_Cell_dev, Offset_dev, W_RU_dev, 
			    Delta_E_dev, Ntras_ac_dev, Nrot_ac_dev);


  MANEJA_ERROR(cudaMemcpy( Delta_E, Delta_E_dev, Ntot_Cell_dev*sizeof(float), cudaMemcpyDeviceToHost ) );
  MANEJA_ERROR(cudaMemcpy( Ntras_ac, Ntras_ac_dev, Ntot_Cell_dev*sizeof(int), cudaMemcpyDeviceToHost ) );
  MANEJA_ERROR(cudaMemcpy( Nrot_ac, Nrot_ac_dev, Ntot_Cell_dev*sizeof(int), cudaMemcpyDeviceToHost ) );
  MANEJA_ERROR(cudaMemcpy( Sphere_Cell, Sphere_Celldev, Ntot_Cell_dev*Nat_Cell_Max_dev*(Ndim_dev+4)*sizeof(float), cudaMemcpyDeviceToHost ) );
  MANEJA_ERROR(cudaMemcpy( overlap_error, overlap_error_dev, Ntot_Cell_dev*sizeof(int), cudaMemcpyDeviceToHost ) );

  /* actualize energy and number of accepted moves */
  float Total_energy = 0;
  int No_tras_ac_subsweep =0;
  int No_rot_ac_subsweep =0;
  for ( i=0; i < Ntot_Cell_dev; i++) {
     if(overlap_error[i] == 1) {
       printf("\n Error: overlap in old configuration");
       exit(1);
     }
     Total_energy = Total_energy + Delta_E[i];
     No_tras_ac_subsweep = No_tras_ac_subsweep + Ntras_ac[i];
     No_rot_ac_subsweep = No_rot_ac_subsweep + Nrot_ac[i];
  }
  *En_tot = *En_tot + Total_energy;
  *No_tras_ac = *No_tras_ac + No_tras_ac_subsweep;
  *No_rot_ac = *No_rot_ac + No_rot_ac_subsweep;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void shiftGPU(int Ndim, int Nat_Cell_Max, int Ntot_Cell, int *Nl_Cell, float *Sphere_Cell, 
                         float *Sphere_Cell_new, int *List_Cell, int *List_Cell_new, 
                         int *Nat_Cell, int *Nat_Cell_new, double *W_RU, int *ivec, float *disp)
{

  int Nl0, Nl1, Nl2, Indx_Cell, Nat_Cell_ijk;
  int Indice, i, Indx, Indice2, Nat_Neigh;
  int Neigh_x, Neigh_y, Neigh_z, Indx_Cell2;
  int ICell_Neigh, j, Nat_Cell_ijk_new, jcell_x, jcell_y, jcell_z;
  int Indicep, Indx2, Indice3, ICell_Neigh2;

  float Dx, Dy, Dz, Q0, Q1, Q2, Q3;
  double WR0, WR1, WR2;

  Indx_Cell = threadIdx.x + blockIdx.x * blockDim.x;
 
  if (Indx_Cell < Ntot_Cell) {

  Nl0 = Nl_Cell[0];
  Nl1 = Nl_Cell[1];
  Nl2 = Nl_Cell[2];

  jcell_z = Indx_Cell/(Nl0*Nl1)  ;
  jcell_y = ((Indx_Cell - Nl0*Nl1*jcell_z)/Nl0)  ;
  if (jcell_y == 0 && jcell_z == 0) {
     jcell_x = Indx_Cell ;
 }
  else {
     jcell_x = ( Indx_Cell % (jcell_y*Nl0+jcell_z*Nl0*Nl1)) ;
  }


  if ( jcell_x < Nl0  && jcell_y < Nl1 && jcell_z < Nl2 ) {

      Indx_Cell = ( (jcell_x + Nl0) % Nl0) + ((jcell_y+Nl1) % Nl1) * Nl0 + ((jcell_z+Nl2) % Nl2) * Nl0*Nl1;

      Indice = Indx_Cell*Nat_Cell_Max*(Ndim+4);
      Indicep = Indx_Cell*Nat_Cell_Max*2;

      Nat_Cell_ijk = Nat_Cell[Indx_Cell];
      Nat_Cell_ijk_new = 0;

      WR0=W_RU[0]*0.5;
      WR1=W_RU[1]*0.5;
      WR2=W_RU[2]*0.5;

      if ( Nat_Cell_ijk > 0 ) {

          for (  i=0; i < Nat_Cell_ijk; i++ ) {
              
              Indx = Indice + i*(Ndim+4);
              Indx2 = Indicep + i*2;
              if ( Indx+6 <  Ntot_Cell*Nat_Cell_Max*(Ndim+4)  && Indx2+1 < Ntot_Cell*Nat_Cell_Max*2 ) {
              Dx= Sphere_Cell[ Indx ];
              Dy= Sphere_Cell[ Indx+1 ];
              Dz= Sphere_Cell[ Indx+2 ];
              Q0= Sphere_Cell[ Indx+3 ];
              Q1= Sphere_Cell[ Indx+4 ];
              Q2= Sphere_Cell[ Indx+5 ];
              Q3= Sphere_Cell[ Indx+6 ];
              Dx = Dx - ivec[0]*disp[0];
              Dy = Dy - ivec[1]*disp[1];
              Dz = Dz - ivec[2]*disp[2];
              if( ( Dx >= -WR0 ) &&  ( Dx < WR0 ) && ( Dy >= -WR1 ) &&  ( Dy < WR1 ) && ( Dz >= -WR2 ) &&  ( Dz < WR2 ) ) {
                   Nat_Cell_ijk_new=Nat_Cell_ijk_new+1;
                   if (Nat_Cell_ijk_new > Nat_Cell_Max) {
                        printf("\n ERROR in SHIFT Kernel");
                   }
                   Indice2 = Indice + (Nat_Cell_ijk_new-1)*(Ndim+4);
                   if( Indice2+6 > Ntot_Cell*Nat_Cell_Max*(Ndim+4) ) {
                      printf("\n error in shiftcell 1 %d",Indice2); }
                   else {
                     Sphere_Cell_new[ Indice2 ] = Dx;
                     Sphere_Cell_new[ Indice2 + 1 ] = Dy;
                     Sphere_Cell_new[ Indice2 + 2 ] = Dz;
                     Sphere_Cell_new[ Indice2 + 3 ] = Q0;
                     Sphere_Cell_new[ Indice2 + 4 ] = Q1;
                     Sphere_Cell_new[ Indice2 + 5 ] = Q2;
                     Sphere_Cell_new[ Indice2 + 6 ] = Q3;
                     Indice3 = Indicep + (Nat_Cell_ijk_new-1)*2;
                     List_Cell_new[ Indice3 ] = List_Cell[ Indx2 ];
                     List_Cell_new[ Indice3 +1 ] = List_Cell[ Indx2 + 1];

                   }
 //                           printf("\n Indice 2 %d: %f",Indice2+6,Sphere_Cell_new[ Indice2+6 ]);
              }
          }
         } // if Nat_Cell_ijk > 0
         }
          // checking also the neighbour cell in the direction of ivec
          Neigh_x = jcell_x + ivec[0];
          Neigh_y = jcell_y + ivec[1];
          Neigh_z = jcell_z + ivec[2];

          Indx_Cell2 = ( (Neigh_x + Nl0) % Nl0) + ((Neigh_y+Nl1) % Nl1) * Nl0 + ((Neigh_z+Nl2) % Nl2) * Nl0*Nl1;

          ICell_Neigh= Indx_Cell2*Nat_Cell_Max*(Ndim+4);
          ICell_Neigh2= Indx_Cell2*Nat_Cell_Max*2;

          if (Indx_Cell2 > Ntot_Cell-1) {
               printf("\n error in shiftcell X %d",Indx_Cell2); }
          else {
               Nat_Neigh = Nat_Cell[ Indx_Cell2 ];
          }

          if ( Nat_Neigh > 0) {

             for (j=0; j < Nat_Neigh; j++) {

                Indx = ICell_Neigh + j*(Ndim+4);
                Indx2 = ICell_Neigh2 + j*2;

                if( Indx+6 > Ntot_Cell*Nat_Cell_Max*(Ndim+4)-1 ) {
                     printf("\n error in shiftcell 2 %d",Indx); }
                else {
                   Dx= Sphere_Cell[ Indx ];
                   Dy= Sphere_Cell[ Indx+1 ];
                   Dz= Sphere_Cell[ Indx+2 ];

                   Q0= Sphere_Cell[ Indx+3 ];
                   Q1= Sphere_Cell[ Indx+4 ];
                   Q2= Sphere_Cell[ Indx+5 ];
                   Q3= Sphere_Cell[ Indx+6 ];

                   Dx = Dx - ivec[0]*disp[0];
                   Dy = Dy - ivec[1]*disp[1];
                   Dz = Dz - ivec[2]*disp[2];
 
                   if( ( Dx >= -WR0 ) &&  ( Dx < WR0 ) && ( Dy >= -WR1 ) &&  ( Dy < WR1 ) && ( Dz >= -WR2 ) &&  ( Dz < WR2 ) ) {
                   } /* particle stays in the cell: do nothing */ 
                   else {

                       Dx = Dx + ivec[0]*W_RU[0];
                       Dy = Dy + ivec[1]*W_RU[1];
                       Dz = Dz + ivec[2]*W_RU[2];
                       Nat_Cell_ijk_new=Nat_Cell_ijk_new+1;
                       if (Nat_Cell_ijk_new > Nat_Cell_Max) {
                             printf("\n ERROR: number of atoms in cell larger than Nat_Cell_Max %d %d",Nat_Cell_ijk_new,Nat_Cell_Max);
                       } else  {

                           Indice2 = Indice + (Nat_Cell_ijk_new-1)*(Ndim+4);

                           if( Indice2+6 > Ntot_Cell*Nat_Cell_Max*(Ndim+4)-1 ) {
                                     printf("\n error in shiftcell 3 %d",Indx);}
                           else {

                            Sphere_Cell_new[ Indice2 ] = Dx;
                            Sphere_Cell_new[ Indice2 + 1 ] = Dy;
                            Sphere_Cell_new[ Indice2 + 2 ] = Dz;

                            Sphere_Cell_new[ Indice2 + 3 ] = Q0;
                            Sphere_Cell_new[ Indice2 + 4 ] = Q1;
                            Sphere_Cell_new[ Indice2 + 5 ] = Q2;
                            Sphere_Cell_new[ Indice2 + 6 ] = Q3;

                            Indice3 = Indicep + (Nat_Cell_ijk_new-1)*2;

                            List_Cell_new[ Indice3 ] = List_Cell[ Indx2 ];
                            List_Cell_new[ Indice3 +1 ] = List_Cell[ Indx2 + 1];

                           }

                       }
                }
                }
              }
              } 

              if (Indx_Cell > Ntot_Cell ) {
                     printf("\n error in shiftcell 4 %d %d",Indice, Indx_Cell);}
              else {
                     Nat_Cell_new[Indx_Cell] = Nat_Cell_ijk_new;
              }
  }  // if jcell_x, jcell_y, jcell_z within bounds

}
     //__syncthreads();



}

////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

extern "C" void shiftcellsgpu_(int *Ndim, int *Nat_Cell_Max, int *Ntot_Cell, int *Nl_Cell, 
                               float *Sphere_Cell, int *List_Cell, int *Nat_Cell, double *W_RU, int *ivec, 
                               float *disp) {

  int Ndim_dev, Nat_Cell_Max_dev, Ntot_Cell_dev;
  static int *ivec_dev;

  static float *disp_dev;
  static int init_shift =0;

  Ndim_dev = *Ndim;
  Nat_Cell_Max_dev = *Nat_Cell_Max;
  Ntot_Cell_dev = *Ntot_Cell;

  if ( init_shift == 0 ) {

     MANEJA_ERROR(cudaMalloc((void**)&Sphere_Cell_new, Ntot_Cell_dev*Nat_Cell_Max_dev*(Ndim_dev+4)*sizeof(float)));
     MANEJA_ERROR(cudaMalloc((void**)&List_Cell_new, Ntot_Cell_dev*Nat_Cell_Max_dev*2*sizeof(int)));
     MANEJA_ERROR(cudaMalloc((void**)&Nat_Cell_new, Ntot_Cell_dev*sizeof(int)));
     MANEJA_ERROR(cudaMalloc((void**)&ivec_dev, Ndim_dev*sizeof(int)));
     MANEJA_ERROR(cudaMalloc((void**)&disp_dev, Ndim_dev*sizeof(float)));

     init_shift = 1;
  }

  MANEJA_ERROR(cudaMemset( Sphere_Cell_new, -10., Ntot_Cell_dev*Nat_Cell_Max_dev*(Ndim_dev+4)*sizeof(float)));
  MANEJA_ERROR(cudaMemset( List_Cell_new, 0, Ntot_Cell_dev*Nat_Cell_Max_dev*2*sizeof(int)));
  MANEJA_ERROR(cudaMemset( Nat_Cell_new, 0, Ntot_Cell_dev*sizeof(int)));  

  MANEJA_ERROR(cudaMemcpy( Sphere_Celldev, Sphere_Cell, Ntot_Cell_dev*Nat_Cell_Max_dev*(Ndim_dev+4)*sizeof(float), cudaMemcpyHostToDevice));
  MANEJA_ERROR(cudaMemcpy( List_Cell_dev, List_Cell, Ntot_Cell_dev*Nat_Cell_Max_dev*2*sizeof(int), cudaMemcpyHostToDevice));
  MANEJA_ERROR(cudaMemcpy( Nat_Cell_dev, Nat_Cell, Ntot_Cell_dev*sizeof(int), cudaMemcpyHostToDevice));  
  MANEJA_ERROR(cudaMemcpy( ivec_dev, ivec, Ndim_dev*sizeof(int), cudaMemcpyHostToDevice));  
  MANEJA_ERROR(cudaMemcpy( disp_dev, disp, Ndim_dev*sizeof(float), cudaMemcpyHostToDevice));  
  MANEJA_ERROR(cudaMemcpy( W_RU_dev, W_RU, Ndim_dev*sizeof(double), cudaMemcpyHostToDevice));

  // Call GPU kernel for energy calculation
  int nblock = Ntot_Cell_dev/NTHREAD;
  if ( nblock == 0 ) { nblock=1; }
  int Ntot_jobs = nblock* NTHREAD;
  if( Ntot_jobs < Ntot_Cell_dev ) {
       nblock = nblock+1;
  }
  Ntot_jobs = nblock* NTHREAD;

   shiftGPU<<<nblock,NTHREAD>>>( Ndim_dev, Nat_Cell_Max_dev, Ntot_Cell_dev, Nl_Cell_dev, Sphere_Celldev, 
                                   Sphere_Cell_new, List_Cell_dev, List_Cell_new, Nat_Cell_dev, 
                                   Nat_Cell_new, W_RU_dev, ivec_dev, disp_dev);

  MANEJA_ERROR(cudaMemcpy( Sphere_Cell, Sphere_Cell_new, Ntot_Cell_dev*Nat_Cell_Max_dev*(Ndim_dev+4)*sizeof(float), cudaMemcpyDeviceToHost ) );
  MANEJA_ERROR(cudaMemcpy( List_Cell, List_Cell_new, Ntot_Cell_dev*Nat_Cell_Max_dev*2*sizeof(int), cudaMemcpyDeviceToHost ) );
  MANEJA_ERROR(cudaMemcpy( Nat_Cell, Nat_Cell_new, Ntot_Cell_dev*sizeof(int), cudaMemcpyDeviceToHost));  

}


extern "C" void clean_() {

  MANEJA_ERROR(cudaFree(Map_Cell_dev));
  MANEJA_ERROR(cudaFree(Offset_dev));
  MANEJA_ERROR(cudaFree(Nl_Cell_dev));
  MANEJA_ERROR(cudaFree(Nat_Cell_dev));
  MANEJA_ERROR(cudaFree(Npatches_dev));
  MANEJA_ERROR(cudaFree(Ntor_angle_dev));
  MANEJA_ERROR(cudaFree(List_Cell_dev));

  MANEJA_ERROR(cudaFree(W_RU_dev));
  MANEJA_ERROR(cudaFree(patchdev));
  MANEJA_ERROR(cudaFree(hdev));
  MANEJA_ERROR(cudaFree(DR_Cell_dev));
  MANEJA_ERROR(cudaFree(Sphere_Celldev));
  MANEJA_ERROR(cudaFree(energydev));
  MANEJA_ERROR(cudaFree(ref_tor_dev));
  MANEJA_ERROR(cudaFree(tor_angle_dev));
  MANEJA_ERROR(cudaFree(Vpot_Matrix_dev));

}
///////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////


extern "C" void limpio_(  ) {
   cudaError_t reset = cudaDeviceReset( ) ;
}


///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
extern "C" void energyinitialize_(int *Ndim, int *Nat_Cell_Max, int *Ntot_Cell,
                          float *h, float *Sphere_Cell, int *List_Cell, int *Map_Cell,
                          int *Nl_Cell, int *Nat_Cell, float *DR_Cell, double *W_RU) {

  int Ndim_dev, Nat_Cell_Max_dev, Ntot_Cell_dev;

  Ndim_dev = *Ndim;
  Nat_Cell_Max_dev = *Nat_Cell_Max;
  Ntot_Cell_dev = *Ntot_Cell;

  //printf("\n updating Ntot_Cell %d in GPU\n",Ntot_Cell_dev);
  MANEJA_ERROR(cudaFree(DR_Cell_dev));
  MANEJA_ERROR(cudaFree(Map_Cell_dev));
  MANEJA_ERROR(cudaFree(Sphere_Celldev));
  MANEJA_ERROR(cudaFree(Nat_Cell_dev));
  MANEJA_ERROR(cudaFree(Nl_Cell_dev));
  MANEJA_ERROR(cudaFree(List_Cell_dev));
  MANEJA_ERROR(cudaFree(Sphere_Cell_new));
  MANEJA_ERROR(cudaFree(List_Cell_new));
  MANEJA_ERROR(cudaFree(Nat_Cell_new));

  MANEJA_ERROR(cudaFree(Ntras_ac_dev));
  MANEJA_ERROR(cudaFree(Nrot_ac_dev));
  MANEJA_ERROR(cudaFree(Delta_E_dev));
  MANEJA_ERROR(cudaFree(overlap_error_dev));
  MANEJA_ERROR(cudaFree(energydev));

  free(Ntras_ac);
  free(Nrot_ac);
  free(Delta_E);
  free(overlap_error);
  free(energy);

  MANEJA_ERROR(cudaMalloc((void**)&DR_Cell_dev, Ntot_Cell_dev*26*3*sizeof(float)));
  MANEJA_ERROR(cudaMalloc((void**)&Map_Cell_dev, Ntot_Cell_dev*26*sizeof(int)));
  MANEJA_ERROR(cudaMalloc((void**)&Sphere_Celldev, Ntot_Cell_dev*Nat_Cell_Max_dev*(Ndim_dev+4)*sizeof(float)));
  MANEJA_ERROR(cudaMalloc((void**)&Nat_Cell_dev, Ntot_Cell_dev*sizeof(int)));
  MANEJA_ERROR(cudaMalloc((void**)&Nl_Cell_dev, Ndim_dev*sizeof(int)));
  MANEJA_ERROR(cudaMalloc((void**)&List_Cell_dev, Ntot_Cell_dev*Nat_Cell_Max_dev*2*sizeof(int)));

  MANEJA_ERROR(cudaMalloc((void**)&Sphere_Cell_new, Ntot_Cell_dev*Nat_Cell_Max_dev*(Ndim_dev+4)*sizeof(float)));
  MANEJA_ERROR(cudaMalloc((void**)&List_Cell_new, Ntot_Cell_dev*Nat_Cell_Max_dev*2*sizeof(int)));
  MANEJA_ERROR(cudaMalloc((void**)&Nat_Cell_new, Ntot_Cell_dev*sizeof(int)));

  MANEJA_ERROR(cudaMalloc((void**)&Ntras_ac_dev, Ntot_Cell_dev*sizeof(int)));
  MANEJA_ERROR(cudaMalloc((void**)&Nrot_ac_dev, Ntot_Cell_dev*sizeof(int)));
  MANEJA_ERROR(cudaMalloc((void**)&Delta_E_dev, Ntot_Cell_dev*sizeof(float)));
  MANEJA_ERROR(cudaMalloc((void**)&overlap_error_dev,Ntot_Cell_dev*sizeof(int)));
  MANEJA_ERROR(cudaMalloc((void**)&energydev, Ntot_Cell_dev*sizeof(float)));

  energy = (float*)malloc(Ntot_Cell_dev*sizeof(float));
  Delta_E = (float*)malloc(Ntot_Cell_dev*sizeof(float));
  Ntras_ac = (int*)malloc(Ntot_Cell_dev*sizeof(int));
  Nrot_ac = (int*)malloc(Ntot_Cell_dev*sizeof(int));
  overlap_error = (int*)malloc(Ntot_Cell_dev*sizeof(int));

  MANEJA_ERROR(cudaMemcpy( DR_Cell_dev, DR_Cell, Ntot_Cell_dev*26*3*sizeof(float), cudaMemcpyHostToDevice));
  MANEJA_ERROR(cudaMemcpy( Map_Cell_dev, Map_Cell, Ntot_Cell_dev*26*sizeof(int), cudaMemcpyHostToDevice));
  MANEJA_ERROR(cudaMemcpy( Nl_Cell_dev, Nl_Cell, Ndim_dev*sizeof(int), cudaMemcpyHostToDevice));

  MANEJA_ERROR(cudaMemcpy( hdev, h, Ndim_dev*Ndim_dev*sizeof(float), cudaMemcpyHostToDevice));
  MANEJA_ERROR(cudaMemcpy( W_RU_dev, W_RU, Ndim_dev*sizeof(double), cudaMemcpyHostToDevice));

  MANEJA_ERROR(cudaMemcpy( List_Cell_dev, List_Cell, Ntot_Cell_dev*Nat_Cell_Max_dev*2*sizeof(int), cudaMemcpyHostToDevice));
  MANEJA_ERROR(cudaMemcpy( Sphere_Celldev, Sphere_Cell, Ntot_Cell_dev*Nat_Cell_Max_dev*(Ndim_dev+4)*sizeof(float), cudaMemcpyHostToDevice));
  MANEJA_ERROR(cudaMemcpy( Nat_Cell_dev, Nat_Cell, Ntot_Cell_dev*sizeof(int), cudaMemcpyHostToDevice));
}

