/***********************************************/
/* LOOP ALGORITHM FOR BOSONS AT FINITE DENSITY */
/* worm_v7.c                                   */
/* AUTHOR: MICHAEL G. ENDRES *******************/
/* LAST MODIFIED: 05/31/06 *********************/
/***********************************************/

#include <stdio.h>
#include <math.h>
#include <time.h>
#include "bessi.h"
#include "ran2.h"
#include "ran3.h"
#include "irbit2.h"
#include "mod.h"
#include "parameters.h"

#define RNG ran3

/***********************/
/* PROBABILITY WEIGHTS */
/***********************/
float wS(int j, int k, int n, float s[NT][NS][NS], int l[NT][NS][NS][3], float ds)
{
  float ans=0.0;
  float sTrial;

  sTrial=s[j][k][n]+ds;
  if (sTrial<=0.0) {
    return ans;
  } else {
    ans += sTrial/s[j][k][n];
    ans *= exp(-B*(2.0+1.0/(B*B)+MSq/2.0)*(sTrial*sTrial-s[j][k][n]*s[j][k][n]));
    ans *= exp(-B*(LAMBDA/16.0)*(sTrial*sTrial*sTrial*sTrial-s[j][k][n]*s[j][k][n]*s[j][k][n]*s[j][k][n]));
    ans *= bessi(l[j][k][n][0],sTrial*s[mod(j+1,NT)][k][n]/B);
    ans /= bessi(l[j][k][n][0],s[j][k][n]*s[mod(j+1,NT)][k][n]/B);
    ans *= bessi(l[j][k][n][1],sTrial*s[j][mod(k+1,NS)][n]*B);
    ans /= bessi(l[j][k][n][1],s[j][k][n]*s[j][mod(k+1,NS)][n]*B);
    ans *= bessi(l[j][k][n][2],sTrial*s[j][k][mod(n+1,NS)]*B);
    ans /= bessi(l[j][k][n][2],s[j][k][n]*s[j][k][mod(n+1,NS)]*B);
    ans *= bessi(l[mod(j-1,NT)][k][n][0],sTrial*s[mod(j-1,NT)][k][n]/B);
    ans /= bessi(l[mod(j-1,NT)][k][n][0],s[j][k][n]*s[mod(j-1,NT)][k][n]/B);
    ans *= bessi(l[j][mod(k-1,NS)][n][1],sTrial*s[j][mod(k-1,NS)][n]*B);
    ans /= bessi(l[j][mod(k-1,NS)][n][1],s[j][k][n]*s[j][mod(k-1,NS)][n]*B);
    ans *= bessi(l[j][k][mod(n-1,NS)][2],sTrial*s[j][k][mod(n-1,NS)]*B);
    ans /= bessi(l[j][k][mod(n-1,NS)][2],s[j][k][n]*s[j][k][mod(n-1,NS)]*B);
  }
  return ans;
};

float wPlaqY(float s, float spt, float spx, float sptpx, int lpt, int lpx, int lptpx, int lpxpt, int dQ)
{
  float ans=1.0;

  ans *= bessi(lpt+dQ,s*spt/B);
  ans /= bessi(lpt,s*spt/B);
  ans *= bessi(lpx-dQ,s*spx*B);
  ans /= bessi(lpx,s*spx*B);
  ans *= bessi(lptpx+dQ,spt*sptpx*B);
  ans /= bessi(lptpx,spt*sptpx*B);
  ans *= bessi(lpxpt-dQ,spx*sptpx/B);
  ans /= bessi(lpxpt,spx*sptpx/B);
  return ans;
};

float wPlaqX(float s, float spy, float spt, float spypt, int lpy, int lpt, int lpypt, int lptpy, int dQ)
{
  float ans=1.0;

  ans *= bessi(lpy+dQ,s*spy*B);
  ans /= bessi(lpy,s*spy*B);
  ans *= bessi(lpt-dQ,s*spt/B);
  ans /= bessi(lpt,s*spt/B);
  ans *= bessi(lpypt+dQ,spy*spypt/B);
  ans /= bessi(lpypt,spy*spypt/B);
  ans *= bessi(lptpy-dQ,spt*spypt*B);
  ans /= bessi(lptpy,spt*spypt*B);
  return ans;
};

float wPlaqT(float s, float spx, float spy, float spxpy, int lpx, int lpy, int lpxpy, int lpypx, int dQ)
{
  float ans=1.0;

  ans *= bessi(lpx+dQ,s*spx*B);
  ans /= bessi(lpx,s*spx*B);
  ans *= bessi(lpy-dQ,s*spy*B);
  ans /= bessi(lpy,s*spy*B);
  ans *= bessi(lpxpy+dQ,spx*spxpy*B);
  ans /= bessi(lpxpy,spx*spxpy*B);
  ans *= bessi(lpypx-dQ,spy*spxpy*B);
  ans /= bessi(lpypx,spy*spxpy*B);
  return ans;
};

float wQ(int k, int n, float s[NT][NS][NS], int l[NT][NS][NS][3], int dQ)
{
  float ans=1.0;
  int j;

  for (j=0;j<NT;j++) {
    ans *= bessi(l[j][k][n][0]+dQ,s[j][k][n]*s[mod(j+1,NT)][k][n]/B);
    ans /= bessi(l[j][k][n][0],s[j][k][n]*s[mod(j+1,NT)][k][n]/B);
    ans *= exp(B*MU*dQ);
  };
  return ans;
};

float wJ1(int j, int n, float s[NT][NS][NS], int l[NT][NS][NS][3], int dQ)
{
  float ans=1.0;
  int k;

  for (k=0;k<NS;k++) {
    ans *= bessi(l[j][k][n][1]+dQ,s[j][k][n]*s[j][mod(k+1,NS)][n]*B);
    ans /= bessi(l[j][k][n][1],s[j][k][n]*s[j][mod(k+1,NS)][n]*B);
  };
  return ans;
};

float wJ2(int j, int k, float s[NT][NS][NS], int l[NT][NS][NS][3], int dQ)
{
  float ans=1.0;
  int n;

  for (n=0;n<NS;n++) {
    ans *= bessi(l[j][k][n][2]+dQ,s[j][k][n]*s[j][k][mod(n+1,NS)]*B);
    ans /= bessi(l[j][k][n][2],s[j][k][n]*s[j][k][mod(n+1,NS)]*B);
  };
  return ans;
};

/***************/
/* OBSERVABLES */
/***************/
float s2(float s[NT][NS][NS]) {
  float ans=0.0;
  int j;
  int k;
  int n;

  for (j=0;j<NT;j++) {
    for (k=0;k<NS;k++) {
      for (n=0;n<NS;n++) {
        ans += s[j][k][n]*s[j][k][n];
      };
    };
  };
  ans /= NT*NS*NS;
  return ans;
};

float s4(float s[NT][NS][NS]) {
  float ans=0.0;
  int j;
  int k;
  int n;

  for (j=0;j<NT;j++) {
    for (k=0;k<NS;k++) {
      for (n=0;n<NS;n++) {
        ans += s[j][k][n]*s[j][k][n]*s[j][k][n]*s[j][k][n];
      };
    };
  };
  ans /= NT*NS*NS;
  return ans;
};

/*
 *  Simulation
 */

int main(void)
{
  int i;  /*  MC step counter  */
  int j;  /*  Site counter for t direction */
  int k;  /*  Site counter for x direction  */
  int n;  /*  Site counter for y direction  */
  int r;  /*  Extra counter  */
  long seed;  /*  Random number generator seed  */
  long iseed;  /*  Random bit generator seed  */
  float s[NT][NS][NS];  /*  Site variables  */
  int l[NT][NS][NS][3];  /*  Link variables  */
  float ds;  /*  Trial change in site variable  */
  int dQ;  /*  Trial change in charge  */
  int Q;  /*  Charge of system  */
  int J1;  /*  Global current in x-direction  */
  int J2;  /*  Global charge in y-direction  */
  float w;  /* Transition probability used for sites and links  */
  int rAcc;  /*  Acceptance rate counter  */
  FILE *binaryp;

  /********************/
  /* INITIALIZE SEEDS */
  /********************/
  seed = -time((time_t *)NULL);  /* Random generator */
  iseed = -time((time_t *)NULL);  /* Bit generator */

 /******************************/
  /* INITIALIZE STATE OF SYSTEM */
  /******************************/
  binaryp=fopen("system.conf","rb");
  if (binaryp==NULL) {
    /******************/
    /* USE COLD START */
    /******************/
    for (j=0;j<NT;j++) {
      for (k=0;k<NS;k++) {
        for (n=0;n<NS;n++) {
          s[j][k][n]=SINIT;
          l[j][k][n][0]=0;
          l[j][k][n][1]=0;
          l[j][k][n][2]=0;
        };
      };
      l[j][0][0][0]=QINIT;
    };
    Q=QINIT;
    J1=0;
    J2=0;
  }else{
    /******************************/
    /* USE PREVIOUS CONFIGURATION */
    /******************************/
    fread(&Q, sizeof (int), 1, binaryp);
    fread(&J1, sizeof (int), 1, binaryp);
    fread(&J2, sizeof (int), 1, binaryp);
    fread(s, sizeof (float), NT*NS*NS, binaryp);
    fread(l, sizeof (int), NT*NS*NS*3, binaryp);
    fclose(binaryp);
  };

/***  BEGIN MC SIMULATIONS  ***/
  for (i=0;i<NIT;i++) {
    rAcc=0;  /* Accceptance rate counter  */
    for (j=0;j<NT;j++) {  /*  Sum over t  */
      for (k=0;k<NS;k++) {  /*  Sum over x  */
        for (n=0;n<NS;n++) {  /*  Sum over y  */
/***  UPDATE SITE  ***/
          ds=D*(2.0*RNG(&seed)-1.0);  /*  Generate a trial step  */
          if (wS(j,k,n,s,l,ds)>RNG(&seed)) {  /* Update site variable  */
            s[j][k][n] += ds;
            rAcc++;
          };
/***  UPDATE PLAQUETTE Y  ***/
          dQ=(2*irbit2(&iseed)-1);  /*  Generate trial plaquette step for y oriented plaquettes */
          w=wPlaqY(s[j][k][n],s[mod(j+1,NT)][k][n],s[j][mod(k+1,NS)][n],s[mod(j+1,NT)][mod(k+1,NS)][n],
              l[j][k][n][0],l[j][k][n][1],l[mod(j+1,NT)][k][n][1],l[j][mod(k+1,NS)][n][0],dQ);
          if (w>RNG(&seed)) {  /*  Update plaquette  */
            l[j][k][n][0] += dQ;
            l[j][k][n][1] -= dQ;
            l[mod(j+1,NT)][k][n][1] += dQ;
            l[j][mod(k+1,NS)][n][0] -= dQ;
          };
/***  UPDATE PLAQUETTE T  ***/
          dQ=(2*irbit2(&iseed)-1);  /*  Generate trial plaquette step for t oriented plaquettes */
          w=wPlaqT(s[j][k][n],s[j][mod(k+1,NS)][n],s[j][k][mod(n+1,NS)],s[j][mod(k+1,NS)][mod(n+1,NS)],
            l[j][k][n][1],l[j][k][n][2],l[j][mod(k+1,NS)][n][2],l[j][k][mod(n+1,NS)][1],dQ);
          if (w>RNG(&seed)) {  /*  Update plaquette  */
            l[j][k][n][1] += dQ;
            l[j][k][n][2] -= dQ;
            l[j][mod(k+1,NS)][n][2] += dQ;
            l[j][k][mod(n+1,NS)][1] -= dQ;
          };
/***  UPDATE PLAQUETTE X  ***/
          dQ=(2*irbit2(&iseed)-1);  /*  Generate trial plaquette step for x oriented plaquettes */
          w=wPlaqX(s[j][k][n],s[j][k][mod(n+1,NS)],s[mod(j+1,NT)][k][n],s[mod(j+1,NT)][k][mod(n+1,NS)],
            l[j][k][n][2],l[j][k][n][0],l[j][k][mod(n+1,NS)][0],l[mod(j+1,NT)][k][n][2],dQ);
          if (w>RNG(&seed)) {  /*  Update plaquette  */
            l[j][k][n][2] += dQ;
            l[j][k][n][0] -= dQ;
            l[j][k][mod(n+1,NS)][0] += dQ;
            l[mod(j+1,NT)][k][n][2] -= dQ;
          };
/***  UPDATE TOTAL CHARGE, GLOBAL CURRENTS ***/
          dQ=(2*irbit2(&iseed)-1);
          if (wQ(k,n,s,l,dQ)>RNG(&seed)) {
            Q += dQ;
            for (r=0;r<NT;r++) {
              l[r][k][n][0] += dQ;
            };
          };

          dQ=(2*irbit2(&iseed)-1);
          if (wJ1(j,n,s,l,dQ)>RNG(&seed)) {
            J1 += dQ;
            for (r=0;r<NS;r++) {
              l[j][r][n][1] += dQ;
            };
          };
          dQ=(2*irbit2(&iseed)-1);
          if (wJ2(j,k,s,l,dQ)>RNG(&seed)) {
            J2 += dQ;
            for (r=0;r<NS;r++) {
              l[j][k][r][2] += dQ;
            };
          };

/***  PRINT LOCAL PARAMETERS  ***/
//          printf("\n%f",s[j][k][n]);  /*  Print value of site variable  */
//          printf(" ");  /* Add a space. */
//          printf("\n%d %d %d",l[j][k][n][0],l[j][k][n][1],l[j][k][n][2]);  /*  Print value of link variable  */
//          printf("\n%d",div(j,k,n,l));  /*  Print value of divergence about site (j,k,n)  */

        }  /* End of n for loop  */
      }  /*  End of k for loop  */
    }  /*  End of j for loop  */
/***  PRINT AVERAGES  ***/
    if ((i>NEQ)&&(i%TCORR==0)) {
      printf("%d\n",Q);  /*  Print total charge corresponding to field configuration  */
//      printf("%f\n",s2(s));  /*  Print average |\phi|^2  */
//      printf("%f",Q*1.0/(NS*NS));  /*  Print total charge corresponding to field configuration  */
//      printf(" ");  /* Add a space. */
//      printf("%d",J1);  /*  Print currents  */
//      printf(" ");  /* Add a space. */
//      printf("%d",J2);  /*  Print currents  */
//      printf(" ");  /* Add a space. */
//      printf("%d",rAcc);  /*  Print acceptance rate  */
//      printf("%f",rAcc*1.0/(NT*NS*NS));  /*  Print acceptance rate  */
//      printf("\n");  /* Add new line. */
    };
  };  /* End of i for loop */

  /*****************************/
  /* WRITE SYSTEM STATE T DISK */
  /*****************************/
  binaryp=fopen("system.conf","wb");
  fwrite(&Q, sizeof (int), 1, binaryp);
  fwrite(&J1, sizeof (int), 1, binaryp);
  fwrite(&J2, sizeof (int), 1, binaryp);
  fwrite(s, sizeof (float), NT*NS*NS, binaryp);
  fwrite(l, sizeof (int), NT*NS*NS*3, binaryp);
  fclose(binaryp);
  /************************/
  /* END OF MC SIMULATION */
  /************************/
  return 0;
}
