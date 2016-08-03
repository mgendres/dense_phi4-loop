/*  Number of time lattice sites. */
#define NT 2
/*  Number of space lattice sites */
#define NS 2
/* Ratio of time and space lattice constants */
#define B 1.0
/*  Particle mass  */
#define MSq -26.05
/*  Chemical potential  */
#define MU 0.1
/*  Quartic coupling  */
#define LAMBDA 192.0
/*  Initial system charge  */
#define QINIT 0
/*  Initial value for radial mode  */
#define SINIT 0.5
/*  Number of iterations  */
#define NIT 100000000
/*  Number of iterations before equilibrium.
 *  Data collection begins one counter reaches this value.  */
#define NEQ 1
/* Estimate of correlation time.
 *  Data is recorded once every TCORR steps.  */ 
#define TCORR 100
/*  Maximum step size of random walker  */
#define D 0.45
