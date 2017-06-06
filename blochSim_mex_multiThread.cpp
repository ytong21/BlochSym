
#include "mex.h" 
#include <stdio.h>
#include <math.h>
#include <vector>
#include <omp.h>


#define TWOPI	6.28318530717959
#define GAMMA   TWOPI

#define DEBUG

// CTR:
#define DEBUG_printf( ... ) \
	if (debugflag) mexPrintf(__VA_ARGS__); \
	else (void)0
// (The (void)0 gives an error if DEBUG_printf is missing a terminating ;.
static bool debugflag = false;
// End CTR.

void multmatvec(double *mat, double *vec, double *matvec)

	/* Multiply 3x3 matrix by 3x1 vector. */

{
*matvec++ = mat[0]*vec[0] + mat[3]*vec[1] + mat[6]*vec[2];
*matvec++ = mat[1]*vec[0] + mat[4]*vec[1] + mat[7]*vec[2];
*matvec++ = mat[2]*vec[0] + mat[5]*vec[1] + mat[8]*vec[2];
}



void addvecs(double *vec1, double *vec2, double *vecsum)

	/* Add two 3x1 Vectors */

{
*vecsum++ = *vec1++ + *vec2++;
*vecsum++ = *vec1++ + *vec2++;
*vecsum++ = *vec1++ + *vec2++;
}




void adjmat(double *mat, double *adj)

/* ======== Adjoint of a 3x3 matrix ========= */

{
*adj++ = (mat[4]*mat[8]-mat[7]*mat[5]);	
*adj++ =-(mat[1]*mat[8]-mat[7]*mat[2]);
*adj++ = (mat[1]*mat[5]-mat[4]*mat[2]);
*adj++ =-(mat[3]*mat[8]-mat[6]*mat[5]);
*adj++ = (mat[0]*mat[8]-mat[6]*mat[2]);
*adj++ =-(mat[0]*mat[5]-mat[3]*mat[2]);
*adj++ = (mat[3]*mat[7]-mat[6]*mat[4]);
*adj++ =-(mat[0]*mat[7]-mat[6]*mat[1]);
*adj++ = (mat[0]*mat[4]-mat[3]*mat[1]);
}


void zeromat(double *mat)

/* ====== Set a 3x3 matrix to all zeros	======= */

{
*mat++=0;
*mat++=0;
*mat++=0;
*mat++=0;
*mat++=0;
*mat++=0;
*mat++=0;
*mat++=0;
*mat++=0;
}


void eyemat(double *mat)

/* ======== Return 3x3 Identity Matrix  ========= */

{
zeromat(mat);
mat[0]=1;
mat[4]=1;
mat[8]=1;

}

double detmat(double *mat)

/* ======== Determinant of a 3x3 matrix ======== */

{
double det;

det = mat[0]*mat[4]*mat[8];
det+= mat[3]*mat[7]*mat[2];
det+= mat[6]*mat[1]*mat[5];
det-= mat[0]*mat[7]*mat[5];
det-= mat[3]*mat[1]*mat[8];
det-= mat[6]*mat[4]*mat[2];

return det;
}


void scalemat(double *mat, double scalar)

/* ======== multiply a matrix by a scalar ========= */

{
*mat++ *= scalar;
*mat++ *= scalar;
*mat++ *= scalar;
*mat++ *= scalar;
*mat++ *= scalar;
*mat++ *= scalar;
*mat++ *= scalar;
*mat++ *= scalar;
*mat++ *= scalar;
}


void invmat(double *mat, double *imat)

/* ======== Inverse of a 3x3 matrix ========= */
/*	DO NOT MAKE THE OUTPUT THE SAME AS ONE OF THE INPUTS!! */

{
double det;
int count;

det = detmat(mat);	/* Determinant */
adjmat(mat, imat);	/* Adjoint */

for (count=0; count<9; count++)
	*imat = *imat / det;
    imat++;
}


void addmats(double *mat1, double *mat2, double *matsum)

/* ====== Add two 3x3 matrices.	====== */

{
*matsum++ = *mat1++ + *mat2++;
*matsum++ = *mat1++ + *mat2++;
*matsum++ = *mat1++ + *mat2++;
*matsum++ = *mat1++ + *mat2++;
*matsum++ = *mat1++ + *mat2++;
*matsum++ = *mat1++ + *mat2++;
*matsum++ = *mat1++ + *mat2++;
*matsum++ = *mat1++ + *mat2++;
*matsum++ = *mat1++ + *mat2++;
}


void multmats(double *mat1, double *mat2, double *matproduct)

/* ======= Multiply two 3x3 matrices. ====== */
/*	DO NOT MAKE THE OUTPUT THE SAME AS ONE OF THE INPUTS!! */

{
*matproduct++ = mat1[0]*mat2[0] + mat1[3]*mat2[1] + mat1[6]*mat2[2];
*matproduct++ = mat1[1]*mat2[0] + mat1[4]*mat2[1] + mat1[7]*mat2[2];
*matproduct++ = mat1[2]*mat2[0] + mat1[5]*mat2[1] + mat1[8]*mat2[2];
*matproduct++ = mat1[0]*mat2[3] + mat1[3]*mat2[4] + mat1[6]*mat2[5];
*matproduct++ = mat1[1]*mat2[3] + mat1[4]*mat2[4] + mat1[7]*mat2[5];
*matproduct++ = mat1[2]*mat2[3] + mat1[5]*mat2[4] + mat1[8]*mat2[5];
*matproduct++ = mat1[0]*mat2[6] + mat1[3]*mat2[7] + mat1[6]*mat2[8];
*matproduct++ = mat1[1]*mat2[6] + mat1[4]*mat2[7] + mat1[7]*mat2[8];
*matproduct++ = mat1[2]*mat2[6] + mat1[5]*mat2[7] + mat1[8]*mat2[8];
}


void calcrotmat(double nx, double ny, double nz, double *rmat)

	/* Find the rotation matrix that rotates |n| radians about
		the vector given by nx,ny,nz				*/
{
double ar, ai, br, bi, hp, cp, sp;
double arar, aiai, arai2, brbr, bibi, brbi2, arbi2, aibr2, arbr2, aibi2;
double phi;

phi = sqrt(nx*nx+ny*ny+nz*nz);

if (phi == 0.0)
	{
	*rmat++ = 1;
	*rmat++	= 0;
	*rmat++ = 0;
	*rmat++ = 0;
	*rmat++ = 1;
	*rmat++	= 0;
	*rmat++ = 0;
	*rmat++ = 0;
	*rmat++ = 1;
	}

/*DEBUG_printf("calcrotmat(%6.3f,%6.3f,%6.3f) -> phi = %6.3f\n",nx,ny,nz,phi);*/

else
	{
	/* First define Cayley-Klein parameters 	*/
	hp = phi/2;		
	cp = cos(hp);
	sp = sin(hp)/phi;	/* /phi because n is unit length in defs. */
	ar = cp;
	ai = -nz*sp;
	br = ny*sp;
	bi = -nx*sp;

 	/* Make auxiliary variables to speed this up	*/

	arar = ar*ar;
	aiai = ai*ai;
	arai2 = 2*ar*ai;
	brbr = br*br;
	bibi = bi*bi;
	brbi2 = 2*br*bi;
	arbi2 = 2*ar*bi;
	aibr2 = 2*ai*br;
	arbr2 = 2*ar*br;
	aibi2 = 2*ai*bi;


	/* Make rotation matrix.	*/

	*rmat++ = arar-aiai-brbr+bibi;
	*rmat++ = -arai2-brbi2;
	*rmat++ = -arbr2+aibi2;
	*rmat++ =  arai2-brbi2; 
	*rmat++ = arar-aiai+brbr-bibi;
	*rmat++ = -aibr2-arbi2;
	*rmat++ =  arbr2+aibi2;
	*rmat++ =  arbi2-aibr2;
	*rmat++ = arar+aiai-brbr-bibi;
	}
}



void zerovec(double *vec)

/*	Set a 3x1 vector to all zeros	*/

{
*vec++=0;
*vec++=0;
*vec++=0;
}


int times2intervals( double *endtimes, double *intervals, long n)
/* ------------------------------------------------------------
	Function takes the given endtimes of intervals, and
	returns the interval lengths in an array, assuming that
	the first interval starts at 0.

	If the intervals are all greater than 0, then this
	returns 1, otherwise it returns 0.
   ------------------------------------------------------------ */

{
int allpos;
int count;
double lasttime;

allpos=1;
lasttime = 0.0;

for (count = 0; count < n; count++)
	{
	intervals[count] = endtimes[count]-lasttime;
	lasttime = endtimes[count];
	if (intervals[count] <= 0)
		allpos =0;
	}

return (allpos);
}






void blochsim(double *b1real, double *b1imag, 
        std::vector<double> sensr, std::vector<double> sensi,
		double *xgrad, double *ygrad, double *zgrad, double *tsteps, 
		int ntime, double df, 
		double dx, double dy, double dz, 
		std::vector<double> &mx, std::vector<double> &my, std::vector<double> &mz, int mode)

	/* Go through time one position.		*/

{
//int count;
int tcount;
double gammadx;
double gammady;
double gammadz;
double rotmat[9];
double amat[9], bvec[3];	/* A and B propagation matrix and vector */
//double arot[9], brot[3];	/* A and B after rotation step. */
// double decmat[9];		/* Decay matrix for each time step. */
// double decvec[3];		/* Recovery vector for each time step. */
double rotx,roty,rotz;		/* Rotation axis coordinates. */
//double mstart[3];
//double mfinish[3];
double imat[9];//, mvec[3];
double mcurr0[3];		/* Current magnetization before rotation. */
double mcurr1[3];		/* Current magnetization before decay. */

eyemat(amat); 		/* A is the identity matrix.	*/
eyemat(imat); 		/* I is the identity matrix.	*/

zerovec(bvec);
// zerovec(decvec);
// zeromat(decmat);

gammadx = dx*GAMMA;	/* Convert to Hz/cm */
gammady = dy*GAMMA;	/* Convert to Hz/cm */
gammadz = dz*GAMMA;	/* Convert to Hz/cm */


mcurr0[0] = mx[0];		/* Set starting x magnetization */
mcurr0[1] = my[0];		/* Set starting y magnetization */
mcurr0[2] = mz[0];		/* Set starting z magnetization */

int ncha = sensr.size();

for (tcount = 0; tcount < ntime; tcount++)
	{
   
	/*	Rotation 	*/

    // Calculate B1 as combination of channels
//     Two complex numbers x=a+ib and y=c+id are multiplied as follows:
// xy	=	(a+ib)(c+id)	(1)
// 	=	ac+ibc+iad-bd	(2)
// 	=	(ac-bd)+i(ad+bc).
    //   (sensr*b1r - sensi*b1i)+i(sensr*b1i+sensi*b1r)
    
    double allChaB1R =0.0, allChaB1I =0.0;
    double * currentB1r = b1real++;
    double * currentB1i = b1imag++;    
    for (int cDx = 0; cDx<ncha; cDx++)
    {   
        allChaB1R += sensr[cDx]* *currentB1r - sensi[cDx]* *currentB1i;
        allChaB1I += sensr[cDx]* *currentB1i + sensi[cDx]* *currentB1r; 
        currentB1r +=ntime;
        currentB1i +=ntime;
    }
    //DEBUG_printf("allChaB1R + 1i * allChaB1I = %f + %f i \n",allChaB1R,allChaB1I); 
	// N.B. The SENSE of ROTATION was changed in code on the B Hargreaves'
	// website. Use NEW (2013) convention here, but keep z-rotation
    // following the convention in M. Levitt. "Spin Dynamics" for the
    // (observable) -1 coherence order.
    rotz = (*xgrad++ * gammadx + *ygrad++ * gammady + *zgrad++ * gammadz +
								df*TWOPI ) * *tsteps;
//     	DEBUG_printf("rotz = %f\n",rotz); 
	rotx = ( allChaB1R * GAMMA * *tsteps);
	roty = ( allChaB1I * GAMMA * *tsteps++);
    // End of change.

	calcrotmat(rotx, roty, rotz, rotmat);

	multmatvec(rotmat,mcurr0,mcurr1);

    //copy back to mcurr0
    mcurr0[0] = mcurr1[0];
    mcurr0[1] = mcurr1[1];
    mcurr0[2] = mcurr1[2];
// 		/* 	Decay	*/
// 
// 	decvec[2]= 1- *e1;
// 	decmat[0]= *e2;
// 	decmat[4]= *e2++;
// 	decmat[8]= *e1++;
	
// 	if (mode == 1)
// 		{
// 		multmats(decmat,arot,amat);
// 		multmatvec(decmat,brot,bvec);
// 		addvecs(bvec,decvec,bvec);
// 		}
// 	else
// 		{
// 		multmatvec(decmat,mcurr1,mcurr0);
// 		addvecs(mcurr0,decvec,mcurr0);
// 		}

	
// 	DEBUG_printf("rotmat = [%6.3f  %6.3f  %6.3f ] \n",rotmat[0],rotmat[3],
// 	  			rotmat[6]);
// 	DEBUG_printf("         [%6.3f  %6.3f  %6.3f ] \n",rotmat[1],rotmat[4],
// 				rotmat[7]);
// 	DEBUG_printf("         [%6.3f  %6.3f  %6.3f ] \n",rotmat[2],rotmat[5],
// 				rotmat[8]);
//     	DEBUG_printf("mcurr0 = [%6.3f  %6.3f  %6.3f ] \n",mcurr0[0],mcurr0[1],
// 	  			mcurr0[2]);
//         	DEBUG_printf("mcurr1 = [%6.3f  %6.3f  %6.3f ] \n",mcurr1[0],mcurr1[1],
// 	  			mcurr1[2]);
    /*
	DEBUG_printf("A = [%6.3f  %6.3f  %6.3f ] \n",amat[0],amat[3],amat[6]);
	DEBUG_printf("    [%6.3f  %6.3f  %6.3f ] \n",amat[1],amat[4],amat[7]);
	DEBUG_printf("    [%6.3f  %6.3f  %6.3f ] \n",amat[2],amat[5],amat[8]);
	DEBUG_printf(" B = <%6.3f,%6.3f,%6.3f> \n",bvec[0],bvec[1],bvec[2]);
	DEBUG_printf("<mx,my,mz> = <%6.3f,%6.3f,%6.3f> \n",
		amat[6] + bvec[0], amat[7] + bvec[1], amat[8] + bvec[2]);

	DEBUG_printf("\n");
	*/

	if (mode == 2)		/* Sample output at times.  */
					/* Only do this if transient! */
		{
		mx[tcount] = mcurr1[0];
		my[tcount] = mcurr1[1];
		mz[tcount] = mcurr1[2];
	
		}	
	}



	/* If only recording the endpoint, either store the last
		point */

if (mode==0)		/* Indicates start at given m, or m0. */
	{
	mx[0] = mcurr1[0];
	my[0] = mcurr1[1];
	mz[0] = mcurr1[2];
	}


}


//blochsimfz(b1r,b1i,sensr,sensi,ncha,gx,gy,gz,tp,ntime,df,dx,dy,dz,npos,mx,my,mz,md);
void blochsimfz(double *b1real, double *b1imag,
        double *sensreal, double *sensimag,const int ncha,
        double *xgrad, double *ygrad, double *zgrad, 
		double *tsteps, 
		int ntime, double *dfreq,
		double *dxpos, double *dypos, double *dzpos, int npos, 
		double *mx, double *my, double *mz, int mode)


{
int count;
int poscount;

int ntout;

double *tstepsptr;
//double *dxptr, *dyptr, *dzptr;
//double *sensrptr, *sensiptr;
//double *dfptr;

if (mode & 2)
	ntout = ntime;
else
	ntout = 1;

tstepsptr = tsteps;


// CTR: If I wanted to make this multi-threaded, this would be the place to build a worker queue...
//     dxptr = dxpos;
//     dyptr = dypos;
//     dzptr = dzpos;
//     sensrptr = sensreal;
//     sensiptr = sensimag;
//     dfptr = dfreq;
    
#pragma omp parallel for
    for (poscount=0; poscount < npos; poscount++)

	{
        // Pull out the bits which have previously been done by pointer arithmetic
        // Sensitivities (1 timepoint x number of channels)
        std::vector<double> sensRVec(ncha), sensIVec(ncha);
        for (int cDx=0; cDx < ncha; cDx++)
        {
            sensRVec[cDx] = sensreal[poscount + (npos*cDx)];
            sensIVec[cDx] = sensimag[poscount + (npos*cDx)];
        }
        
        //df - one double
        double thisDf = dfreq[poscount];
        
        // positions (dx, dy, dz  = 3 x 1 double)
        double thisdx = dxpos[poscount];
        double thisdy = dypos[poscount];
        double thisdz = dzpos[poscount];
       
        // Mx My Mz (starting values
        std::vector<double> mxVec(ntout), myVec(ntout), mzVec(ntout);
        for (int tDx=0; tDx < ntout; tDx++)
        {
            mxVec[tDx] = mx[poscount + tDx];
            myVec[tDx] = my[poscount + tDx];
            mzVec[tDx] = mz[poscount + tDx];
        }
        
		blochsim(b1real, b1imag, sensRVec, sensIVec, xgrad, ygrad, zgrad, tsteps, ntime, 
			thisDf, thisdx, thisdy, 
			thisdz, mxVec, myVec, mzVec, mode);
	
        // Write Mx My Mz back
        for (int tDx=0; tDx < ntout; tDx++)
        {
            mx[poscount + tDx] = mxVec[tDx];
            my[poscount + tDx] = myVec[tDx];
            mz[poscount + tDx] = mzVec[tDx];
        }

    }


}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])

	/* bloch(b1,sens,gradxyz,dt,df,dxyz,mode) */
{
double *b1r;	/* Real-part of B1 field.	*/
double *b1i;	/* Imag-part of B1 field.	*/
double *sensr; // Sensitivity array (real)
double *sensi; // Sensitivity array (imaginary)
double *gx;	/* X-axis gradient. 		*/
double *gy;	/* Y-axis gradient. 		*/
double *gz;	/* Z-axis gradient. 		*/
double *tp;	/* Time steps (s)		*/
double *ti;	/* Time intervals (s) 		*/
double *df;	/* Off-resonance Frequencies (Hz)	*/
double *dx;	/* X Positions (cm)			*/
double *dy;	/* Y Positions (cm)			*/
double *dz;	/* Z Positions (cm)			*/
int md;		/* Mode - 0=from M0, 1=steady-state	*/
double tstep;	/* Time step, if single parameter */

double *mxin;	/* Input points */
double *myin;
double *mzin;

double *mxout;	/* Input points */
double *myout;
double *mzout;

double *mx;	/* Output Arrays */
double *my;
double *mz;

int gyaflag=0;	/* 1 if gy was allocated. */ 
int gzaflag=0;	/* 1 if gz was allocated. */ 
int dyaflag=0;	/* 1 if dy was allocated. */ 
int dzaflag=0;	/* 1 if dz was allocated. */ 

int ntime;	/* Number of time points. 	 */
int ntout;	/* Number of time poitns at output. */
int outsize[2];	/* Output matrix sizes		*/

int ncha;

int ngrad;	/* Number of gradient dimensions */
int nf;
int npos;	/* Number of positions.  Calculated from nposN and nposM, depends on them. */
int nposM;	/* Height of passed position matrix. */
int nposN;	/* Width of passed position matrix. */
//int nfnpos;	/* Number of frequencies * number of positions. */
int ntnpos;	/* Number of output times *number of positions. */
int count;

// #pragma omp parallel
//     mexPrintf("Hello from thread %d, nthreads %d\n", omp_get_thread_num(), omp_get_num_threads());
// 

// CTR: ERROR CHECKING. Test number of inputs.

// Special case - allow bloch('debug',true) or bloch('debug',false) syntax.
if (nrhs == 2) {
	char str1[1024];
	str1[0] = '\0';

	if (mxGetString(prhs[0],str1,sizeof(str1)-1)==0)
	{
		if (strcmp(str1,"debug") == 0)
		{
			debugflag = mxIsLogicalScalarTrue(prhs[1]);
			mexPrintf("Setting debug flag to %s.\n", debugflag ? "true" : "false" );
			return;
		}
	}
	}


// Special case - allow bloch('gamma') syntax.
if (nrhs == 1) {
	char str1[1024];
	str1[0] = '\0';

	if (mxGetString(prhs[0],str1,sizeof(str1)-1)==0)
	{
		if (strcmp(str1,"gamma") == 0)
		{
            if (nlhs < 1)
            {
                mexPrintf("GAMMA = %g.\n", GAMMA );
            }
            else
            {
                double *gamma_return;
                
                plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
                gamma_return = mxGetPr(plhs[0]);
                *gamma_return = GAMMA;
            }
			return;
		}
	}
	}

// Check number of inputs and outputs:
if (nrhs < 7) {
	mexPrintf("Hint: Type 'doc %s' for help on input parameters.\n",mexFunctionName());
    mexErrMsgIdAndTxt("bloch:BadNInput","At least 7 inputs required.");
 }

// Check all inputs are of type "double" (or char for first input):
for (count = 0; count < nrhs; count++)
{
    if (!(
            mxIsDouble(prhs[count]) // Any param can be double.
            || (count == 0 && mxIsChar(prhs[count])) // Or 1st can be char.
       ))
    {
        mexPrintf("Hint: Type 'doc %s' for help on input parameters.\n",mexFunctionName());
        mexErrMsgIdAndTxt("bloch:BadNInput","All inputs must be of type double.");
    }
}

// End CTR.

#ifdef DEBUG
  DEBUG_printf("---------------------------------------\n");
  DEBUG_printf("3D-position + frequency Bloch Simulator\n");
  DEBUG_printf("---------------------------------------\n\n");
#endif

ntime = mxGetM(prhs[0]);	/* Number of Time, RF, and Grad points */
ncha = mxGetN(prhs[0]);  // number of RF channels;
                          
/* ====================== RF (in V) =========================
 * :  If complex, split up.  If real, allocate an imaginary part. ==== */
if (mxIsComplex(prhs[0]))
	{
	b1r = mxGetPr(prhs[0]);
	b1i = mxGetPi(prhs[0]);
	}
else
	{
		mexErrMsgIdAndTxt("bloch:realB1","Rf Voltage must be complex.");

	}
#ifdef DEBUG
  DEBUG_printf("%d B1 points.\n",ntime);
#endif

  /* ====================== Sensitivities (B1 in Hz/V) =========================
 * :  If complex, split up.  If real, allocate an imaginary part. ==== */
if (mxIsComplex(prhs[1]))
	{
	sensr = mxGetPr(prhs[1]);
	sensi = mxGetPi(prhs[1]);
	}
else
	{
		mexErrMsgIdAndTxt("bloch:realSens","Sensitivities must be complex.");
	}
#ifdef DEBUG
  DEBUG_printf("%d B1 points.\n",ntime);
#endif


/* ======================= Gradients ========================= */

ngrad = mxGetM(prhs[2]) * mxGetN(prhs[2]);	/* Number of Time, RF, and Grad points */
gx = mxGetPr(prhs[2]);				/* X-gradient is first N points. */

if (ngrad < 2*ntime)		/* Need to allocate Y-gradient. */
	{
	#ifdef DEBUG
	  DEBUG_printf("Assuming 1-Dimensional Gradient\n");
	#endif
	gy = (double *)malloc(ntime * sizeof(double));
	gyaflag=1;
	for (count=0; count < ntime; count++)
		gy[count]=0.0;
	}
else
	{
	#ifdef DEBUG
	  DEBUG_printf("Assuming (at least) 2-Dimensional Gradient\n");
	#endif
	gy = gx + ntime;	/* Assign from Nx3 input array. */
	}

if (ngrad < 3*ntime)		/* Need to allocate Z-gradient. */
	{
	gz = (double *)malloc(ntime * sizeof(double));
	gzaflag=1;
	for (count=0; count < ntime; count++)
		gz[count]=0.0;
	}
else
	{
	#ifdef DEBUG
	  DEBUG_printf("Assuming 3-Dimensional Gradient\n");
	#endif
	gz = gx + 2*ntime; 	/* Assign from Nx3 input array. */
	}

	/* Warning if Gradient length is not 1x, 2x, or 3x RF length. */

	
#ifdef DEBUG
  DEBUG_printf("%d Gradient Points (total) \n",ngrad);
#endif
if ( (ngrad != ntime) && (ngrad != 2*ntime) && (ngrad != 3*ntime) )
		mexErrMsgIdAndTxt("bloch:BadGradientLength","Gradient length differs from B1 length.");


if (gx == NULL) 
	mexErrMsgIdAndTxt("bloch:BadGradientLength","gx is not allocated.");
if (gy == NULL) 
	mexErrMsgIdAndTxt("bloch:BadGradientLength","gy is not allocated.");
if (gz == NULL) 
	mexErrMsgIdAndTxt("bloch:BadGradientLength","gz is not allocated.");



/* === Time points ===== */

/*	THREE Cases:
		1) Single value given -> this is the interval length for all.
		2) List of intervals given.
		3) Monotonically INCREASING list of end times given.

	For all cases, the goal is for tp to have the intervals.
*/

ti = NULL;
tp = mxGetPr(prhs[3]);

if (mxGetM(prhs[3]) * mxGetN(prhs[3]) == 1)	/* === Case 1 === */
	{
	tp = (double *)malloc(ntime * sizeof(double));
	tstep = *(mxGetPr(prhs[3]));
	for (count =0; count < ntime; count++)
		tp[count]=tstep;
    }

else if (mxGetM(prhs[3]) * mxGetN(prhs[3]) != ntime)
	mexErrMsgIdAndTxt("bloch:BadB1Length","Time-point length differs from B1 length.");

else	
	{
	tp = mxGetPr(prhs[3]);
	ti = (double *)malloc(ntime * sizeof(double));
	if (( times2intervals( tp, ti, ntime )))
		{
		DEBUG_printf("Times are monotonically increasing. \n");
		tp = ti;
		}
	}

/* === Frequency Points ===== */

df = mxGetPr(prhs[4]);
nf = mxGetM(prhs[4]) * mxGetN(prhs[4]);
	
#ifdef DEBUG
  DEBUG_printf("%d Frequency points.\n",nf);
#endif


/* === Position Points ===== */

nposM = mxGetM(prhs[5]);
nposN = mxGetN(prhs[5]);

#ifdef DEBUG
  DEBUG_printf("Position vector is %d x %d. \n",nposM,nposN);
#endif

if (nposN==3)			/* Assume 3 position dimensions given */
	{
	npos = nposM;
	#ifdef DEBUG
	  DEBUG_printf("Assuming %d 3-Dimensional Positions\n",npos);
	#endif
	dx = mxGetPr(prhs[5]);
	dy = dx + npos;
	dz = dy + npos;
	}

else if (nposN==2)		/* Assume only 2 position dimensions given */
	{
	npos = nposM;
	#ifdef DEBUG
	  DEBUG_printf("Assuming %d 2-Dimensional Positions\n",npos);
	#endif
	dx = mxGetPr(prhs[5]);
	dy = dx + npos;
	dz = (double *)malloc(npos * sizeof(double));
	dzaflag=1;
	for (count=0; count < npos; count++)
		dz[count]=0.0;
	}

else				/* Either 1xN, Nx1 or something random.  In all these
				   cases we assume that 1 position is given, because it
				   is too much work to try to figure out anything else! */
	{
	npos = nposM * nposN;
	#ifdef DEBUG
	  DEBUG_printf("Assuming %d 1-Dimensional Positions\n",npos);
	#endif
	dx = mxGetPr(prhs[5]);
	dy = (double *)malloc(npos * sizeof(double));
	dz = (double *)malloc(npos * sizeof(double));
	dyaflag=1;
	dzaflag=1;
	for (count=0; count < npos; count++)
		{
		dy[count]=0.0;
		dz[count]=0.0;
		}
	#ifdef DEBUG
	  if ((nposM !=1) && (nposN!=1))		
		{
		DEBUG_printf("Position vector should be 1xN, Nx1, Nx2 or Nx3. \n");
		DEBUG_printf(" -> Assuming 1 position dimension is given. \n");
		}	
	#endif
	}

if (dx == NULL) 
	mexErrMsgIdAndTxt("bloch:BadPosition","dx is not allocated.");
if (dy == NULL) 
	mexErrMsgIdAndTxt("bloch:BadPosition","dy is not allocated.");
if (dz == NULL) 
	mexErrMsgIdAndTxt("bloch:BadPosition","dz is not allocated.");

if (npos != nf)
{
   mexErrMsgIdAndTxt("bloch:BadNumPos","npos must equal df size");
}
  

/* ===== Mode, defaults to 0 (simulate single endpoint, transient). ==== */

if (nrhs > 6)
	md = (int)(*mxGetPr(prhs[6]));		
else
	md = 0;


if (md & 2)
	ntout = ntime;		/* Include time points.	*/
else
	ntout = 1;

#ifdef DEBUG
  DEBUG_printf("Mode = %d, %d Output Time Points \n",md,ntout);
#endif

ntnpos = ntout*npos;


#ifdef DEBUG
if ((md & 1)==0)
	DEBUG_printf("Simulation from Initial Condition.\n");
else
	DEBUG_printf("Simulation of Steady-State.\n");


if ((md & 2)==0)
	DEBUG_printf("Simulation to Endpoint. \n");
else
	DEBUG_printf("Simulation over Time.\n");
#endif


/* ===== Allocate Output Magnetization vectors arrays.	*/

plhs[0] = mxCreateDoubleMatrix(ntnpos,1,mxREAL);	/* Mx, output. */
plhs[1] = mxCreateDoubleMatrix(ntnpos,1,mxREAL);	/* My, output. */
plhs[2] = mxCreateDoubleMatrix(ntnpos,1,mxREAL);	/* Mz, output. */

mx = mxGetPr(plhs[0]);
my = mxGetPr(plhs[1]);
mz = mxGetPr(plhs[2]);

mxout = mx;
myout = my;
mzout = mz;

// Initial magnetisation
	#ifdef DEBUG

	  DEBUG_printf(" --> Using [0; 0; 1] for initial magnetization. \n");
	#endif
	for (count =0; count < npos; count++)
		{
		*mxout = 0;	/* Set magnetization to Equilibrium */
		*myout = 0;
		*mzout = 1;
		mxout += ntout;
		myout += ntout;
		mzout += ntout;
		}
		

/* ======= Do The Simulation! ====== */

#ifdef DEBUG
  DEBUG_printf("Calling blochsimfz() function in Mex file.\n");
//     DEBUG_printf("ntime = %i.\n",ntime);
//     DEBUG_printf("tp[0] = %f, [1] = %f.\n",tp[0],tp[1]);

#endif

blochsimfz(b1r,b1i,sensr,sensi,ncha,gx,gy,gz,tp,ntime,df,dx,dy,dz,npos,mx,my,mz,md);

/* ======= Reshape Output Matrices ====== */

// if ((ntout > 1) && (nf > 1) && (npos > 1))
// 	{
// 	outsize[0]=ntout;
// 	outsize[1]=npos;
// 	mxSetDimensions(plhs[0],outsize,2);  /* Set to 3D array. */
// 	mxSetDimensions(plhs[1],outsize,2);  /* Set to 3D array. */
// 	mxSetDimensions(plhs[2],outsize,2);  /* Set to 3D array. */
// 	}
// else			/* Basically "squeeze" the matrix. */
// 	{
	if (ntout > 1)	
		{
		outsize[0]=ntout;
		outsize[1]=npos;
		}
	else
		{
		outsize[0]=npos;
		outsize[1]=1;
		}
	mxSetDimensions(plhs[0],outsize,2);  /* Set to 2D array. */
	mxSetDimensions(plhs[1],outsize,2);  /* Set to 2D array. */
	mxSetDimensions(plhs[2],outsize,2);  /* Set to 2D array. */
// 	}


/* ====== Free up allocated memory, if necessary. ===== */

if (!mxIsComplex(prhs[0]))
	free(b1i);	/* We had to allocate this before. */

if (mxGetM(prhs[2]) * mxGetN(prhs[2]) == 1)
	free(tp);	/* We had to allocate this. */

if (ti != NULL)
	free(ti);

if (dyaflag==1)
	free(dy);
if (dzaflag==1)
	free(dz);
if (gyaflag==1)
	free(gy);
if (gzaflag==1)
	free(gz);

}






