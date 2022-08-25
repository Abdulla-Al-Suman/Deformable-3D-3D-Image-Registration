#include "mex.h"
#include "matrix.h"

/*
 * CalcEPD.c
 * Calculates the similarity measure D = I.*E0
 *
 * This is a MEX-file for MATLAB.
 */



// void CalcEPD3D(double *EL, int M_EL, int N_EL, double *DR, double *S, int Dims)
// {
// 	mwSize i;
//     //mwSize j;
//     double a, b, c;
// 	double sumD;
// 	sumD = 0;
//     a=EL[0]+EL[1];
// 	printf("EL[0] = %d\n", a);	
// // 	for (i=0; i<M_EL; i++) {
// //         //for (j=0; j<N_EL; j++){
// //          a=EL[i]-1; // A[m + M*n] 0  m  M?1 and 0  n  N ? 1)
// //          b=EL[i+M_EL]-1;
// //          c=EL[i+M_EL*2]-1;
// // 
// //   		//sumD += DR[a + b * M_DR + c * (M_DR * N_DR)];
// //        // }
// // 	}
// // 	
// // 	*S = sumD/M_EL*3;
// }

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
		int nrhs, const mxArray *prhs[])
{
	/* variable declarations here */
	
	double *EL;
	double *DR;
    double W = 0;
    int i, a, b, c; 
    mwSize N, MDR, NDR;
    const int *pDims;
// 	double *S;
//     double d;
//     mwSize M_EL;
//     mwSize N_EL;
//     mwSize Dims;
//     mwSize M_DR;
//     mwSize N_DR;	
	
	/* code here */
	
	/* check for proper number of arguments */
	if(nrhs!=2) {
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs",
				"Two inputs required.");
	}
	if(nlhs!=1) {
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs",
				"One output required.");
	}
    
    /* Check data type of first input argument */
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) {
        mexErrMsgTxt("Input argument must be a real double.");
    }
		
	/* create pointers to the real data in the input matrices  */	
    EL = mxGetPr(prhs[0]);
    DR = mxGetPr(prhs[1]);
    pDims = mxGetDimensions(prhs[1]);
    MDR = pDims[0];
    NDR = pDims[1];
  
    // do the calculations
    N =(int) mxGetM(prhs[0]);
    for (i=0; i<N; i++){
        a= *(EL+i)-1;//row,column // A[m + M*n] 0  m  M?1 and 0  n  N ? 1)
        b= *(EL+N+i)-1;//column, row
        c= *(EL+(2*N)+i)-1;//slice number
        // W += a+b+c;
        W += *(DR+b+(a*MDR)+(c*(MDR*NDR))); 
    }
    W = W/(N*3);
    
    //W = *(DR+65536) * *(DR+65536);
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxGetPr(plhs[0])[0] = W;
	
	/* get dimensions of the input matrix */
// 	M_EL = mxGetM(prhs[0]);
//     N_EL = mxGetN(prhs[0]);
//     Dims = mxGetDimensions(prhs[1]);
//     M_DR = mxGetM(prhs[1]);
//     N_DR = mxGetN(prhs[1]);
	
	/* create the output variables */
// 	plhs[0] = mxCreateDoubleScalar(0);
// 	
// 	/* get a pointer to the real data in the output variables */
// 	S = mxGetPr(plhs[0]);
	
	/* call the computational routine */
	// CalcEPD3D(EL,M_EL,N_EL,DR,S,Dims); 
	
}