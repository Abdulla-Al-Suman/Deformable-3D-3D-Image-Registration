#include "mex.h"
#include "matrix.h"
#include "math.h"

//#define INF (unsigned)!((int)0)

/*
 * This is a MEX-file for MATLAB.
 */

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    /* variable declarations here */
    
    int *ER;
    double *DR, FA[14], FA2[11], FA3[10], FA4[8], FA5[5], FA6[8], FA7[6],FA8[9], FA9[7], FA10[6], BA[14], BA2[11], BA3[10], BA4[8], minFA, minBA, *IV;
    double BA5[5], BA6[8], BA7[6], BA8[9], BA9[7], BA10[6];
    mwSize M, N, O;
    const int *pDims;
    double d1 = 3;  //1;
    double d2 = 4;       //1.314;
    double d3 = 5;       //1.628;
    int i, ii, jj, j, k, nDimNum, c;
    //double pinf=mxGetInf;
    // mxArray   *Data;
    
    
    
    /* code here */
    
    /* check for proper number of arguments */
    if(nrhs!=2) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs",
                "Two input required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs",
                "One output required.");
    }
    
    /* create pointers to the real data in the input matrices  */
    
    ER = mxGetPr(prhs[0]);
    IV = mxGetPr(prhs[1]);
    pDims = mxGetDimensions(prhs[0]);
    nDimNum = mxGetNumberOfDimensions(prhs[0]);
    M =pDims[0];
    N =pDims[1];
    O =pDims[2];
    //mexPrintf("3rd dim=%d\n", O);
    
// 	/* create the output matrix */
    plhs[0] = mxCreateNumericArray(nDimNum, pDims, mxDOUBLE_CLASS, mxREAL);
    DR =  mxGetPr( plhs[0]);
    
    for (k=0; k<O; k++){
        for (i=0; i<M; i++){
            for (j=0; j<N; j++){
                if (*(ER+i+j*M+k*(M*N))==0){
                    *(DR+i+j*M+k*(M*N))=*IV;
                    //*(DR+i+j*M+k*(M*N))=1.0/0.0;
                    //*(DR+i+j*M+k*(M*N))= double mxGetInf(void);
                }
                else if (*(ER+i+j*M+k*(M*N))==1) {
                    *(DR+i+j*M+k*(M*N))=0;
                }
            }
            
        }
        
    }
    
    // Forward
    for (k=1; k<O; k++){     
        for (i=1; i<M; i++){
            for (j=1; j<N; j++){
                if (i<M-1 && j<N-1){
                FA[0]=*(DR+i+j*M+k*(M*N));
                FA[1]=*(DR+i+(j-1)*M+k*(M*N))+d1;
                FA[2]=*(DR+i-1+(j-1)*M+k*(M*N))+d2;
                FA[3]=*(DR+i-1+j*M+k*(M*N))+d1;
                FA[4]=*(DR+i-1+(j+1)*M+k*(M*N))+d2;
                FA[5]=*(DR+i-1+(j-1)*M+(k-1)*(M*N))+d3;
                FA[6]=*(DR+i-1+j*M+(k-1)*(M*N))+d2;
                FA[7]=*(DR+i-1+(j+1)*M+(k-1)*(M*N))+d3;
                FA[8]=*(DR+i+(j-1)*M+(k-1)*(M*N))+d2;
                FA[9]=*(DR+i+j*M+(k-1)*(M*N))+d1;
                FA[10]=*(DR+i+(j+1)*M+(k-1)*(M*N))+d2;
                FA[11]=*(DR+i+1+(j-1)*M+(k-1)*(M*N))+d3;
                FA[12]=*(DR+i+1+j*M+(k-1)*(M*N))+d2;
                FA[13]=*(DR+i+1+(j+1)*M+(k-1)*(M*N))+d3;
                
                minFA = FA[0];
                for ( c = 1 ; c < 14 ; c++ )
                {
                    if ( FA[c] < minFA )
                    {
                        minFA= FA[c];
                    }
                }
                *(DR+i+j*M+k*(M*N))=minFA;               
               }
                
                
                else if (k==O-1) {  // for last slice
                    for (ii=0; ii<M; ii++){
                        for (jj=0; jj<N; jj++){
                            
                            if (ii==0 && jj==0){
                                FA5[0]=*(DR+ii+jj*M+k*(M*N));
                                FA5[1]=*(DR+ii+jj*M+(k-1)*(M*N))+d1;
                                FA5[2]=*(DR+ii+(jj+1)*M+(k-1)*(M*N))+d2;
                                FA5[3]=*(DR+ii+1+jj*M+(k-1)*(M*N))+d2;
                                FA5[4]=*(DR+ii+1+(jj+1)*M+(k-1)*(M*N))+d3;
                                
                                minFA = FA5[0];
                                for ( c = 1 ; c < 5 ; c++ )
                                {
                                    if ( FA5[c] < minFA )
                                    {
                                        minFA= FA5[c];
                                    }
                                }
                                *(DR+ii+jj*M+k*(M*N))=minFA;
                            }
                            else if (ii==0 && 0<jj<N-1) {
                                FA6[0]=*(DR+ii+jj*M+k*(M*N));
                                FA6[1]=*(DR+ii+(jj-1)*M+k*(M*N))+d1;
                                FA6[2]=*(DR+ii+(jj-1)*M+(k-1)*(M*N))+d2;
                                FA6[3]=*(DR+ii+jj*M+(k-1)*(M*N))+d1;
                                FA6[4]=*(DR+ii+(jj+1)*M+(k-1)*(M*N))+d2;
                                FA6[5]=*(DR+ii+1+(jj-1)*M+(k-1)*(M*N))+d3;
                                FA6[6]=*(DR+ii+1+jj*M+(k-1)*(M*N))+d2;
                                FA6[7]=*(DR+ii+1+(jj+1)*M+(k-1)*(M*N))+d3;
                                
                                minFA = FA6[0];
                                for ( c = 1 ; c < 8 ; c++ )
                                {
                                    if ( FA6[c] < minFA )
                                    {
                                        minFA= FA6[c];
                                    }
                                }
                                *(DR+ii+jj*M+k*(M*N))=minFA;  
                            }
                            else if (ii==0 && jj==N-1) {
                                FA7[0]=*(DR+ii+jj*M+k*(M*N));
                                FA7[1]=*(DR+ii+(jj-1)*M+k*(M*N))+d1;
                                FA7[2]=*(DR+ii+(jj-1)*M+(k-1)*(M*N))+d2;
                                FA7[3]=*(DR+ii+jj*M+(k-1)*(M*N))+d1;
                                FA7[4]=*(DR+ii+1+(jj-1)*M+(k-1)*(M*N))+d3;
                                FA7[5]=*(DR+ii+1+jj*M+(k-1)*(M*N))+d2;
                                
                                minFA = FA7[0];
                                for ( c = 1 ; c < 6 ; c++ )
                                {
                                    if ( FA7[c] < minFA )
                                    {
                                        minFA= FA7[c];
                                    }
                                }
                                  *(DR+ii+jj*M+k*(M*N))=minFA;                               
                            }
                            else if (0<ii<M-1 && jj==0) {
                                FA8[0]=*(DR+ii+jj*M+k*(M*N));
                                FA8[1]=*(DR+ii-1+jj*M+k*(M*N))+d1;
                                FA8[2]=*(DR+ii-1+(jj+1)*M+k*(M*N))+d2;
                                FA8[3]=*(DR+ii-1+jj*M+(k-1)*(M*N))+d2;
                                FA8[4]=*(DR+ii-1+(jj+1)*M+(k-1)*(M*N))+d3;
                                FA8[5]=*(DR+ii+jj*M+(k-1)*(M*N))+d1;
                                FA8[6]=*(DR+ii+(jj+1)*M+(k-1)*(M*N))+d2;
                                FA8[7]=*(DR+ii+1+jj*M+(k-1)*(M*N))+d2;
                                FA8[8]=*(DR+ii+1+(jj+1)*M+(k-1)*(M*N))+d3;
                                
                                minFA = FA8[0];
                                for ( c = 1 ; c < 9 ; c++ )
                                {
                                    if ( FA8[c] < minFA )
                                    {
                                        minFA= FA8[c];
                                    }
                                }
                                *(DR+ii+jj*M+k*(M*N))=minFA;                            
                            }
                            else if (ii==M-1 && jj==0) {
                                FA9[0]=*(DR+ii+jj*M+k*(M*N));
                                FA9[1]=*(DR+ii-1+jj*M+k*(M*N))+d1;
                                FA9[2]=*(DR+ii-1+(jj+1)*M+k*(M*N))+d2;
                                FA9[3]=*(DR+ii-1+jj*M+(k-1)*(M*N))+d2;
                                FA9[4]=*(DR+ii-1+(jj+1)*M+(k-1)*(M*N))+d3;
                                FA9[5]=*(DR+ii+jj*M+(k-1)*(M*N))+d1;
                                FA9[6]=*(DR+ii+(jj+1)*M+(k-1)*(M*N))+d2;
                                
                                minFA = FA9[0];
                                for ( c = 1 ; c < 7 ; c++ )
                                {
                                    if ( FA9[c] < minFA )
                                    {
                                        minFA= FA9[c];
                                    }
                                }
                                *(DR+ii+jj*M+k*(M*N))=minFA;                               
                            }
                            
                            
                        }
                    }                   
                }

                
            }           
        }      
    }
    
    

    
    // Backward
    for (k=O-2; k>-1; k--){
        for (i=M-2; i>-1; i--){
            for (j=N-2; j>-1; j--){
                if (i>0 && j>0){
                    BA[0]=*(DR+i+j*M+k*(M*N));
                    BA[1]=*(DR+i+(j+1)*M+k*(M*N))+d1;
                    BA[2]=*(DR+i+1+(j+1)*M+k*(M*N))+d2;
                    BA[3]=*(DR+i+1+j*M+k*(M*N))+d1;
                    BA[4]=*(DR+i+1+(j-1)*M+k*(M*N))+d2;
                    BA[5]=*(DR+i-1+(j-1)*M+(k+1)*(M*N))+d3;
                    BA[6]=*(DR+i-1+j*M+(k+1)*(M*N))+d2;
                    BA[7]=*(DR+i-1+(j+1)*M+(k+1)*(M*N))+d3;
                    BA[8]=*(DR+i+(j-1)*M+(k+1)*(M*N))+d2;
                    BA[9]=*(DR+i+j*M+(k+1)*(M*N))+d1;
                    BA[10]=*(DR+i+(j+1)*M+(k+1)*(M*N))+d2;
                    BA[11]=*(DR+i+1+(j-1)*M+(k+1)*(M*N))+d3;
                    BA[12]=*(DR+i+1+j*M+(k+1)*(M*N))+d2;
                    BA[13]=*(DR+i+1+(j+1)*M+(k+1)*(M*N))+d3;
                    
                    minBA = BA[0];
                    for ( c = 1 ; c < 14 ; c++ )
                    {
                        if ( BA[c] < minBA )
                        {
                            minBA= BA[c];
                        }
                    }
                    *(DR+i+j*M+k*(M*N))=minBA;
                }
                
                
                else if (k==0) {  // for First slice
                    for (ii=M-1; ii>-1; ii--){
                        for (jj=N-1; jj>-1; jj--){
                            if (ii==M-1 && jj==N-1){
                                BA5[0]=*(DR+ii+jj*M+k*(M*N));
                                BA5[1]=*(DR+ii-1+(jj-1)*M+(k+1)*(M*N))+d3;
                                BA5[2]=*(DR+ii-1+jj*M+(k+1)*(M*N))+d2;
                                BA5[3]=*(DR+ii+(jj-1)*M+(k+1)*(M*N))+d2;
                                BA5[4]=*(DR+ii+jj*M+(k+1)*(M*N))+d1;
                                
                                minBA = BA5[0];
                                for ( c = 1 ; c < 5 ; c++ )
                                {
                                    if ( BA5[c] < minBA )
                                    {
                                        minBA= BA5[c];
                                    }
                                }
                                *(DR+ii+jj*M+k*(M*N))=minBA;
                            }
                            else if (ii==M-1 && N-1>jj>0){
                                BA6[0]=*(DR+ii+jj*M+k*(M*N));
                                BA6[1]=*(DR+ii+(jj+1)*M+k*(M*N))+d1;
                                BA6[2]=*(DR+ii-1+(jj-1)*M+(k+1)*(M*N))+d3;
                                BA6[3]=*(DR+ii-1+jj*M+(k+1)*(M*N))+d2;
                                BA6[4]=*(DR+ii-1+(jj+1)*M+(k+1)*(M*N))+d3;
                                BA6[5]=*(DR+ii+(jj-1)*M+(k+1)*(M*N))+d2;
                                BA6[6]=*(DR+ii+jj*M+(k+1)*(M*N))+d1;
                                BA6[7]=*(DR+ii+(jj+1)*M+(k+1)*(M*N))+d2;
                                
                                minBA = BA6[0];
                                for ( c = 1 ; c < 8 ; c++ )
                                {
                                    if ( BA6[c] < minBA )
                                    {
                                        minBA= BA6[c];
                                    }
                                }
                                *(DR+ii+jj*M+k*(M*N))=minBA;
                            }
                            else if (ii==M-1 && jj==0){
                                BA7[0]=*(DR+ii+jj*M+k*(M*N));
                                BA7[1]=*(DR+ii+(jj+1)*M+k*(M*N))+d1;
                                BA7[2]=*(DR+ii-1+jj*M+(k+1)*(M*N))+d2;
                                BA7[3]=*(DR+ii-1+(jj+1)*M+(k+1)*(M*N))+d3;
                                BA7[4]=*(DR+ii+jj*M+(k+1)*(M*N))+d1;
                                BA7[5]=*(DR+ii+(jj+1)*M+(k+1)*(M*N))+d2;
                                
                                minBA = BA7[0];
                                for ( c = 1 ; c < 6 ; c++ )
                                {
                                    if ( BA7[c] < minBA )
                                    {
                                        minBA= BA7[c];
                                    }
                                }
                                *(DR+ii+jj*M+k*(M*N))=minBA;
                            }
                            else if (M-1>ii>0 && jj==N-1){
                                BA8[0]=*(DR+ii+jj*M+k*(M*N));
                                BA8[1]=*(DR+ii+1+jj*M+k*(M*N))+d1;
                                BA8[2]=*(DR+ii+1+(jj-1)*M+k*(M*N))+d2;
                                BA8[3]=*(DR+ii-1+(jj-1)*M+(k+1)*(M*N))+d3;
                                BA8[4]=*(DR+ii-1+jj*M+(k+1)*(M*N))+d2;
                                BA8[5]=*(DR+ii+(jj-1)*M+(k+1)*(M*N))+d2;
                                BA8[6]=*(DR+ii+jj*M+(k+1)*(M*N))+d1;
                                BA8[7]=*(DR+ii+1+(jj-1)*M+(k+1)*(M*N))+d3;
                                BA8[8]=*(DR+ii+1+jj*M+(k+1)*(M*N))+d2;
                                
                                minBA = BA8[0];
                                for ( c = 1 ; c < 9 ; c++ )
                                {
                                    if ( BA8[c] < minBA )
                                    {
                                        minBA= BA8[c];
                                    }
                                }
                                *(DR+ii+jj*M+k*(M*N))=minBA;
                            }
                            else if (ii==0 && jj==N-1){
                                BA9[0]=*(DR+ii+jj*M+k*(M*N));
                                BA9[1]=*(DR+ii+1+jj*M+k*(M*N))+d1;
                                BA9[2]=*(DR+ii+1+(jj-1)*M+k*(M*N))+d2;
                                BA9[3]=*(DR+ii+(jj-1)*M+(k+1)*(M*N))+d2;
                                BA9[4]=*(DR+ii+jj*M+(k+1)*(M*N))+d1;
                                BA9[5]=*(DR+ii+1+(jj-1)*M+(k+1)*(M*N))+d3;
                                BA9[6]=*(DR+ii+1+jj*M+(k+1)*(M*N))+d2;
                                
                                minBA = BA9[0];
                                for ( c = 1 ; c < 7 ; c++ )
                                {
                                    if ( BA9[c] < minBA )
                                    {
                                        minBA= BA9[c];
                                    }
                                }
                                *(DR+ii+jj*M+k*(M*N))=minBA;
                            }
                            
                        }
                    }
                }
                
            }
        }      
    }
    
 
    // 2nd Forward
    for (k=1; k<O; k++){
        for (i=1; i<M; i++){
            for (j=1; j<N; j++){
                
                
                if (i==M-1 && j<N-1) {
                    FA2[0]=*(DR+i+j*M+k*(M*N));
                    FA2[1]=*(DR+i+(j-1)*M+k*(M*N))+d1;
                    FA2[2]=*(DR+i-1+(j-1)*M+k*(M*N))+d2;
                    FA2[3]=*(DR+i-1+j*M+k*(M*N))+d1;
                    FA2[4]=*(DR+i-1+(j+1)*M+k*(M*N))+d2;
                    FA2[5]=*(DR+i-1+(j-1)*M+(k-1)*(M*N))+d3;
                    FA2[6]=*(DR+i-1+j*M+(k-1)*(M*N))+d2;
                    FA2[7]=*(DR+i-1+(j+1)*M+(k-1)*(M*N))+d3;
                    FA2[8]=*(DR+i+(j-1)*M+(k-1)*(M*N))+d2;
                    FA2[9]=*(DR+i+j*M+(k-1)*(M*N))+d1;
                    FA2[10]=*(DR+i+(j+1)*M+(k-1)*(M*N))+d2;
                    
                    minFA = FA2[0];
                    for ( c = 1 ; c < 11 ; c++ )
                    {
                        if ( FA2[c] < minFA )
                        {
                            minFA= FA2[c];
                        }
                    }
                    *(DR+i+j*M+k*(M*N))=minFA;
                    
                }
                else if (i<M-1 && j==N-1) {
                    FA3[0]=*(DR+i+j*M+k*(M*N));
                    FA3[1]=*(DR+i+(j-1)*M+k*(M*N))+d1;
                    FA3[2]=*(DR+i-1+(j-1)*M+k*(M*N))+d2;
                    FA3[3]=*(DR+i-1+j*M+k*(M*N))+d1;
                    FA3[4]=*(DR+i-1+(j-1)*M+(k-1)*(M*N))+d3;
                    FA3[5]=*(DR+i-1+j*M+(k-1)*(M*N))+d2;
                    FA3[6]=*(DR+i+(j-1)*M+(k-1)*(M*N))+d2;
                    FA3[7]=*(DR+i+j*M+(k-1)*(M*N))+d1;
                    FA3[8]=*(DR+i+1+(j-1)*M+(k-1)*(M*N))+d3;
                    FA3[9]=*(DR+i+1+j*M+(k-1)*(M*N))+d2;
                    
                    
                    minFA = FA3[0];
                    for ( c = 1 ; c < 10 ; c++ )
                    {
                        if ( FA3[c] < minFA )
                        {
                            minFA= FA3[c];
                        }
                    }
                    *(DR+i+j*M+k*(M*N))=minFA;
                    
                }
                else if (i==M-1 && j==N-1) {
                    FA4[0]=*(DR+i+j*M+k*(M*N));
                    FA4[1]=*(DR+i+(j-1)*M+k*(M*N))+d1;
                    FA4[2]=*(DR+i-1+(j-1)*M+k*(M*N))+d2;
                    FA4[3]=*(DR+i-1+j*M+k*(M*N))+d1;
                    FA4[4]=*(DR+i-1+(j-1)*M+(k-1)*(M*N))+d3;
                    FA4[5]=*(DR+i-1+j*M+(k-1)*(M*N))+d2;
                    FA4[6]=*(DR+i+(j-1)*M+(k-1)*(M*N))+d2;
                    FA4[7]=*(DR+i+j*M+(k-1)*(M*N))+d1;
                    
                    minFA = FA4[0];
                    for ( c = 1 ; c < 8 ; c++ )
                    {
                        if ( FA4[c] < minFA )
                        {
                            minFA= FA4[c];
                        }
                    }
                    *(DR+i+j*M+k*(M*N))=minFA;
                    
                }
                
            }
        }
        // for top right corner voxels
        FA10[0]=*(DR+0+(N-1)*M+k*(M*N));
        FA10[1]=*(DR+0+(N-2)*M+k*(M*N))+d1;
        FA10[2]=*(DR+0+(N-1)*M+(k-1)*(M*N))+d1;
        FA10[3]=*(DR+0+(N-2)*M+(k-1)*(M*N))+d2;
        FA10[4]=*(DR+0+1+(N-2)*M+(k-1)*(M*N))+d3;
        FA10[5]=*(DR+0+1+(N-1)*M+(k-1)*(M*N))+d2;
        
        minFA = FA10[0];
        for ( c = 1 ; c < 6 ; c++ )
        {
            if ( FA10[c] < minFA )
            {
                minFA= FA10[c];
            }
        }
        *(DR+0+(N-1)*M+k*(M*N))=minFA;
        
    }
    
    // 2nd Backward
    
    for (k=O-2; k>-1; k--){
        for (i=M-2; i>-1; i--){
            for (j=N-2; j>-1; j--){
                
                if (i==0 && j>0) {
                    BA2[0]=*(DR+i+j*M+k*(M*N));
                    BA2[1]=*(DR+i+(j+1)*M+k*(M*N))+d1;
                    BA2[2]=*(DR+i+1+(j+1)*M+k*(M*N))+d2;
                    BA2[3]=*(DR+i+1+j*M+k*(M*N))+d1;
                    BA2[4]=*(DR+i+1+(j-1)*M+k*(M*N))+d2;
                    BA2[5]=*(DR+i+(j-1)*M+(k+1)*(M*N))+d2;
                    BA2[6]=*(DR+i+j*M+(k+1)*(M*N))+d1;
                    BA2[7]=*(DR+i+(j+1)*M+(k+1)*(M*N))+d2;
                    BA2[8]=*(DR+i+1+(j-1)*M+(k+1)*(M*N))+d3;
                    BA2[9]=*(DR+i+1+j*M+(k+1)*(M*N))+d2;
                    BA2[10]=*(DR+i+1+(j+1)*M+(k+1)*(M*N))+d3;
                    
                    minBA = BA2[0];
                    for ( c = 1 ; c < 11 ; c++ )
                    {
                        if ( BA2[c] < minBA )
                        {
                            minBA= BA2[c];
                        }
                    }
                    *(DR+i+j*M+k*(M*N))=minBA;
                    
                }
                else if (i>0 && j==0) {
                    BA3[0]=*(DR+i+j*M+k*(M*N));
                    BA3[1]=*(DR+i+(j+1)*M+k*(M*N))+d1;
                    BA3[2]=*(DR+i+1+(j+1)*M+k*(M*N))+d2;
                    BA3[3]=*(DR+i+1+j*M+k*(M*N))+d1;
                    BA3[4]=*(DR+i-1+j*M+(k+1)*(M*N))+d2;
                    BA3[5]=*(DR+i-1+(j+1)*M+(k+1)*(M*N))+d3;
                    BA3[6]=*(DR+i+j*M+(k+1)*(M*N))+d1;
                    BA3[7]=*(DR+i+(j+1)*M+(k+1)*(M*N))+d2;
                    BA3[8]=*(DR+i+1+j*M+(k+1)*(M*N))+d2;
                    BA3[9]=*(DR+i+1+(j+1)*M+(k+1)*(M*N))+d3;
                    
                    minBA = BA3[0];
                    for ( c = 1 ; c < 10 ; c++ )
                    {
                        if ( BA3[c] < minBA )
                        {
                            minBA= BA3[c];
                        }
                    }
                    *(DR+i+j*M+k*(M*N))=minBA;
                    
                }
                else if (i==0 && j==0) {
                    BA4[0]=*(DR+i+j*M+k*(M*N));
                    BA4[1]=*(DR+i+(j+1)*M+k*(M*N))+d1;
                    BA4[2]=*(DR+i+1+(j+1)*M+k*(M*N))+d2;
                    BA4[3]=*(DR+i+1+j*M+k*(M*N))+d1;
                    BA4[4]=*(DR+i+j*M+(k+1)*(M*N))+d1;
                    BA4[5]=*(DR+i+(j+1)*M+(k+1)*(M*N))+d2;
                    BA4[6]=*(DR+i+1+j*M+(k+1)*(M*N))+d2;
                    BA4[7]=*(DR+i+1+(j+1)*M+(k+1)*(M*N))+d3;
                    
                    minBA = BA4[0];
                    for ( c = 1 ; c < 8 ; c++ )
                    {
                        if ( BA4[c] < minBA )
                        {
                            minBA= BA4[c];
                        }
                    }
                    *(DR+i+j*M+k*(M*N))=minBA;
                }
                
            }
        }
        
        // for bottom left corner voxels
        BA10[0]=*(DR+(M-1)+(0)*M+k*(M*N));
        BA10[1]=*(DR+(M-1)+(1)*M+k*(M*N))+d1;
        BA10[2]=*(DR+(M-1)+(0)*M+(k+1)*(M*N))+d1;
        BA10[3]=*(DR+(M-1)+(1)*M+(k+1)*(M*N))+d2;
        BA10[4]=*(DR+(M-2)+(0)*M+(k+1)*(M*N))+d2;
        BA10[5]=*(DR+(M-2)+(1)*M+(k+1)*(M*N))+d3;
        
        minBA = BA10[0];
        for ( c = 1 ; c < 6 ; c++ )
        {
            if ( BA10[c] < minBA )
            {
                minBA= BA10[c];
            }
        }
        *(DR+(M-1)+(0)*M+k*(M*N))=minBA;
    }    
    
}