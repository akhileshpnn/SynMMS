/*==========================================================
 * maltesCA.cpp - Cellular Automaton for simulation of 
 *                reaction-diffusion systems
 *
 * The calling syntax is:
 *
 *      outMatrix = 
 *              maltesCA(initImg,configFile,outFile,iters)
 *
 * This is a MEX-file for MATLAB.
 *
 *========================================================*/
/* $Revision: 0.9 $ */

#define _MATLAB_MEX

#include "mex.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>

#include "cellField_StatTub.h"

using namespace std;

//cellField knows the imageData variable via 'extern',
//DO NOT CHANGE THIS VARIABLE NAME!
unsigned char *imageData;               /* image char input matrix */
double *initVals;
mwSize ncols,nrows,nstacks;
mxArray *Params;
double *outMatrix;
double *outMean;
int Reps;
int OutFreq;

const int nSpecies=9;

/* The computational routine */
void maltesCA()
{
/*    for (int x=0; x<nrows; x++)
        for (int y=0; y<ncols; y++) {
            outMatrix[x+nrows*y]=imageData[x+nrows*y];
        }
    // Initialize from passed image-data-matrix
*/
    int FieldSize[3]={(int)ncols,(int)nrows,(int)nstacks};
    int cellVol=ncols*nrows*nstacks*nSpecies;
    cellField<double> TestField(FieldSize);

    TestField.showField(outMatrix);
    for (int i=1;i<=Reps;i++) {
        for (int j=0; j<10; ++j)
            TestField.interact();
        TestField.diffuse_varGeom();
        if (i%OutFreq==0) {
            TestField.showField(&(outMatrix[(i/OutFreq)*cellVol]));
        }
    }
//    TestField.returnMean(outMean);
//  second left-hand variable defunct atm
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    
    /* check for proper number of arguments */
    if( nrhs<1 || ! mxIsStruct(prhs[0]) ) {
        const char *fields[]={"StochVal","DiffRad","Reps","Stacks","OutFreq","IntPar"};
        plhs[0] = mxCreateStructMatrix(1,1,6,fields);
        mxSetField(plhs[0], 0, "StochVal",mxCreateDoubleScalar(0.0));
        mxSetField(plhs[0], 0, "DiffRad",mxCreateDoubleScalar(1.0));
        mxSetField(plhs[0], 0, "Reps",mxCreateDoubleScalar(100));
        mxSetField(plhs[0], 0, "Stacks",mxCreateDoubleScalar(1));
        mxSetField(plhs[0], 0, "OutFreq",mxCreateDoubleScalar(1));
        mxSetField(plhs[0], 0, "IntPar",mxCreateDoubleScalar(0));
        mexPrintf("parameter-struct like this as first argument needed\n");
        return;
    }
    if(nrhs<3 ) {
        mexErrMsgIdAndTxt("MyToolbox:maltesCA:nrhs","More arguments needed: \
                           params first, compartmental geometry second, initMatrix third.");
    }
    if(nlhs!=2) {
        mexErrMsgIdAndTxt("MyToolbox:maltesCA:nlhs","Two output (5D-stack and Mean-Value) required.");
    }
    Params=mxDuplicateArray(prhs[0]);
    Reps=(int)(*(static_cast<double*>(mxGetData(mxGetField(Params,0,"Reps")))));
    if (! Reps) Reps=1;
    nstacks=(mwSize)(*(static_cast<double*>(mxGetData(mxGetField(Params,0,"Stacks")))));
    if (!nstacks) nstacks=1;
    
//    imageData = (char*)(mxGetPr(prhs[0]));
    imageData = reinterpret_cast<unsigned char*>(mxGetChars(prhs[1]));
    initVals = reinterpret_cast<double*>(mxGetChars(prhs[2]));


    /* create a pointer to the real data in the input matrix
    inMatrix = mxGetPr(prhs[1]);
     */

    /* get dimensions of the input matrix */
    ncols = mxGetM(prhs[1]);
    nrows = mxGetN(prhs[1])/nstacks;

    /* create the output matrix */
//    mwSize Dims[3]={ncols,nrows,Reps+1};
//    plhs[0] = mxCreateNumericArray(3,Dims,mxDOUBLE_CLASS,mxREAL);
    
    OutFreq=(int)(*(static_cast<double*>(mxGetData(mxGetField(Params,0,"OutFreq")))));
    mwSize Dim5=(mwSize)(Reps/OutFreq)+1;

    mwSize Dims[5]={ncols,nrows,nstacks,nSpecies,Dim5};
    plhs[0] = mxCreateNumericArray(5,Dims,mxDOUBLE_CLASS,mxREAL);
    
/*
    mwSize Dims_2[1]={Reps/OutFreq+1};
    plhs[1] = mxCreateNumericArray(1,Dims_2,mxDOUBLE_CLASS,mxREAL);
*/    
    mwSize Dims_2[1]={Dim5};
    plhs[1] = mxCreateNumericArray(1,Dims_2,mxDOUBLE_CLASS,mxREAL);
    /* get a pointer to the real data in the output matrix */
    outMatrix = mxGetPr(plhs[0]);
    outMean = mxGetPr(plhs[1]);

    /* call the computational routine */
    maltesCA();
}
