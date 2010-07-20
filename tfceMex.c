/*
 * Christian Gaser
 * $Id: tfceMex.c 226 2009-12-03 07:22:53Z gaser $ 
 *
 */

#include "math.h"
#include "mex.h"
#include <stdlib.h>

#ifndef MAX
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#endif

#ifndef MIN
#define MIN(A,B) ((A) > (B) ? (B) : (A))
#endif

void tfce(double *inData, double *outData, int numSteps, const int *dims)
{
   double fmax = 0.0f, thresh, valToAdd, delta;
   double E = 0.5, H = 2.0;
   int i, j, k, ti, tj, tk, maxi, maxj, maxk, mini, minj, mink, temp, growingInd, growingCur;
   int numVoxels = dims[0] * dims[1] * dims[2];
   bool* flagUsed;
   short* growing;
   
   flagUsed = (bool*)malloc(numVoxels*sizeof(bool));
   growing  = (short*)malloc(numVoxels*3*sizeof(short));
   
   for (temp = 0; temp < numVoxels; ++temp)
   {
      if (inData[temp] > fmax) fmax = inData[temp];
      outData[temp] = 0.0f;
   }
   
   delta = fmax/numSteps;
   for (thresh = delta / 2.0; thresh < fmax; thresh += delta)
   {
      for (temp = 0; temp < numVoxels; ++temp) flagUsed[temp] = false;
      for (k = 0; k < dims[2]; ++k)
      {
         for (j = 0; j < dims[1]; ++j)
         {
            for (i = 0; i < dims[0]; ++i)
            {
               temp = k*(dims[0]*dims[1])+(j*dims[0])+i;
               if (!flagUsed[temp] && inData[temp] >= thresh)
               {
                  flagUsed[temp] = true;
                  growingInd = 3;
                  growingCur = 0;
                  growing[0] = i;
                  growing[1] = j;
                  growing[2] = k;
                  while (growingCur < growingInd)
                  {
                     maxi = MIN(dims[0], growing[growingCur    ] + 2);
                     maxj = MIN(dims[1], growing[growingCur + 1] + 2);
                     maxk = MIN(dims[2], growing[growingCur + 2] + 2);
                     mini = MAX(0, growing[growingCur    ] - 1);
                     minj = MAX(0, growing[growingCur + 1] - 1);
                     mink = MAX(0, growing[growingCur + 2] - 1);
                     for (tk = mink; tk < maxk; ++tk)
                     {
                        for (tj = minj; tj < maxj; ++tj)
                        {
                           for (ti = mini; ti < maxi; ++ti)
                           {
                              temp = tk*(dims[0]*dims[1])+(tj*dims[0])+ti;
                              if (!flagUsed[temp] && inData[temp] >= thresh)
                              {
                                 flagUsed[temp] = true;
                                 growing[growingInd    ] = ti;
                                 growing[growingInd + 1] = tj;
                                 growing[growingInd + 2] = tk;
                                 growingInd += 3;
                              }
                           }
                        }
                     }
                     growingCur += 3;
                  }
                  growingCur = 0;
                  valToAdd = pow(growingInd / 3.0f, E) * pow(thresh, H) * delta;

                  while (growingCur < growingInd)
                  {
                     outData[growing[growingCur + 2]*(dims[0]*dims[1])+(growing[growingCur + 1]*dims[0])+growing[growingCur]] += valToAdd;
                     growingCur += 3;
                  }
               }
            }
         }
      }
   }
   
   free(flagUsed);
   free(growing);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

/* Declarations */
double *inData, *outData;
int numSteps, ndim;
const int *dims;

/* check inputs */
if (nrhs!=2)
  mexErrMsgTxt("2 inputs required.");
else if (nlhs>2)
  mexErrMsgTxt("Too many output arguments.");
  
if (!mxIsDouble(prhs[0]))
	mexErrMsgTxt("First argument must be double.");

/* get input inDatage */
inData = (double*)mxGetPr(prhs[0]);

ndim = mxGetNumberOfDimensions(prhs[0]);
if (ndim!=3)
  mxErrMsgTxt("Images does not have 3 dimensions.");
  
dims = mxGetDimensions(prhs[0]);

/* get parameters */
numSteps = (int)(mxGetScalar(prhs[1]));

/*Allocate memory and assign output pointer*/
plhs[0] = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS, mxREAL);

/*Get a pointer to the data space in our newly allocated memory*/
outData = mxGetPr(plhs[0]);

tfce(inData, outData, numSteps, dims); 

return;
}

