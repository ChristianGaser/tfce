/*
 * TFCE estimation
 * t         - T map 
 * dh        - steps size (e.g. dh = max(abs(t))/100)
 * E         - TFCE parameter
 * H         - TFCE parameter
 * calc_neg  - also calc neg. TFCE values (default 1)
 * single_threaded - use single thread only (default 0)
 *
 * Christian Gaser
 * $Id: tfceMex_pthread.c 78 2015-09-10 10:23:01Z gaser $ 
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

/* Multithreading stuff */
#include <pthread.h>

typedef struct{
    double *inData;   
    double *outData;
    double thresh;    
    const int *dims;    
    double E;
    double H;  
    int calc_neg;  
    int threaded;  
} myargument;

/*
  clustering based on BrainModelVolumeTFCE.cxx from caret
*/

void* ThreadFunc( void* pArguments )
{
  double valToAdd;
  int i, j, k, ti, tj, tk, maxi, maxj, maxk, mini, minj, mink, ind1, growingInd, growingCur;
  long ind;
  int numVoxels, calc_neg;
  char *flagUsed;
  short *growing;
  double *inData, *outData;
  double thresh, E, H;    
  const int *dims;    
  pthread_mutex_t mutex;
   
  myargument arg;
  arg=*(myargument *) pArguments;

  inData   = arg.inData;    
  outData  = arg.outData;    
  thresh   = arg.thresh;    
  dims     = arg.dims;    
  E        = arg.E;    
  H        = arg.H;    
  calc_neg = arg.calc_neg;
  
  numVoxels = dims[0] * dims[1] * dims[2];

  if (arg.threaded)
    pthread_mutex_lock(&mutex);

  flagUsed = (char*) malloc(numVoxels*sizeof(char));
  growing  = (short*)malloc(numVoxels*3*sizeof(short));

  for (i = 0; i < numVoxels; ++i) flagUsed[i] = 0;

  for (k = 0; k < (int)dims[2]; ++k) for (j = 0; j < (int)dims[1]; ++j) for (i = 0; i < (int)dims[0]; ++i)
  {
    ind = (long)k*(dims[0]*dims[1]) + (j*dims[0]) + i;
            
    /* estimate positive tfce values */ 
    if (!flagUsed[ind] && (inData[ind] >= thresh))
    {
      flagUsed[ind] = 1;
      growingInd = 3;
      growingCur = 0;
      growing[0] = i;
      growing[1] = j;
      growing[2] = k;
      
      while (growingCur < growingInd)
      {
        maxi = MIN((int)dims[0], growing[growingCur    ] + 2);
        maxj = MIN((int)dims[1], growing[growingCur + 1] + 2);
        maxk = MIN((int)dims[2], growing[growingCur + 2] + 2);
        
        mini = MAX(0, growing[growingCur    ] - 1);
        minj = MAX(0, growing[growingCur + 1] - 1);
        mink = MAX(0, growing[growingCur + 2] - 1);
        
        for (tk = mink; tk < maxk; ++tk) for (tj = minj; tj < maxj; ++tj) for (ti = mini; ti < maxi; ++ti)
        {
          ind1 = tk*(dims[0]*dims[1]) + (tj*dims[0]) + ti;
          
          if (!flagUsed[ind1] && inData[ind1] >= thresh)
          {
            flagUsed[ind1] = 1;
            growing[growingInd    ] = ti;
            growing[growingInd + 1] = tj;
            growing[growingInd + 2] = tk;
            growingInd += 3;
          }
        }
        growingCur += 3;
      }
      
      growingCur = 0;
      valToAdd = pow(growingInd / 3.0, E) * pow(thresh, H);

      while (growingCur < growingInd)
      {
        outData[growing[growingCur + 2]*(dims[0]*dims[1])+(growing[growingCur + 1]*dims[0])+growing[growingCur]] += valToAdd;
        growingCur += 3;
      }
    }    

    /* estimate negative tfce values */ 
    if (!flagUsed[ind] && (-inData[ind] >= thresh) && calc_neg)
    {
      flagUsed[ind] = 1;
      growingInd = 3;
      growingCur = 0;
      growing[0] = i;
      growing[1] = j;
      growing[2] = k;
      
      while (growingCur < growingInd)
      {
        maxi = MIN((int)dims[0], growing[growingCur    ] + 2);
        maxj = MIN((int)dims[1], growing[growingCur + 1] + 2);
        maxk = MIN((int)dims[2], growing[growingCur + 2] + 2);
        
        mini = MAX(0, growing[growingCur    ] - 1);
        minj = MAX(0, growing[growingCur + 1] - 1);
        mink = MAX(0, growing[growingCur + 2] - 1);
        
        for (tk = mink; tk < maxk; ++tk) for (tj = minj; tj < maxj; ++tj) for (ti = mini; ti < maxi; ++ti)
        {
          ind1 = tk*(dims[0]*dims[1]) + (tj*dims[0]) + ti;
          
          if (!flagUsed[ind1] && -inData[ind1] >= thresh)
          {
            flagUsed[ind1] = 1;
            growing[growingInd    ] = ti;
            growing[growingInd + 1] = tj;
            growing[growingInd + 2] = tk;
            growingInd += 3;
          }
        }
        growingCur += 3;
      }
      
      growingCur = 0;
      valToAdd = pow(growingInd / 3.0, E) * pow(thresh, H);

      while (growingCur < growingInd)
      {
        outData[growing[growingCur + 2]*(dims[0]*dims[1])+(growing[growingCur + 1]*dims[0])+growing[growingCur]] -= valToAdd;
        growingCur += 3;
      }
    }    
  }
   
  free(flagUsed);
  free(growing);
  
  if (arg.threaded) {
    pthread_mutex_unlock(&mutex);
    pthread_exit(0);    
  }
  
}

void tfce(double *inData, double *outData, double dh, const int *dims, double E, double H, int calc_neg)
{
  double fmax = 0.0, curThr, tmp_value;
  int i, Nthreads;
  int numVoxels = dims[0] * dims[1] * dims[2];
  pthread_mutex_t mutex;   
  myargument *ThreadArgs;  

  pthread_t * ThreadList;

  for (i = 0; i < numVoxels; ++i)
  {
     tmp_value = fabs(inData[i]);
     if (tmp_value > fmax) fmax = tmp_value;
     outData[i] = 0.0;
  }
   
  /* get # of steps = # of threads */
  Nthreads = (int)ceil(fmax/dh);
  
  /* Reserve room for handles of threads in ThreadList*/
  ThreadList = (pthread_t *) calloc(Nthreads,sizeof(pthread_t));
  ThreadArgs = (myargument*) calloc( Nthreads,sizeof(myargument));
  pthread_mutex_init(&mutex, NULL);
  
  for (i=0; i<Nthreads; i++)
  {         
    curThr = (i+1)*dh;
    
	/* Make Thread Structure   */
    ThreadArgs[i].inData = inData;
    ThreadArgs[i].outData = outData;
    ThreadArgs[i].thresh = curThr;
    ThreadArgs[i].dims = dims;   
    ThreadArgs[i].E = E;
    ThreadArgs[i].H = H;  
    ThreadArgs[i].calc_neg = calc_neg;  
    ThreadArgs[i].threaded = 1;  
      
  }

  for (i=0; i<Nthreads; i++)
  {
    if(pthread_create(&ThreadList[i], NULL, ThreadFunc,&ThreadArgs[i]))
    {
       printf("Threads cannot be created\n");
       exit(1);
    }        
  }
  
  for (i=0; i<Nthreads; i++)
    pthread_join(ThreadList[i],NULL);

  pthread_mutex_destroy(&mutex);

  free(ThreadList);
  free(ThreadArgs);

}

void tfce_singlethreaded(double *inData, double *outData, double dh, const int *dims, double E, double H, int calc_neg)
{
  double fmax = 0.0, tmp_value;
  int i, n_steps;
  int numVoxels = dims[0] * dims[1] * dims[2];
  myargument arg;
  
  for (i = 0; i < numVoxels; ++i)
  {
     tmp_value = fabs(inData[i]);
     if (tmp_value > fmax) fmax = tmp_value;
     outData[i] = 0.0;
  }
   
  /* get # of steps = # of threads */
  n_steps = (int)ceil(fmax/dh);
    
  arg.inData = inData;
  arg.outData = outData;
  arg.dims = dims;   
  arg.E = E;
  arg.H = H;  
  arg.calc_neg = calc_neg;  
  arg.threaded = 0;
  
  for (i=0; i<n_steps; i++)
  {         
    arg.thresh = (i+1)*dh;
    ThreadFunc(&arg);
  }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

/* Declarations */
double *inData, *outData, dh, E, H;
int ndim, calc_neg;
const int *dims;

/* check inputs */
if (nrhs<4)
  mexErrMsgTxt("4 inputs required.");
else if (nlhs>2)
  mexErrMsgTxt("Too many output arguments.");
  
if (!mxIsDouble(prhs[0]))
	mexErrMsgTxt("First argument must be double.");

/* get input */
inData = (double*)mxGetPr(prhs[0]);

ndim = mxGetNumberOfDimensions(prhs[0]);
if (ndim!=3)
  mexErrMsgTxt("Images does not have 3 dimensions.");
  
dims = mxGetDimensions(prhs[0]);

/* get parameters */
dh = (double)(mxGetScalar(prhs[1]));
E = (double)(mxGetScalar(prhs[2]));
H = (double)(mxGetScalar(prhs[3]));

if (nrhs>4)
     calc_neg = (int)(mxGetScalar(prhs[4]));
else calc_neg = 1;

/* Allocate memory and assign output pointer */
plhs[0] = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS, mxREAL);

/* Get a pointer to the data space in our newly allocated memory */
outData = mxGetPr(plhs[0]);

/* use single-threaded version on request */
if (nrhs>5)
  if ((int)mxGetScalar(prhs[5]))
    tfce_singlethreaded(inData, outData, dh, dims, E, H, calc_neg); 
  else tfce(inData, outData, dh, dims, E, H, calc_neg);
else tfce(inData, outData, dh, dims, E, H, calc_neg); 

return;
}

