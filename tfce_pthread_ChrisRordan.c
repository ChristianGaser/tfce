#include <stdio.h>
#include <immintrin.h> //optional: cross-platform memory alignment
#include <string.h> //optional, memset

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
 * $Id: tfceMex_pthread.c 182 2019-08-09 15:02:06Z gaser $ 
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
    const  mwSize *dims;    
    double E;
    double H;  
    int    calc_neg;  
    int    threaded;  
} myargument;

pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

/*
  clustering based on BrainModelVolumeTFCE.cxx from caret
*/
#define ThreadCR
#ifdef ThreadCR

static int nifti_tfce(double * inimg, double * outimg, int nx, int ny, int nz, double H, double E, int c)  { 
	//https://www.fmrib.ox.ac.uk/datasets/techrep/tr08ss1/tr08ss1.pdf
	int nvox3D = nx * ny * nz;
	int numk = c;
	int i,j,x,y,z;
	int dx,dy,ii,jj,step,qlo,qhi;
	double valToAdd;

	if ((c != 6) && (c != 18) && (c != 26)) {
		printf("suitable values for c are 6, 18 or 26\n");
		numk = 6;	
	}
	//set up kernel to search for neighbors. Since we already included sides, we do not worry about A<->P and L<->R wrap
	int32_t * k = (int32_t *)_mm_malloc(3*numk*sizeof(int32_t), 64); //kernel: offset, x, y
	int mxDx = 1; //connectivity 6: faces only
	if (numk == 18) mxDx = 2; //connectivity 18: faces+edges
	if (numk == 26) mxDx = 3; //connectivity 26: faces+edges+corners
	j = 0;
	for (z = -1; z <= 1; z++ )
		for (y = -1; y <= 1; y++ )
			for (x = -1; x <= 1; x++ ) {
				int dx = abs(x)+abs(y)+abs(z);
				if ((dx > mxDx) || (dx == 0)) continue; 
				k[j] = x + (y * nx) + (z * nx * ny);
				k[j+numk] = x; //avoid left-right wrap
				k[j+numk+numk] = x; //avoid anterior-posterior wrap				
				j++;
			} //for x
	//identify clusters
	double mx = (inimg[0]);
	for (ii = 0; ii < nvox3D; ii++ )
		mx = MAX((inimg[ii]),mx);	
	double dh = mx/100.0;
	for (i = 0; i < nvox3D; i++ )
		outimg[i] = 0.0;	
	int32_t * q = (int32_t *)_mm_malloc(nvox3D*sizeof(int32_t), 64); //queue with untested seed
	uint8_t * vxs = (uint8_t *)_mm_malloc(nvox3D*sizeof(uint8_t), 64); 
	int n_steps = (int)ceil(mx/dh);
	for (step=0; step<n_steps; step++) {
		double thresh = (step+1)*dh;
		memset(vxs, 0, nvox3D*sizeof(uint8_t));
		for (i = 0; i < nvox3D; i++ )
			if (inimg[i] >= thresh)
				vxs[i] = 1; //survives, unclustered
		i = 0;
		for (z = 0; z < nz; z++ )
			for (y = 0; y < ny; y++ )
				for (x = 0; x < nx; x++ ) {
					if (vxs[i] == 0) {
						i++;
						continue;
					} //voxel did not survive or already clustered
					qlo = 0;
					qhi = 0; 
					q[qhi] = i; //add starting voxel as seed in queue
					vxs[i] = 0; //do not find again!
					while (qhi >= qlo) { //first in, first out queue
						//retire one seed, add 0..6, 0..18 or 0..26 new ones (depending on connectivity)
						for (j = 0; j < numk; j++) {
							jj = q[qlo] + k[j];
							if ((jj < 0) || (jj >= nvox3D)) continue; //voxel in volume
							if (vxs[jj] == 0) continue; //already found or did not survive threshold
							dx = x+k[j+numk];
							if ((dx < 0) || (dx >= nx)) continue; //wrapped left-right
							dy = y+k[j+numk+numk];
							if ((dy < 0) || (dy >= ny)) continue; //wrapped anterior-posterior
							//add new seed:
							vxs[jj] = 0; //do not find again!
							qhi++;
							q[qhi] = jj;	 
						}        	
						qlo++;
					} //while qhi >= qlo: continue until all seeds tested
					valToAdd = pow(qhi+1, E) * pow(thresh, H); //"supporting section", Dark Gray in Figure 1
					for (j = 0; j <= qhi; j++)
						outimg[q[j]] += valToAdd;
					i++;
				} //for each voxel
	} //for each step
	_mm_free (q);
	_mm_free (vxs);
	_mm_free (k);	
	return 0;
} //nifti_tfce()


void* ThreadFuncCR( void *pArguments ) {
	myargument arg;
	arg=*(myargument *) pArguments;
	double *inData   = arg.inData;    
	double *outData  = arg.outData;
	double thresh   = arg.thresh;    
	double E        = arg.E;    
	double H        = arg.H; 
	const mwSize *dims     = arg.dims; 
	int c = 26; //to do: set
	nifti_tfce(inData,outData, (long)dims[0], (long)dims[1], (long)dims[2], H, E, c);
	//to do: calc_neg
	return NULL;
}
#endif

void* ThreadFunc( void *pArguments )
{
  double valToAdd;
  int i, j, k, ti, tj, tk, maxi, maxj, maxk, mini, minj, mink, growingInd, growingCur;
  long ind, ind1;
  int calc_neg;
  long numVoxels;
  char *flagUsed;
  short *growing;
  double *inData, *outData;
  double thresh, E, H;    
  const mwSize *dims;    

  myargument arg;
  arg=*(myargument *) pArguments;

  inData   = arg.inData;    
  outData  = arg.outData;    
  thresh   = arg.thresh;    
  dims     = arg.dims;    
  E        = arg.E;    
  H        = arg.H;    
  calc_neg = arg.calc_neg;
  
  numVoxels = (long)dims[0] * (long)dims[1] * (long)dims[2];

  flagUsed = (char*) malloc(numVoxels*sizeof(char));
  growing  = (short*)malloc(numVoxels*3*sizeof(short));

  for (i = 0; i < numVoxels; ++i) flagUsed[i] = 0;

  for (k = 0; k < (int)dims[2]; ++k) for (j = 0; j < (int)dims[1]; ++j) for (i = 0; i < (int)dims[0]; ++i)
  {
    ind = (long)k*((long)dims[0]*(long)dims[1]) + ((long)j*(long)dims[0]) + (long)i;
            
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
          ind1 = (long)tk*((long)dims[0]*(long)dims[1]) + ((long)tj*(long)dims[0]) + (long)ti;
          
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
        if (arg.threaded) pthread_mutex_lock(&mutex);
        outData[growing[growingCur + 2]*(dims[0]*dims[1])+(growing[growingCur + 1]*dims[0])+growing[growingCur]] += valToAdd;
        growingCur += 3;
        if (arg.threaded) pthread_mutex_unlock(&mutex);
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
          ind1 = (long)tk*((long)dims[0]*(long)dims[1]) + ((long)tj*(long)dims[0]) + (long)ti;
          
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
        if (arg.threaded) pthread_mutex_lock(&mutex);
        outData[growing[growingCur + 2]*(dims[0]*dims[1])+(growing[growingCur + 1]*dims[0])+growing[growingCur]] -= valToAdd;
        growingCur += 3;
        if (arg.threaded) pthread_mutex_unlock(&mutex);
      }
    }    
  }
   
  free(flagUsed);
  free(growing);
  
  if (arg.threaded) pthread_exit((void*) 0);

  return NULL;
}

void tfce(double *inData, double *outData, double dh, const mwSize *dims, double E, double H, int calc_neg)
{
  double fmax = 0.0, curThr, tmp_value;
  int i, n_steps;
  long numVoxels = (long)dims[0] * (long)dims[1] * (long)dims[2];
  myargument *ThreadArgs;  
  pthread_t *ThreadList; 

  for (i = 0; i < numVoxels; ++i)
  {
     tmp_value = fabs(inData[i]);
     if (tmp_value > fmax) fmax = tmp_value;
     outData[i] = 0.0;
  }
   
  /* get # of steps = # of threads */
  n_steps = (int)ceil(fmax/dh);

  /* Reserve room for handles of threads in ThreadList*/
  ThreadList = (pthread_t *) calloc(n_steps,sizeof(pthread_t));
  ThreadArgs = (myargument*) calloc( n_steps,sizeof(myargument));
  if (pthread_mutex_init(&mutex, NULL) != 0)
	{
			printf("\n mutex init failed\n");
			exit(1);
	}

  for (i=0; i<n_steps; i++)
  {         
    curThr = (double)(i+1)*dh;
    
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
      
  for (i=0; i<n_steps; i++)
  {         
    if(pthread_create(&(ThreadList[i]), NULL, ThreadFunc, &ThreadArgs[i]))
    {
       printf("Threads cannot be created\n");
       exit(1);
    }        
  }

  for (i=0; i<n_steps; i++)
    pthread_join(ThreadList[i],NULL);
  
      
  pthread_mutex_destroy(&mutex);

  free(ThreadList);
  free(ThreadArgs);

}

void tfce_singlethreaded(double *inData, double *outData, double dh, const mwSize *dims, double E, double H, int calc_neg)
{
  double fmax = 0.0, tmp_value;
  int i, n_steps;
  long numVoxels = (long)dims[0] * (long)dims[1] * (long)dims[2];
  myargument arg;
  arg.inData = inData;
  arg.outData = outData;
  arg.dims = dims;   
  arg.E = E;
  arg.H = H;  
  arg.calc_neg = calc_neg;  
  arg.threaded = 0;
  #ifdef ThreadCR
  ThreadFuncCR(&arg);
  #else

  for (i = 0; i < numVoxels; ++i)
  {
     tmp_value = fabs(inData[i]);
     if (tmp_value > fmax) fmax = tmp_value;
     outData[i] = 0.0;
  }
   
  /* get # of steps = # of threads */
  n_steps = (int)ceil(fmax/dh);
    
  
  for (i=0; i<n_steps; i++)
  {         
    arg.thresh = (double)(i+1)*dh;
    ThreadFunc(&arg);
  }
  #endif
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

/* Declarations */
double *inData, *outData, dh, E, H;
int ndim, calc_neg;
const mwSize *dims;

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



