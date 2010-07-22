#ifndef lint
static char sccsid[]="@(#)cg_glm_get_Beta_ResSS.c  1.01 Christian Gaser 07/09/09";
#endif

#include <math.h>
#include "mex.h"

#include "spm_mapping.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int n_subj, n_beta, n_slices, n_values;
  int n, i, j, k, z, ind1, ind2, n_slices_x_values;
  double *Beta, *ResSS, *image, *estimates, *mask, *X, *pKX, *TH, *W;
  double sum, W2, ival;
  MAPTYPE *maps, *map_mask, *get_maps();
  static double mat[] = {1, 0, 0, 0, 0, 1, 0, 0,  0, 0, 1, 0, 0, 0, 0, 1};

  if (nrhs != 6) mexErrMsgTxt("Six input arguments required.");
  if (nlhs != 2) mexErrMsgTxt("Two output arguments required.");
 
  n_subj = mxGetM(prhs[2]);
  n_beta = mxGetN(prhs[2]);

  map_mask = get_maps(prhs[1], &n);
  if (n!=1)
  {
    free_maps(map_mask, n);
    mexErrMsgTxt("Only single file as mask allowed.");
  }

  maps = get_maps(prhs[0], &n);
  if (n!=n_subj)
  {
    free_maps(maps, n);
    free_maps(map_mask, 1);
    mexErrMsgTxt("Different number of scans in design matrix.");
  }

  for(i=1; i<n_subj; i++)
  {
    if (  maps[i].dim[0] != maps[0].dim[0] ||
      maps[i].dim[1] != maps[0].dim[1] ||
      maps[i].dim[2] != maps[0].dim[2])
      {
        free_maps(maps, n_subj);
        mexErrMsgTxt("Incompatible image dimensions.");
      }
  }

  n_slices = maps[0].dim[2];
  n_values = maps[0].dim[0]*maps[0].dim[1];
  
  if ((n_slices!=map_mask[0].dim[2]) || (n_values!=map_mask[0].dim[0]*map_mask[0].dim[1]))
  {
    free_maps(maps, n);
    free_maps(map_mask, 1);
    mexErrMsgTxt("Incompatible dimensions between mask and images.");
  }
  
  X = mxGetPr(prhs[2]);
  pKX = mxGetPr(prhs[3]);
  TH = mxGetPr(prhs[4]);
  W = mxGetPr(prhs[5]);
  
  n_subj = mxGetM(prhs[2]);
  if (n_subj!=mxGetM(prhs[4]))
  {
    free_maps(maps, n);
    free_maps(map_mask, 1);
    mexErrMsgTxt("Incompatible dimensions of thresholds.");
  }

  for(i=0; i<n_subj; i++)
    if (mxIsInf(TH[i]))
      TH[i] = -1e15;
  
  plhs[0] = mxCreateDoubleMatrix(n_slices*n_values, n_beta, mxREAL);
  Beta = mxGetPr(plhs[0]);
  plhs[1] = mxCreateDoubleMatrix(n_slices*n_values, 1, mxREAL);
  ResSS = mxGetPr(plhs[1]);

  image  = (double *)mxCalloc(n_values, sizeof(double));
  estimates  = (double *)mxCalloc(n_values*n_subj, sizeof(double));
  mask = (double *)mxCalloc(n_values, sizeof(double));
  
  n_slices_x_values = n_slices*n_values;
  
  /* initialize Beta and ResSS with zeros */
  for(z=0; z<n_slices; z++)
  {
    ind2 = z*n_values;  

    for(j=0; j<n_values; j++)
    {
      ResSS[j + ind2] = 0;
      for(k=0; k<n_beta; k++)
        Beta[j + ind2 + (k*n_slices_x_values)] = 0;
    }
  }


  for(z=0; z<n_slices; z++)
  {
    mat[14] = z + 1.0;  
    ind2 = z*n_values;  

    /* load mask */  
    slice(mat, mask, map_mask[0].dim[0],map_mask[0].dim[1], &map_mask[0], 0, 0.0);
    
    for(i=0; i<n_subj; i++)
    {
      W2 = W[i];
      slice(mat, image, maps[i].dim[0],maps[i].dim[1], &maps[i], 0, 0.0);
      for(j=0; j<n_values; j++)
      {
        ival = W2*image[j];
        if ((mask[j] > 0) & (ival>TH[i]))
        {
          /* initialize estimates with image values */
          estimates[j + (i*n_values)] = ival;
          /* calculate betas */
          for(k=0; k<n_beta; k++)
            Beta[j + ind2 + (k*n_slices_x_values)] += pKX[k + (i*n_beta)] * ival;
        }
      }
    }

    /* get estimates */
    for(i=0; i<n_subj; i++)
    {
      ind1 = i*n_values;
      for(j=0; j<n_values; j++)
      {
        ival = W2*image[j];
        if ((mask[j] > 0) & (ival>TH[i]))
        {
          sum = 0.0;
          /* calculate difference between estimates and original values */
          for(k=0; k<n_beta; k++)
            sum += (X[i + (k*n_subj)] * Beta[j + ind2 + (k*n_slices_x_values)]);
          estimates[j + ind1] -= sum;
        }
      }
    }

    /* calculate sum of residual squares */
    for(j=0; j<n_values; j++)
    {
      if (mask[j] > 0)
      {
        sum = 0.0;
        for(i=0; i<n_subj; i++)
        {
          ind1 = j + (i*n_values);  
          sum += estimates[ind1] * estimates[ind1];
        }
        ResSS[j + ind2] = sum;
      }
    }
  }
  
  mxFree((char *)mask);
  mxFree((char *)image);
  mxFree((char *)estimates);
  free_maps(maps, n_subj);
  free_maps(map_mask, 1);

}
