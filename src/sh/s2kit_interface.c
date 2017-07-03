/*
Short wrapper function to fascilitate spherical harmonics transformation and inverse SH trafos
*/

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "fftw3.h"
#include "s2kit10/makeweights.h"
#include "s2kit10/cospmls.h"
#include "s2kit10/FST_semi_memo.h"
#include "s2kit10/csecond.h"

/**************************************************/
/**************************************************/

double *seminaive_naive_tablespace;
double *trans_seminaive_naive_tablespace;
double **seminaive_naive_table;
double **trans_seminaive_naive_table;

double *weights;
double *workspace;

int oldbw = 0;

void release_data();

void initialize(int bw)
{
	int cutoff = bw ;
	int size = 2*bw;

	if(oldbw != bw || weights == NULL || workspace == NULL)
	{
		if(oldbw != 0)
			release_data();

		oldbw = bw;
		//printf("Initializing data for spherical harmonics!\n");
		weights = (double *) malloc(sizeof(double) * 4 * bw);
		workspace = (double *) malloc(sizeof(double) * ((12 * (bw*bw)) + (16 * bw)));

		seminaive_naive_tablespace = (double *) malloc(sizeof(double) * (Reduced_Naive_TableSize(bw,cutoff) + Reduced_SpharmonicTableSize(bw,cutoff)));
		seminaive_naive_table= SemiNaive_Naive_Pml_Table(bw, cutoff, seminaive_naive_tablespace, workspace);

		trans_seminaive_naive_tablespace =	(double *) malloc(sizeof(double) * (Reduced_Naive_TableSize(bw,cutoff) + Reduced_SpharmonicTableSize(bw,cutoff)));
		trans_seminaive_naive_table = Transpose_SemiNaive_Naive_Pml_Table(seminaive_naive_table, bw, cutoff,trans_seminaive_naive_tablespace, workspace);

		/* now make the weights */
		makeweights( bw, weights );
	}
}

void release_data()
{
	//printf("Releasing old data structures...\n");
	 /* clean up */
	free(workspace);
	free(trans_seminaive_naive_table);
	free(trans_seminaive_naive_tablespace);
	free(seminaive_naive_table);
	free(seminaive_naive_tablespace);

	free(weights);
}

void transform_sample(double *real_data, double *imag_data,
					  double *real_output, double *imag_output,
					  int bw, int forward)
{

  FILE *fp ;
  int i;
  int l, m, dummy;
  int  order ;
  int rank, howmany_rank ;

  fftw_plan dctPlan, fftPlan ;
  fftw_iodim dims[1], howmany_dims[1];

  /*** ASSUMING WILL SEMINAIVE ALL ORDERS ***/
  int cutoff = bw;
  int size = 2*bw;

  double *rdata, *idata;
  double *rcoeffs, *icoeffs;

	//Initialize the spherical harmonics framework, if already done, do nothing
	initialize(bw);

  if(forward)
  {
	  //Input
	  rdata = real_data;
	  idata = imag_data;
	  //Output
	  rcoeffs = real_output;
	  icoeffs = imag_output;
  }
  else
  {
	  //Input
	  rcoeffs = real_data;
	  icoeffs = imag_data;
	  //Output
	  rdata = real_output;
	  idata = imag_output;
  }

  /* forward DCT */
  if(forward)
  {

	  /* make DCT plan -> note that I will be using the GURU
		 interface to execute these plans within the routines*/
	dctPlan = fftw_plan_r2r_1d( 2*bw, weights, rdata, FFTW_REDFT10, FFTW_ESTIMATE ) ;
	rank = 1 ;
	dims[0].n = 2*bw ;
	dims[0].is = 1 ;
	dims[0].os = 2*bw ;
	howmany_rank = 1 ;
	howmany_dims[0].n = 2*bw ;
	howmany_dims[0].is = 2*bw ;
	howmany_dims[0].os = 1 ;
  }
  else
  {
	  dctPlan = fftw_plan_r2r_1d( 2*bw, weights, rdata, FFTW_REDFT01, FFTW_ESTIMATE ) ;
	  rank = 1 ;
	  dims[0].n = 2*bw ;
	  dims[0].is = 2*bw ;
	  dims[0].os = 1 ;
	  howmany_rank = 1 ;
	  howmany_dims[0].n = 2*bw ;
	  howmany_dims[0].is = 1 ;
	  howmany_dims[0].os = 2*bw ;
  }

  /* forward fft */
  fftPlan = fftw_plan_guru_split_dft( rank, dims, howmany_rank, howmany_dims, rdata, idata, workspace, workspace+(4*bw*bw), FFTW_ESTIMATE );

  /* now do the forward spherical transform */
  if(forward)
  {

	  FST_semi_memo(rdata, idata,
		  rcoeffs, icoeffs,
		  bw,
		  seminaive_naive_table,
		  workspace,
		  0,
		  cutoff,
		  &dctPlan,
		  &fftPlan,
		  weights );

  }
  else
  {
	  InvFST_semi_memo(rcoeffs,icoeffs,
		  rdata, idata,
		  bw,
		  trans_seminaive_naive_table,
		  workspace,
		  0,
		  cutoff,
		  &dctPlan,
		  &fftPlan );
  }

	fftw_destroy_plan( fftPlan );
	fftw_destroy_plan( dctPlan );
}
