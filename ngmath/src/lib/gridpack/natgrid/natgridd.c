/*
 *  Copyright (C) 1991-1996
 *          by
 *    David F. Watson
 */

#include "nngheadd.h"

double **c_natgridd(int n, double x[], double y[], double z[],
                   int nxi, int nyi, double xi[], double yi[], int *ier)
{  
   double **data_out;

   *ier = 0;

   if (single_point == 0)
   {
      asflag = 1;
      Initialized(n, x, y, nxi, nyi, xi, yi);

      if (ReadData(n,x,y,z) != 0)
      {
         *ier = error_status;
         return ( (double **) NULL);
      }
   }

   if (adf)
   {
      CircOut();
      if (error_status)
      {
         *ier = error_status;
         return ( (double **) NULL);
      }
   }
   if (igrad)
   {
      Gradient();
      if (error_status)
      {
         *ier = error_status;
         return ( (double **) NULL);
      }
   }

   data_out = MakeGrid(nxi, nyi, xi, yi);
   if (error_status)
   {
      *ier = error_status;
      return ( (double **) NULL);
   }

   if (single_point == 0)
   {
      Terminate();
   }

   return (data_out);
}
void Initialized(int n, double x[], double y[], int nxi, int nyi,
                double xi[], double yi[])
{

   double xil, xir, yib, yit;

   error_status = 0;
   datcnt       = 0;
   magx_orig    = magx;
   magy_orig    = magy;
   magz_orig    = magz;
   iscale       = 0;
   magx_auto    = 1.;
   magy_auto    = 1.;
   magz_auto    = 1.;

/*
 *  Find the limits of the output array.
 */
   xstart       = armind(nxi, xi);
   xend         = armaxd(nxi, xi);
   ystart       = armind(nyi, yi);
   yend         = armaxd(nyi, yi);

/*
 *  Find the limits of the input array.
 */
   xil          = armind(n, x);
   xir          = armaxd(n, x);
   yib          = armind(n, y);
   yit          = armaxd(n, y);

/*
 *  As the default (that is, unless horizontal and vertical overlaps
 *  have been specifically set by the user) choose the overlap values
 *  as the smallest values that will make all input data points included
 *  in the overlap region.
 */
   if (horilap EQ -1.) {
     if ( (xstart >= xil) && (xend <= xir) ) {
       horilap = 1.01 * (((xstart-xil) < (xir-xend)) ?
                          (xir-xend) : (xstart-xil));
     }
     else if ( (xstart >= xil) && (xend >= xir) ) {
       horilap = 1.01 * (xstart-xil);
     }
     else if ( (xstart <= xil) && (xend <= xir) ) {
       horilap = 1.01 * (xir-xend);
     }
     else if ( (xstart <= xil) && (xir <= xend) ) {
       horilap = 0.;
     }
   }
   if (horilap <= EPSILON) {
     horilap = 0.01 * (xend - xstart);
   }
   if (vertlap EQ -1.) {
     if ( (yib <= ystart) && (yend <= yit) ) {
       vertlap = 1.01 * (((ystart-yib) < (yit-yend)) ?
                          (yit-yend) : (ystart-yib));
     }
     else if ( (ystart <= yib) && (yend <= yit) ) {
       vertlap = 1.01 * (yit-yend);
     }
     else if ( (yib <= ystart) && (yit <= yend) ) {
       vertlap = 1.01 * (ystart-yib);
     }
     else if ( (ystart <= yib) && (yit <= yend) ) {
       vertlap = 0.;
     }
   }
   if (vertlap <= EPSILON) {
     vertlap = 0.01 * (yend - ystart);
   }
}

void Terminate()
{
   magx = magx_orig;
   magy = magy_orig;
   magz = magz_orig;
}

double armind(int num, double *x)
{
   int i;
   float amin;
   amin = x[0];
   for (i = 1 ; i < num ; i++)
      if (x[i] < amin) amin = x[i];
   return(amin);
}
double armaxd(int num, double *x)
{
   int i;
   float amax;
   amax = x[0];
   for (i = 1 ; i < num ; i++)
      if (x[i] > amax) amax = x[i];
   return(amax);
}
