#include <stdio.h>
#include <string.h>
#include "wrapper.h"
#include <ncarg/ngmath.h>

char shmsg[61];

NhlErrorTypes shgrid_W(void)
{
  int ier = 0;

  float *xi;
  int dsizes_xi[NCL_MAX_DIMENSIONS];
  float *yi;
  int dsizes_yi[NCL_MAX_DIMENSIONS];
  float *zi;
  int dsizes_zi[NCL_MAX_DIMENSIONS];
  float *fval;
  int dsizes_fval[NCL_MAX_DIMENSIONS];
  float *xo;
  int dsizes_xo[NCL_MAX_DIMENSIONS];
  float *yo;
  int dsizes_yo[NCL_MAX_DIMENSIONS];
  float *zo;
  int dsizes_zo[NCL_MAX_DIMENSIONS];

  int       has_missing_xi, has_missing_yi, has_missing_zi, has_missing_fval;
  int       has_missing_xo, has_missing_yo, has_missing_zo;
  NclScalar missing_xi, missing_yi, missing_zi, missing_fval;
  NclScalar missing_xo, missing_yo, missing_zo;

  float *uo;
  int i, j, ndims_uo, dsizes_uo[3], num_missing, num_points;

/*
 * Retrieve argument #0 (x coordinates).
 */
  xi = (float *) NclGetArgValue(
       0,
       7,
       NULL,
       dsizes_xi,
       &missing_xi,
       &has_missing_xi,
       NULL,
       2);
/*
 * Retrieve argument #1 (y coordinates).
 */
  yi = (float *) NclGetArgValue(
       1,
       7,
       NULL,
       dsizes_yi,
       &missing_yi,
       &has_missing_yi,
       NULL,
       2);
/*
 * Retrieve argument #2 (z coordinates).
 */
  zi = (float *) NclGetArgValue(
       2,
       7,
       NULL,
       dsizes_zi,
       &missing_zi,
       &has_missing_zi,
       NULL,
       2);
/*
 *  Retrieve argument #3 (input functional values)
 */
  fval = (float *) NclGetArgValue(
         3,
         7,
         NULL,
         dsizes_fval,
         &missing_fval,
         &has_missing_fval,
         NULL,
         2);
/*
 *  Retrieve argument #4 (output x coordinates)
 */
  xo = (float *) NclGetArgValue(
       4,
       7,
           NULL,
       dsizes_xo,
       &missing_xo,
       &has_missing_xo,
       NULL,
       2);
/*
 *  Retrieve argument #5 (output y coordinates)
 */
  yo = (float *) NclGetArgValue(
       5,
       7,
           NULL,
       dsizes_yo,
       &missing_yo,
       &has_missing_yo,
       NULL,
       2);
/*
 *  Retrieve argument #6 (output z coordinates)
 */
  zo = (float *) NclGetArgValue(
       6,
       7,
           NULL,
       dsizes_zo,
       &missing_zo,
       &has_missing_zo,
       NULL,
       2);

/*
 *  Check that all of the first four arguments are of the same size.
 */
  if( (dsizes_xi[0] != dsizes_yi[0]) ||
      (dsizes_yi[0] != dsizes_zi[0]) ||
      (dsizes_zi[0] != dsizes_fval[0]) ) {
        NhlPError(NhlFATAL,NhlEUNKNOWN,
           "shgrid: all arrays for the input functional values"
           " must be of the same size");
        return(NhlFATAL);
  }
  else {
    num_points = dsizes_xi[0];
    num_missing = 0;
  }

/*
 *  Check for missing values.  Each argument xi, yi, zi, and fval must 
 *  be checked separately, since it may be that _FillValue may be
 *  set for just one or two of the arguments.
 *
 *  If a missing value is found, ignore that input point by collapsing
 *  all of the input arrays.
 */
  if (has_missing_xi) {
    i = 0;
    lab0:
    if ((xi[i] == missing_xi.floatval) && (i < num_points)) {
      for (j = i; j < num_points - 1; j++) {
        xi[j] = xi[j+1];
        yi[j] = yi[j+1];
        zi[j] = zi[j+1];
        fval[j] = fval[j+1];
      }
      num_points--;
      if (num_points < 3) {
        NhlPError(NhlFATAL,NhlEUNKNOWN,
          "shgrid: missing values in the input data have reduced the number "
          "of valid points to less than 3.\n");
        return(NhlFATAL);
        
      }
      num_missing++;
      goto lab0;
    }
    i++;
    if (i < num_points) goto lab0;
  }
  if (has_missing_yi) {
    i = 0;
    lab1:
    if ((yi[i] == missing_yi.floatval) && (i < num_points)) {
      for (j = i; j < num_points - 1; j++) {
        xi[j] = xi[j+1];
        yi[j] = yi[j+1];
        zi[j] = zi[j+1];
        fval[j] = fval[j+1];
      }
      num_points--;
      if (num_points < 3) {
        NhlPError(NhlFATAL,NhlEUNKNOWN,
          "shgrid: missing values in the input data have reduced the number "
          "of valid points to less than 3.\n");
        return(NhlFATAL);
        
      }
      num_missing++;
      goto lab1;
    }
    i++;
    if (i < num_points) goto lab1;
  }
  if (has_missing_zi) {
    i = 0;
    lab2:
    if ((zi[i] == missing_zi.floatval) && (i < num_points)) {
      for (j = i; j < num_points - 1; j++) {
        xi[j] = xi[j+1];
        yi[j] = yi[j+1];
        zi[j] = zi[j+1];
        fval[j] = fval[j+1];
      }
      num_points--;
      if (num_points < 3) {
        NhlPError(NhlFATAL,NhlEUNKNOWN,
          "shgrid: missing values in the input data have reduced the number "
          "of valid points to less than 3.\n");
        return(NhlFATAL);
        
      }
      num_missing++;
      goto lab2;
    }
    i++;
    if (i < num_points) goto lab2;
  }
  if (has_missing_fval) {
    i = 0;
    lab3:
    if ((fval[i] == missing_fval.floatval) && (i < num_points)) {
      for (j = i; j < num_points - 1; j++) {
        xi[j] = xi[j+1];
        yi[j] = yi[j+1];
        zi[j] = zi[j+1];
        fval[j] = fval[j+1];
      }
      num_points--;
      if (num_points < 3) {
        NhlPError(NhlFATAL,NhlEUNKNOWN,
          "shgrid: missing values in the input data have reduced the number "
          "of valid points to less than 3.\n");
        return(NhlFATAL);
      }
      num_missing++;
      goto lab3;
    }
    i++;
    if (i < num_points) goto lab3;
  }

/*
 *  Issue warning if missing values have been detected.
 */
  if (num_missing > 0) {
    NhlPError(NhlWARNING,NhlEUNKNOWN,
      "shgrid: missing values in %d input points - those points ignored.",
      num_missing);
  }

/*
 *  Missing values not allowed in the output arrays.
 */
  if (has_missing_xo) {
    for (i = 0; i < dsizes_xo[0]; i++) {
      if (xo[i] == missing_xo.floatval) {
        NhlPError(NhlFATAL,NhlEUNKNOWN,
              "shgrid: output arrays not allowed to contain a missing value.");
        return(NhlFATAL);
      }
    }
  }
  if (has_missing_yo) {
    for (i = 0; i < dsizes_yo[0]; i++) {
      if (yo[i] == missing_yo.floatval) {
        NhlPError(NhlFATAL,NhlEUNKNOWN,
              "shgrid: output arrays not allowed to contain a missing value.");
        return(NhlFATAL);
      }
    }
  }
  if (has_missing_zo) {
    for (i = 0; i < dsizes_zo[0]; i++) {
      if (zo[i] == missing_zo.floatval) {
        NhlPError(NhlFATAL,NhlEUNKNOWN,
              "shgrid: output arrays not allowed to contain a missing value.");
        return(NhlFATAL);
      }
    }
  }

/*
 *  Call the C procedure.
 */
  uo = c_shgrid(dsizes_xi[0], xi, yi, zi, fval, dsizes_xo[0], 
                dsizes_yo[0], dsizes_zo[0], xo, yo, zo, &ier);
  if (ier != 0) {
    sprintf(shmsg, "shgrid: Error number %d.", ier);
    NhlPError(NhlFATAL, NhlEUNKNOWN, shmsg);
    return(NhlFATAL);
  }

  ndims_uo = 3;
  dsizes_uo[0] = dsizes_xo[0];
  dsizes_uo[1] = dsizes_yo[0];
  dsizes_uo[2] = dsizes_zo[0];
  return(NclReturnValue( (void *) uo, ndims_uo, dsizes_uo, NULL, NCL_float, 0));
}

NhlErrorTypes shgetnp_W(void)
{
  int i, j, ier = 0, num_missing, num_points, *k, ll;

  float *px;
  int dsizes_px[NCL_MAX_DIMENSIONS];
  float *py;
  int dsizes_py[NCL_MAX_DIMENSIONS];
  float *pz;
  int dsizes_pz[NCL_MAX_DIMENSIONS];

  float *xi;
  int dsizes_xi[NCL_MAX_DIMENSIONS];
  float *yi;
  int dsizes_yi[NCL_MAX_DIMENSIONS];
  float *zi;
  int dsizes_zi[NCL_MAX_DIMENSIONS];
 
  int *flag;
  int dsizes_flag[NCL_MAX_DIMENSIONS];

  int       has_missing_px, has_missing_py, has_missing_pz;
  NclScalar missing_px, missing_py, missing_pz;
  int       has_missing_xi, has_missing_yi, has_missing_zi;
  NclScalar missing_xi, missing_yi, missing_zi;

  static    int index, index_dims = 1;

/*
 * Retrieve argument #0 (x coordinate of reference point).
 */
  px = (float *) NclGetArgValue(
       0,
       7,
       NULL,
       dsizes_px,
       &missing_px,
       &has_missing_px,
       NULL,
       2);
/*
 * Retrieve argument #1 (y coordinate of reference point).
 */
  py = (float *) NclGetArgValue(
       1,
       7,
       NULL,
       dsizes_py,
       &missing_py,
       &has_missing_py,
       NULL,
       2);
/*
 * Retrieve argument #2 (z coordinate of reference point).
 */
  pz = (float *) NclGetArgValue(
       2,
       7,
       NULL,
       dsizes_pz,
       &missing_pz,
       &has_missing_pz,
       NULL,
       2);
/*
 * Retrieve argument #3 (x coordinates of input data).
 */
  xi = (float *) NclGetArgValue(
       3,
       7,
       NULL,
       dsizes_xi,
       &missing_xi,
       &has_missing_xi,
       NULL,
       2);
/*
 * Retrieve argument #4 (y coordinates of input data).
 */
  yi = (float *) NclGetArgValue(
       4,
       7,
       NULL,
       dsizes_yi,
       &missing_yi,
       &has_missing_yi,
       NULL,
       2);
/*
 * Retrieve argument #5 (z coordinates of input data).
 */
  zi = (float *) NclGetArgValue(
       5,
       7,
       NULL,
       dsizes_zi,
       &missing_zi,
       &has_missing_zi,
       NULL,
       2);
/*
 *  Retrieve argument #6 (call flag)
 */
  flag = (int *) NclGetArgValue(
         6,
         7,
         NULL,
         dsizes_flag,
         NULL,
         NULL,
         NULL,
         2);

/*
 *  Check that the input arrays are the same size.
 */
  if( (dsizes_xi[0] != dsizes_yi[0]) ||
      (dsizes_yi[0] != dsizes_zi[0])) {
        NhlPError(NhlFATAL,NhlEUNKNOWN,
           "shgetnp: all arrays for the input functional values"
           " must be of the same size");
        return(NhlFATAL);
  }
  else {
    num_points = dsizes_xi[0];
    num_missing = 0;
  }

/*
 *  Check for missing values.  Each argument xi, yi, and zi must 
 *  be checked separately, since it may be that _FillValue may be
 *  set for just one or two of the arguments.
 *
 *  If a missing value is found, ignore that input point by collapsing
 *  all of the input arrays.
 */
  if (has_missing_xi) {
    i = 0;
    lab0:
    if ((xi[i] == missing_xi.floatval) && (i < num_points)) {
      for (j = i; j < num_points - 1; j++) {
        xi[j] = xi[j+1];
        yi[j] = yi[j+1];
        zi[j] = zi[j+1];
      }
      num_points--;
      if (num_points < 3) {
        NhlPError(NhlFATAL,NhlEUNKNOWN,
          "shgrid: missing values in the input data have reduced the number "
          "of valid points to less than 3.\n");
        return(NhlFATAL);
        
      }
      num_missing++;
      goto lab0;
    }
    i++;
    if (i < num_points) goto lab0;
  }
  if (has_missing_yi) {
    i = 0;
    lab1:
    if ((yi[i] == missing_yi.floatval) && (i < num_points)) {
      for (j = i; j < num_points - 1; j++) {
        xi[j] = xi[j+1];
        yi[j] = yi[j+1];
        zi[j] = zi[j+1];
      }
      num_points--;
      if (num_points < 3) {
        NhlPError(NhlFATAL,NhlEUNKNOWN,
          "shgrid: missing values in the input data have reduced the number "
          "of valid points to less than 3.\n");
        return(NhlFATAL);
        
      }
      num_missing++;
      goto lab1;
    }
    i++;
    if (i < num_points) goto lab1;
  }
  if (has_missing_zi) {
    i = 0;
    lab2:
    if ((zi[i] == missing_zi.floatval) && (i < num_points)) {
      for (j = i; j < num_points - 1; j++) {
        xi[j] = xi[j+1];
        yi[j] = yi[j+1];
        zi[j] = zi[j+1];
      }
      num_points--;
      if (num_points < 3) {
        NhlPError(NhlFATAL,NhlEUNKNOWN,
          "shgrid: missing values in the input data have reduced the number "
          "of valid points to less than 3.\n");
        return(NhlFATAL);
        
      }
      num_missing++;
      goto lab2;
    }
    i++;
    if (i < num_points) goto lab2;
  }

/*
 *  Issue warning if missing values have been detected.
 */
  if (num_missing > 0) {
    NhlPError(NhlWARNING,NhlEUNKNOWN,
      "shgrid: missing values in %d input points - those points ignored.",
      num_missing);
  }

/*
 *  Call the C procedure.
 */
  index = c_shgetnp(*px, *py, *pz, dsizes_xi[0], xi, yi, zi, *flag, &ier);
  if (ier != 0) {
    sprintf(shmsg, "shgetnp: Error number %d.", ier);
    NhlPError(NhlFATAL, NhlEUNKNOWN, shmsg);
    return(NhlFATAL);
  }
  return(NclReturnValue( (void *) (&index), 1, &index_dims, NULL, NCL_int, 0));
}

NhlErrorTypes shsetp_W(void)
{

  char  *arg1, *cval;
  int   numpi, numpf, numpc, i;

/*
 *  List the integer parameter names.  To add new ones,
 *  all that needs to be done is add the names to this list.
 */
  char *params_i[] = {"nls", "NLS", "nfl", "NFL", "ncl", "NCL"};

/*
 * Input array variables
 */
  string *pname;
  int dsizes_pname[NCL_MAX_DIMENSIONS];
  void *pvalue;
  int dsizes_pvalue[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_pname, type_pvalue;

/*
 * Retrieve argument #1
 */
  pname = (string *) NclGetArgValue(
          0,
          2,
          NULL,
          dsizes_pname,
          NULL,
          NULL,
          &type_pname,
          2);

  arg1 = NrmQuarkToString(*pname);
 
/*
 *  Check to see if the parameter name is valid.
 */
  numpi = sizeof(params_i)/sizeof(void *);
  for (i = 0; i < numpi; i++) {
    if (!strncmp(arg1, params_i[i], strlen(params_i[i]))) {
      goto OK_NAME;
    }
  }
  NhlPError(NhlFATAL, NhlEUNKNOWN, "shsetp: unrecognized parameter name");
  return(NhlFATAL);

/*
 * Retrieve argument #2
 */
OK_NAME: pvalue = (void *) NclGetArgValue(
           1,
           2,
           NULL,
           dsizes_pvalue,
           NULL,
           NULL,
           &type_pvalue,
           2);

/*
 *  Process the parameter if it has an integer value.
 */
  if (type_pvalue == NCL_int) {
    for (i = 0; i < numpi; i++) {
      if (!strncmp(arg1, params_i[i], strlen(params_i[i]))) {
        c_shseti(arg1, *((int *) pvalue));
        return(NhlNOERROR);
      }
    }
    NhlPError(NhlFATAL, NhlEUNKNOWN, "The specified value for the parameter "
              "has an invalid type");
    return(NhlFATAL);
  }
  else {
    NhlPError(NhlFATAL, NhlEUNKNOWN, "The specified value for the parameter "
              "has an invalid type");
    return(NhlFATAL);
  }
}

NhlErrorTypes shgetp_W(void)
{
/*
 *  Get values for shgrid parameters.
 */

  char  *arg1, *cval;
  int   numpi, numpf, numpc, i;
  string *pvalue, *qvalue;

/*
 *  List the integer parameter names.  To add new ones,
 *  all that needs to be done is add the names to this list.
 */
  char *params_i[] = {"nls", "NLS", "nfl", "NFL", "ncl", "NCL"};

/*
 * Input array variable
 */
  string *pname;
  int dsizes_pname[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_pname;
  float *fval;
  int *ival;
  int ret_size = 1; 

/*
 * Retrieve argument #1
 */
  pname = (string *) NclGetArgValue(
          0,
          1,
          NULL,
          dsizes_pname,
          NULL,
          NULL,
          &type_pname,
          2);

  arg1 = NrmQuarkToString(*pname);

/*
 *  Check to see if the parameter name is valid.
 */
  numpi = sizeof(params_i)/sizeof(void *);
  for (i = 0; i < numpi; i++) {
    if (!strncmp(arg1, params_i[i], strlen(params_i[i]))) {
      goto OK_NAME;
    }
  }
  NhlPError(NhlFATAL, NhlEUNKNOWN, "shgetp: unrecognized parameter name");
  return(NhlFATAL);

/*
 *  Process the parameter if it has an integer value.
 */
OK_NAME:  
  for (i = 0; i < numpi; i++) {
    if (!strncmp(arg1, params_i[i], strlen(params_i[i]))) {
      ival = (int *) calloc(1,sizeof(int));
      *ival = c_shgeti(arg1);
      return(NclReturnValue( (void *) ival, 1, &ret_size, NULL, NCL_int, 0));
    }
  }
  NhlPError(NhlFATAL, NhlEUNKNOWN, "shgetp: impossible to get this message");
  return(NhlFATAL);
}
