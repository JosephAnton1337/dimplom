#include <stdio.h>
#include <math.h>
#include "consts.h"
#define DOPRITMP "tmp.out"
void dopri8(double x,double *y,double xend,double eps,double hmax,double h);
void dopri8_print(const char *filename,double x,double *y,double xend,double eps,double hmax,double h);
