/* {{{ copyright */

/*  ztop - convert a single z value into p

    Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002 University of Oxford  */

/*  CCOPYRIGHT */

/* }}} */
/* {{{ defines */

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include "libprob.h"
#include <algorithm>

#define MINP 1e-15
using namespace std;
using namespace MISCMATHS;

/* }}} */
/* {{{ usage */

void usage(void)
{
  printf("Usage: ztop <z> [-2] [-g <number_of_resels>]\n");
  printf("-2 : use 2-tailed conversion (default is 1-tailed)\n");
  printf("-g : use GRF maximum-height theory instead of Gaussian pdf\n");
  exit(1);
}

/* }}} */

int main(int argc,char *argv[])
{
  int i, twotailed=0, grf=0;
  double p, z, nresels=0;

  if (argc<2)
    usage();

  z=atof(argv[1]);

  /* {{{ process args */

for (i=2; i<argc; i++)
{
  if (!strcmp(argv[i], "-2"))
    twotailed=1;
  else if (!strcmp(argv[i], "-g"))
    {
      grf=1;
      i++;
      if (argc<i+1)
    {
      printf("Error: no value given following -g\n");
      usage();
    }
      nresels=atof(argv[i]);
    }

  else
    usage();
}

/* }}} */

  if (twotailed)
    z=fabs(z);

  if (grf) {
    if (z<2)
      p=1; /* Below z of 2 E(EC) becomes non-monotonic */
    else
      p = nresels * 0.11694 * exp(-0.5*z*z)*(z*z-1);  /* 0.11694 = (4ln2)^1.5 / (2pi)^2 */
  }  else {
    p=1-ndtr(z);
  }

  if (twotailed)
    p*=2;

  p=min(p,1.0);

  if (p>0.0001) 
    printf("%f\n",p);
  else
    printf("%e\n",p);

  return(0);
}
