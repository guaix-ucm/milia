/*
 * Copyright 2008-2011 Sergio Pascual
 * 
 * This file is part of Milia
 * 
 * Milia is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Milia is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with Milia.  If not, see <http://www.gnu.org/licenses/>.
 * 
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#define HUBBLE 50
#define MATTER 1
#define LAMBDA 0
#include <iostream>
#include "milia/metric.h"
#include <cstdlib>
#include <cmath>

#include <popt.h>

using namespace std;

int main(int argc,const char **argv)
{
  double hubble(HUBBLE);
  double matter(MATTER);
  double lambda(LAMBDA);
  
  int lt(0),age(0),dl(0),vol(0),dc(0),dm(0),da(0),DM(0);
  //  Options option=NONE_OPTION;
  
  poptOption optionsTable[] = {
    {"hubble", 'h',POPT_ARG_DOUBLE,&hubble,0,
     "Hubble constant","50"},
    {"matter", 'm', POPT_ARG_DOUBLE, &matter,0,
     "Matter density","1"},
    {"lambda", 'l', POPT_ARG_DOUBLE, &lambda,0,
     "Lambda density","0"},
    {"dc",'\0',POPT_ARG_NONE,&dc,0,
     "Comoving distance (line of sight)",NULL},
    {"dm",'\0',POPT_ARG_NONE,&dm,0,
     "Comoving distance (transverse)",NULL},
    {"da",'\0',POPT_ARG_NONE,&da,0,
     "Angular diameter distance",NULL},
    {"dl",'\0',POPT_ARG_NONE,&dl,0,
     "Luminosity distance",NULL},
    {"DM",'\0',POPT_ARG_NONE,&DM,0,
     "Distance modulus",NULL},
    {"lt",'\0',POPT_ARG_NONE,&lt,0,
     "Look-back time",NULL},
    {"age",'\0',POPT_ARG_NONE,&age,0,
     "Age of the Universe",NULL},
    {"vol",'\0',POPT_ARG_NONE,&vol,0,
     "Comoving volume",NULL},
    POPT_AUTOHELP
    POPT_TABLEEND
  };
  
#ifdef PACKAGE
  poptContext optCon=poptGetContext(PACKAGE,argc,argv,optionsTable,0);
  poptReadDefaultConfig(optCon,0);
#else
  poptContext optCon=poptGetContext(NULL,argc,argv,optionsTable,0);
#endif
  //poptSetOtherOptionHelp(optCon,"[Cosmology] [cuantities] redshift");
  if(argc<2){
    poptPrintUsage(optCon,stderr,0);
    poptFreeContext(optCon);
    return 1;
  }
  
  /* Now do options processing, get portname */
  int c=poptGetNextOpt(optCon);
  
  if(c<-1){
    /* an error occurred during option processing */
    cerr<<poptBadOption(optCon,POPT_BADOPTION_NOALIAS)
	<<": "<<poptStrerror(c)<<endl;
    poptFreeContext(optCon);
    return 1;
  }
  
  const char* redshift;
  const milia::metric b(hubble,matter,lambda);

  while(redshift=poptGetArg(optCon)){
    const double z=atof(redshift);
    if(lt==1)
      cout<<b.lt(z)<<endl;
    if(age==1)
      cout<<b.age(z)<<endl;
    if(dl==1)
      cout<<b.dl(z)<<endl;
    if(vol==1)
      cout<<b.vol(z)<<endl;
    if(dc==1)
      cout<<b.dc(z)<<endl;
    if(dm==1)
      cout<<b.dm(z)<<endl;
    if(da==1)
      cout<<b.da(z)<<endl;
    if(DM==1)
      cout<<b.DM(z)<<endl;
  }
  
  poptFreeContext(optCon);
  return 0;
}
