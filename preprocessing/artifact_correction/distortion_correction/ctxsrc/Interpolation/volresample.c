///////////////////////////////////////////////////////////////////////
// This source code is Copyright ©2007-2010 CorTechs Labs
// Inc. Permission is granted by CorTechs Labs to copy, distribute and
// modify this source code, provided that it is used for research
// purposes only, and only by not-for-profit organizations, and
// further provided that:
// 
//    * this copyright notice appears in all copies.
//    * you indemnify and hold harmless CorTechs Labs and its
//      successors and associates from any and all liability from any
//      use of the information.
//    * you ensure that these conditions apply to any further
//      redistribution of the source code or derived software
// 
// For commercial entities and for commercial applications, the right
// to copy, distribute or modify this source code is retained by
// CorTechs Labs, and must be obtained through completion of a written
// license agreement with CorTechs Labs Inc.
///////////////////////////////////////////////////////////////////////

/**
 * Volume resample 
*/

#include <time.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "interpolation.h"

void usage()
{
    printf("\n\n");
    printf("\n Volume Resample");
    printf("\n Author: Gen-Nan Chen, CorTechs Labs Inc.");
    printf("\n email: gnchen@cortechs.net");
    printf("\n \n");
    printf("\n usage:");
    printf("\n      volresample [options] volinfn voloutfn");
    printf("\n");
    printf("\n options:");
    printf("\n      -h: Help");
    printf("\n      -M: 4X4 Regsitration Matrix file (default=I)");
    printf("\n      -t: Output volume template mgh file");
    printf("\n          (default is 256^3, 1mm^3 in coronal View, .i.e atlas space)");
    printf("\n      -m: Interpolation Method (0: Nearest neighbor, 1: Linear, 2: Cubic (default))");
    printf("\n      -e: Edge effect  (0: Mirror , 1: Zeropadding(default))");
    printf("\n      -n: Allow Negative value");
    printf("\n\n");
}


int main(int argc, char* argv[])
{
    //Input
    char  *infn=NULL, *outfn=NULL, *templatefn=NULL, *regfn=NULL;
    hmat mr2o;
    fvol *vol, *volr;
    int method =2;
    int edge =1;
    int bpositive = 1;
    char c;
    time_t ti1, ti2, tt1, tt2;
    tt1 = time(NULL);
    //getoption
    while ((c = getopt (argc, argv, "M:t:m:e:nh")) != -1)
    {
       switch(c)
       { 
       case 'M':
            regfn = optarg;
            break;
       case 't':
            templatefn = optarg;
            break;
       case 'm':
            method = atoi(optarg);
            break;
       case 'e':
            edge = atoi(optarg);
            break;
       case 'n':
            bpositive = 0;
            break;
       case 'h':
       default:
            usage();
            return 0;
       }         
    }
    
    //check input file
    infn = argv[optind];
    if (infn == NULL)
    {
        printf("\n Error : Need input file !!!");
        usage();
        return 1;
    }
    else
    {
        printf("\n Reading input file %s ....", infn);
        ti1=time(NULL);
        vol = fVolFromMGH(infn, 1);
        ti2=time(NULL);
        printf("DONE in %10.2f secs", difftime(ti2,ti1));
        if (vol ==NULL)
        {
            printf("\n Error : CANNOT read input file %s !!!\n", infn);
            return 1;
        }
    }
    
    outfn = argv[optind+1];
    if (outfn == NULL)
    {
        printf("\n Error : Need output filename !!!");
        usage();
        return 1;
    }
    
    //check template file
    if (templatefn == NULL)
    {
       printf("\n Allocate atlas-like volume");
       volr = fVolNewAtlVol();
       if (volr ==NULL)
        {
            printf("\n Error : CANNOT allocate resample volume !!!\n");
            return 1;
        } 
    }
    else
    {
        printf("\n Reading template file %s", templatefn);
        volr = fVolFromMGH(templatefn, 0);
        if (volr ==NULL)
        {
            printf("\n Error : CANNOT read template file %s !!!", templatefn);
            return 1;
        }
    }
    
    //check regstration file
    if (regfn)
    {
        if (hMatFromFile(regfn, &mr2o) == 0)
        {
            printf("\n Error : CANNOT read regsitration file %s !!! \n", regfn);
            return 1;
        }
    }
    
    hMatPrint(&mr2o);
    
    // check interp method
    if (method >2)
    {
        printf("\n Error : Wrong interpolation method !!! \n");
        usage ();
        return 1;
    }
    
    //check edge effect
    if (edge >1)
    {
        printf("\n Error : Wrong Edge effect !!! \n");
        usage ();
        return 1;
    }
    
    INTERPMETHOD iMethod = NEAREST;
    switch(method) {
    case 0: iMethod = NEAREST;
    case 1: iMethod = LINEAR;
    case 2: iMethod = CUBIC;
    }
    
    INTERPEDGE iEdge = MIRROR;
    switch(edge) {
    case 0: iEdge = MIRROR;
    case 1: iEdge = ZEROPAD;
    case 2: iEdge = THROWAWAY;
    }
    
    printf("\n Reslice Volume....");
    ti1=time(NULL);
    //fVolResample(vol, volr, &mr2o, method, edge, bpositive);
    fVolResample(vol, volr, &mr2o, iMethod, iEdge, bpositive);
    ti2=time(NULL);
    printf("DONE in %10.2f secs", difftime(ti2,ti1));
    
    printf("\n Writing Volume to %s....", outfn);
    ti1=time(NULL);
    fVol2MGH(volr, outfn, FLOAT);
    ti2=time(NULL);
    printf("DONE in %10.2f secs", difftime(ti2,ti1));
    
    
    fVolDelete(vol);
    fVolDelete(volr);  
    
    tt2 = time(NULL);  
    printf("\n volresample took  %10.2f secs totally\n\n", difftime(tt2,tt1));
    return 0;
}
