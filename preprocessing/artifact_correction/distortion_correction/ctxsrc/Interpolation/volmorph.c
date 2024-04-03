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
 * Volume Nonlinear Morph
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
    printf("\n Volume Nonlinear Morph");
    printf("\n Author: Gen-Nan Chen, CorTechs Labs Inc.");
    printf("\n email: gnchen@cortechs.net");
    printf("\n \n");
    printf("\n usage:");
    printf("\n      volmorph [options] volinfn dLfn dPfn dHfn voloutfn\n");
    printf("\n          volinfn: input volume file");
    printf("\n          dLfn: vector field in L coordinate volume file");
    printf("\n          dHfn: vector field in P coordinate volume file");
    printf("\n          dHfn: vector field in H coordinate volume file");
    printf("\n          voloufn: output volume file");
    printf("\n");
    printf("\n options:");
    printf("\n      -h: Help");
    printf("\n      -M: Mask for vector field file. Need to be the same as volr");
    printf("\n          (default: every voxel count)");
    printf("\n      -t: Output volume template mgh file");
    printf("\n          (default is in the same as dLfn)");
    printf("\n      -m: Interpolation Method (0: Nearest neighbor, 1: Linear, 2: Cubic (default))");
    printf("\n      -e: Edge effect  (0: Mirror , 1: Zeropadding(default))");
    printf("\n      -j: Apply jacobian");
    printf("\n      -n: Allow Negative value");
    printf("\n\n");
}


int main(int argc, char* argv[])
{
     //Input
    char  *infn=NULL, *outfn=NULL, *templatefn=NULL, *dLfn=NULL, *dHfn=NULL, *dPfn=NULL, *maskfn=NULL;
    fvol *vol, *volr, *dL, *dP, *dH;
    ucvol *mask;
    int method =2;
    int edge =1;
    int bpositive = 1;
    int bjacbian = 0;
    int bsame = 1;
    char c;
    time_t ti1, ti2, tt1, tt2;
    tt1 = time(NULL);
    //getoption
    while ((c = getopt (argc, argv, "M:t:m:e:njh")) != -1)
    {
       switch(c)
       { 
       case 'M':
            maskfn = optarg;
            break;
       case 'j':
            bjacbian = 1;
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
    
     //check input dL file
    dLfn = argv[optind+1];
    if (dLfn == NULL)
    {
        printf("\n Error : Need vector field in L coordinate volume file !!!");
        usage();
        return 1;
    }
    else
    {
        printf("\n Reading vector field in L coordinate file %s ....", dLfn);
        ti1=time(NULL);
        dL = fVolFromMGH(dLfn, 1);
        ti2=time(NULL);
        printf("DONE in %10.2f secs", difftime(ti2,ti1));
        if (dL ==NULL)
        {
            printf("\n Error : CANNOT read ivector field in L coordinate file %s !!!\n", dLfn);
            return 1;
        }
    }
    
    //check input dL file
    dPfn = argv[optind+2];
    if (dPfn == NULL)
    {
        printf("\n Error : Need vector field in P coordinate volume file !!!");
        usage();
        return 1;
    }
    else
    {
        printf("\n Reading vector field in P coordinate file %s ....", dPfn);
        ti1=time(NULL);
        dP = fVolFromMGH(dPfn, 1);
        ti2=time(NULL);
        printf("DONE in %10.2f secs", difftime(ti2,ti1));
        if (dP ==NULL)
        {
            printf("\n Error : CANNOT read ivector field in P coordinate file %s !!!\n", dPfn);
            return 1;
        }
    }
    
    //check input dL file
    dHfn = argv[optind+3];
    if (dHfn == NULL)
    {
        printf("\n Error : Need vector field in H coordinate volume file !!!");
        usage();
        return 1;
    }
    else
    {
        printf("\n Reading vector field in H coordinate file %s ....", dHfn);
        ti1=time(NULL);
        dH = fVolFromMGH(dHfn, 1);
        ti2=time(NULL);
        printf("DONE in %10.2f secs", difftime(ti2,ti1));
        if (dH ==NULL)
        {
            printf("\n Error : CANNOT read ivector field in H coordinate file %s !!!\n", dHfn);
            return 1;
        }
    }
    
    outfn = argv[optind+4];
    if (outfn == NULL)
    {
        printf("\n Error : Need output filename !!!");
        usage();
        return 1;
    }
    
    //check template file
    if (templatefn == NULL)
    {
       printf("\n Allocate dL-like volume");
       volr = fVolCopy(dL, 0);
       //fVolPrint(volr);
       bsame =1;
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
        bsame = 0;
        if (volr ==NULL)
        {
            printf("\n Error : CANNOT read template file %s !!!", templatefn);
            return 1;
        }
    } 
    
    //check mask file
    if (maskfn == NULL)
    {
       printf("\n Allocate volume for mask, same size as volr");
       mask = ucVolNew(volr->info.dim, &(volr->info.M_vxl2lph), &(volr->info.M_pat2grad), NULL);
       ucVolFill(mask, 1);
       if (volr ==NULL)
        {
            printf("\n Error : CANNOT allocate mask volume !!!\n");
            return 1;
        } 
    }
    else
    {
        printf("\n Reading mask file %s", maskfn);
        mask = ucVolFromMGH(maskfn, 1);
        if (volr ==NULL)
        {
            printf("\n Error : CANNOT read mask file %s !!!", templatefn);
            return 1;
        }
    } 
    
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
    
    printf("\n Morph Volume....");
    ti1=time(NULL);
    fvols *dfield = fVolsNewVectorField(dL, dP, dH, 0);
    //fVolMorph(vol, volr, bsame, dfield, mask, method, edge, bpositive, bjacbian);
    fVolMorph(vol, volr, bsame, dfield, mask, iMethod, iEdge, bpositive, bjacbian);
    ti2=time(NULL);
    printf("DONE in %10.2f secs", difftime(ti2,ti1));
    
    printf("\n Writing Volume to %s....", outfn);
    ti1=time(NULL);
    fVol2MGH(volr, outfn, FLOAT);
    ti2=time(NULL);
    printf("DONE in %10.2f secs", difftime(ti2,ti1));
    
    
    fVolDelete(vol);
    fVolDelete(dL);
    fVolDelete(dP);
    fVolDelete(dH);
    fVolsDelete(dfield);
    ucVolDelete(mask);
    fVolDelete(volr);  
    
    tt2 = time(NULL);  
    printf("\n volmorph took  %10.2f secs totally\n\n", difftime(tt2,tt1));
    
    return 0;
}
