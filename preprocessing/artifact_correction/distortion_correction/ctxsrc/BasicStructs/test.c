// test.c : Defines the entry point for the console application.

//



#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <time.h>
#include "basicstructs.h"


int main(int argc, char* argv[])
{
	
   dmat *m = dMatfromFile(4, 4, "Matl2vol_rb.txt");
   printf("\n Input:\n");
   dMatPrint(m, NULL);
   printf("\n Inverse:\n");
   dmat *mi = dMatInv(m);
   dMatPrint(mi, NULL);
   
   printf("\n matrix Mutplication \n");
   dmat *mm = dMatMultiply(mi, m, 1., 0.);
   dMatPrint(mm, NULL);
   
   printf("\n matrix add \n");
   dmat *me = dMatEye(4);
   dmat *ma = dMatAdd(me, mm, -1.);
   dMatPrint(ma, NULL);
   
   printf("\n norm2 = %f", dMatNorm2(ma, -1));
   printf("\n Dot = %f", dMatDot(ma, ma, -1));
   
   dmat *c3 = dMatGetCol(m, 3);
   printf("\n colume 3 \n");
   dMatPrint(c3, NULL);
   printf("\n c3 norm2_3 = %f", dMatNorm2(c3, 3));
   printf("\n c3 Dot_3 = %f", dMatDot(c3, c3, 3));
   printf("\n c3 norm2 = %f", dMatNorm2(c3, -1));
   printf("\n c3 Dot = %f", dMatDot(c3, c3, -1));
   
   dmat *r2 = dMatGetRow(m, 2);
   printf("\n row 2 \n");
   dMatPrint(r2, NULL);
   printf("\n r2 norm2 = %f", dMatNorm2(r2, 3));
   printf("\n r2 Dot = %f", dMatDot(r2, r2, 3));
   
    printf("\n scale row 2 \n");
    dmat *sr2 = dMatScal(r2, 2.);
    dMatPrint(sr2, NULL);
   
   dmat *r1 = dMatGetRow(m, 1);
   printf("\n row 1 \n");
   dMatPrint(r1, NULL);
   printf("\n Angle between r1 r2 =%f", dMatComputeAngle(r1, r2, -1));
   printf("\n Angle between r1 r2 =%f", dMatComputeAngle(r1, r2, 3));
   
   dmat *cr12 = dMatCross3(r1, r2);
    printf("\n cross of row 1 and row2 \n");
   dMatPrint(cr12, NULL);
   
  
  dMatDel(cr12);
   dMatDel(r1);
   dMatDel(sr2);
   dMatDel(c3);
   dMatDel(r2);
   dMatDel(ma);
   dMatDel(me);
   dMatDel(mi);
   dMatDel(mm);
   dMatDel(m);
   return 0;	
}



