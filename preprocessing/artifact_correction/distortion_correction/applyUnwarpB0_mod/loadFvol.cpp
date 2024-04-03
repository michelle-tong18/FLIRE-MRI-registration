#include <sys/types.h>           // for the system call stat()
#include <sys/stat.h>            // for the system call stat()
#include <string>
#include <iostream>
#include "imgvol.h"              // CTX stuff. fvol
#include "cstring_copy.h"

using std::endl;
using std::clog;
using std::string;


fvol*
loadFvol(const string fileName) {
  
  if( fileName.empty() )
    return 0;
  
  struct stat  statbuf;                            // To see status of input files
  
  int status ( stat(fileName.c_str(), &statbuf) ); // Check that file fileName exists, exit if not
  if(status != 0 || !(statbuf.st_mode & S_IFREG) )  {
    clog << __FILE__ << ": " << fileName << " does not exist. Exiting...\n" << endl;
    exit(0);
  }
  
  char* fileNameNonConst ( cstring_copy(const_cast<char*>(fileName.c_str())) );
  int bLoadData ( 1 );
  
  fvol* vol ( fVolFromMGH(fileNameNonConst, bLoadData) );
  
  delete fileNameNonConst;
  
  return vol;
  
}
