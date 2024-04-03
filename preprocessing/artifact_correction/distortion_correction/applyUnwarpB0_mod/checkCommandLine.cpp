#include <iostream>
#include <vector>
#include <string>
#include <algorithm>   // for find()

using std::vector;
using std::string;
using std::clog;
using std::endl;
using std::find;


bool
checkCommandLine(char** argv, const vector<string>& commandLine) {
  
  bool errorFlag ( false );
  vector<string>::const_iterator pos;
  
  bool f (false); // forward image
  pos = find(commandLine.begin(), commandLine.end(), "-f");
  if(pos != commandLine.end())
    f = true;
  
  bool r (false); // reverse image
  pos = find(commandLine.begin(), commandLine.end(), "-r");
  if(pos != commandLine.end())
    r = true;

  bool d (false); // reverse image
  pos = find(commandLine.begin(), commandLine.end(), "-d");
  if(pos != commandLine.end())
    d = true;
  
  bool ip (false); // input parameters file
  pos = find(commandLine.begin(), commandLine.end(), "-ip");
  if(pos != commandLine.end())
    ip = true;
  
  if( (!ip && !d) || (!ip && !f && !r) || (r && !d && !ip) || (f && !d && !ip) )
    errorFlag = true;
  
  return errorFlag;
}
