// ==================================================
// Copyright (c) 2010 Dominic Holland and Anders Dale
// All rights reserved
// ==================================================

////////////////////////////////////////////////////////////////////
// Load pareameters from an input-parameters file into an
// InputParams class.
//
// One improvement to be made here would be to include checking that
// ALL REQUIRED parameters were actually loaded...
////////////////////////////////////////////////////////////////////

#include <fstream>        // to open files for reading & writing
//#include <iostream>
//#include <string>
#include <sstream>        // for stringstream
#include <vector>
#include <map>
#include <algorithm>      // for copy() and find()
#include <iterator>       // for ostream_iterator
#include <cctype>         // for isspace()
#include "inputParams.h"
#include "split.h"        // split a string of words into a vector of words
#include "checkCommandLine.h"

using std::ifstream;
using std::ostream;
using std::cerr;
using std::endl;
using std::string;
using std::map;
using std::vector;
using std::stringstream;
using std::exit;
using std::find;

void usage();


void
InputParams::readParameters(int argc, char** argv) {
  
  // Some parameters and names optionally can be supplied on the comand line.
  // Command line values override anyother values.
  vector<string> commandLine;
  
  // Load command line into a string.
  for(int i=0; i<argc; ++i) {
    commandLine.push_back(argv[i]);
    std::clog << commandLine[i] << endl;
  }
  
  if( checkCommandLine(argv, commandLine) ) {
    usage();
    exit(0);
  }
  
  vector<string>::iterator pos;
  
  bool fFlag (false);     // If set, commandline value will override value in file, if provided.
  pos = find(commandLine.begin(), commandLine.end(), "-f");
  if(pos != commandLine.end()) {
    forwardImageInFileName = *(++pos);
    fFlag = true;
  }

  bool rFlag (false);
  pos = find(commandLine.begin(), commandLine.end(), "-r");
  if(pos != commandLine.end()) {
    reverseImageInFileName = *(++pos);
    rFlag = true;
  }
  
  bool foFlag (false);     // If set, commandline value will override value in file, if provided.
  pos = find(commandLine.begin(), commandLine.end(), "-fo");
  if(pos != commandLine.end()) {
    forwardImageOutFileName = *(++pos);
    foFlag = true;
  }
  
  bool roFlag (false);
  pos = find(commandLine.begin(), commandLine.end(), "-ro");
  if(pos != commandLine.end()) {
    reverseImageOutFileName = *(++pos);
    roFlag = true;
  }
    
  bool dFlag (false);      // If set, commandline value will override value in file, if provided.
  pos = find(commandLine.begin(), commandLine.end(), "-d");
  if(pos != commandLine.end()) {
    displacementFieldFileName = *(++pos);
    dFlag = true;
  }
  
  bool odFlag (false);
  pos = find(commandLine.begin(), commandLine.end(), "-od");
  if(pos != commandLine.end()) {
    outDir = *(++pos);
    odFlag = true;
  }
  
  // applyJacobianIntensityCorrection
  bool jicFlag (false);// If set, commandline value will override value in file, if provided.
  pos = find(commandLine.begin(), commandLine.end(), "-jic");
  if(pos != commandLine.end()) {
    jic = true;             // Default set to false
    ++pos;                  // don't forget to iterate on...!
    jicFlag = true;
  }

  
  ////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////
  
  // Open the file referred to by argv[1]
  //ifstream fParamsIn(argv[1]);
  //if(!fParamsIn) errorReadFile(" Cannot open input file", argv[1]);
  
  // Access the file referred to by -ip, or return.
  string inputParamsFileName;
  pos = find(commandLine.begin(), commandLine.end(), "-ip");
  if(pos != commandLine.end())
    inputParamsFileName = *(++pos);
  else
    return; // Use default initializations.
  
  ifstream fParamsIn(inputParamsFileName.c_str());
  if(!fParamsIn) errorReadFile(" Cannot open input file", inputParamsFileName);
  
  map<string, string> strMap;
  string s;
  int count ( 0 );

  while(getline(fParamsIn, s))  {
    
    if( !s.size() ) continue;              // skip empty lines
    
    vector<string> v = split(s);
    if( !v.size() ) continue;              // skip pure white-space lines
    
    if( v.size() < 3 ) continue;           // skip empty variables (nothing after the "=")
                                           // Yes, I know, this could absorb the previous line.
    if(v[0].substr(0,2) == "//") continue; // skip comment lines
    if(v[2] == "//") continue;             // skip lines whose variables are unset (there's only a comment preceded by "//" after the "=")
    
    // Now, start selecting variables...
    if     (!fFlag  && v[0] == "forwardImageInFileName")       { forwardImageInFileName       = v[2];  count++; }
    else if(!rFlag  && v[0] == "reverseImageInFileName")       { reverseImageInFileName       = v[2];  count++; }
    else if(!foFlag && v[0] == "forwardImageOutFileName")      { forwardImageOutFileName      = v[2];  count++; }
    else if(!roFlag && v[0] == "reverseImageOutFileName")      { reverseImageOutFileName      = v[2];  count++; }
    else if(!dFlag  && v[0] == "displacementFieldFileName")    { displacementFieldFileName    = v[2];  count++; }
    else if(!odFlag && v[0] == "outDir")                       { outDir                       = v[2];  count++; } // Set default to current dir.
    else if(!jicFlag && v[0] == "jic")  {                      // v[2] will be "0" or "1". Translate to true or false.
      stringstream ssn(v[2]);                                  // stringstream --> bool
      ssn >> jic;  count++;                                    // Default set to false
    }
    
#if 0
    else {
      cerr << __FILE__ << ":  invalid input file  " << argv[1] << endl;
      cerr << "line " << s << endl;
      cerr << "v[0] " << v[0] << endl;
      exit(EXIT_FAILURE);
    }
#endif
  }
  
  //  if(!fParamsIn.good())   // check if input failed
  //    errorReadFile(fParamsIn, argv[1], count);
}


ostream&
InputParams::print(ostream& os) const {
  os << "\nHave parameters:\n\n";
  
  os << "forwardImageInFileName       = " << forwardImageInFileName       << '\n';	 
  os << "reverseImageInFileName       = " << reverseImageInFileName       << '\n'; 
  
  os << "forwardImageOutFileName      = " << forwardImageOutFileName      << '\n'; 
  os << "reverseImageOutFileName      = " << reverseImageOutFileName      << '\n';
  
  os << "displacementFieldFileName    = " << displacementFieldFileName    << '\n';
  
  os << "outDir                       = " << outDir                       << '\n';
  
  os << "jic                          = " << jic                          << '\n';
  
  os << '\n';
  
  return os;
}


ostream& operator<<(ostream& os, const InputParams& params) {
  return params.print(os);
}


void
errorReadFile(const char* p, const char* p2 = "") {
  std::cerr << __FILE__ << p << ' ' << p2 << std::endl;
  std::exit(EXIT_FAILURE);
}


void
errorReadFile(const char* p, const string s = "") {
  std::cerr << __FILE__ << p << ' ' << s << std::endl;
  std::exit(EXIT_FAILURE);
}

void
errorReadFile(const string p, const string s = "") {
  std::cerr << __FILE__ << ":\n" << p << ' ' << s << std::endl;
  std::exit(EXIT_FAILURE);
}

void
errorReadFile(ifstream& ifs, const char* p, int c) throw(string) {
  ifs.clear();
  cerr << __FILE__ << "\nInputfile " << p << " appears to contain " 
       << c << " elements.\n\n";
  string s = "Bollocks.\n";
  //  throw(s);    // the gratuitous exception
  std::exit(EXIT_FAILURE);
}
