#ifndef BINUTILS_H
#define BINUTILS_H

//c+cpp
#include <cmath>
#include <iostream>
#include <string>

//Local
#include "include/stringUtil.h"

#define _USE_MATH_DEFINES

inline bool goodBinning(std::string config, int nMaxBins, int nBins, std::string binType)
{
  if(nBins > nMaxBins){
    std::cout << "BINUTILS ERROR IN GOODBINNING - config \'" << config << "\' contains " << binType << "  of size \'" << nBins << "', greater than max val \'" << nMaxBins << "\'. Please fix. return false" << std::endl;
    return false;
  }
  return true;
}

inline double mathStringToNum(std::string inStr, const bool doGlobalDebug=false)//currently just handles fractions and pi
{
  double outVal = 1.0;
  inStr = returnAllCapsString(inStr);
  if(inStr.find("PI") != std::string::npos){
    outVal *= M_PI;
    inStr.replace(inStr.find("PI"), 2, "");
  }
  if(inStr.find("/") != std::string::npos){
    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << inStr << std::endl;

    double num = 1.0;
    if(inStr.find("/") != 0) num = std::stod(inStr.substr(0, inStr.find("/")));

    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << inStr << std::endl;

    double denom = std::stod(inStr.substr(inStr.find("/")+1, inStr.size()));
    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    outVal *= num/denom;
  }
  else outVal *= std::stod(inStr);
  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    
  return outVal;
}

inline std::string mathStringToNameStr(std::string inStr, const bool doGlobalDebug=false)//currently just handles fractions and pi
{
  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  inStr = returnAllCapsString(inStr);
  while(inStr.find("/") != std::string::npos){
    inStr.replace(inStr.find("/"), 1, "Over");
  }
  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  while(inStr.find("*") != std::string::npos){
    inStr.replace(inStr.find("*"), 1, "");
  }
    
  return inStr;
}

#endif
