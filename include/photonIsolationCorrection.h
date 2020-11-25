#ifndef PHOTONISOLATIONCORRECTION_H
#define PHOTONISOLATIONCORRECTION_H

//cpp
#include <iostream>
#include <string>
#include <vector>

//Local
//#include "include/checkMakeDir.h"

class photonIsolationCorrection
{
 public:
  photonIsolationCorrection(){return;}
  photonIsolationCorrection(std::string inTableFile);
  ~photonIsolationCorrection(){};
  void SetTable(std::string inTableFile);
  double GetCorrectedIsolation(double inVal);
  void PrintTableTex();
  
 private:
  std::string m_tableFileName;
  checkMakeDir m_check;
  bool m_isInit;
  bool m_isDescending;
  std::vector<double> m_centVals;
};

#endif
