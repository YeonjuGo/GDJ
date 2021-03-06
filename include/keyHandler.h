//Author: Chris McGinn (2020.03.04)
//Contact at chmc7718@colorado.edu or cffionn on skype for bugs

#ifndef KEYHANDLER_H
#define KEYHANDLER_H

//c+cpp
#include <string>
#include <vector>

class keyHandler{
 public:
  keyHandler(){};
  keyHandler(std::string in_handlerName);
  keyHandler(std::string in_handlerName, std::vector<unsigned long long> in_valMaxes);
  ~keyHandler();

  bool Init(std::string in_handlerName, std::vector<unsigned long long> in_valMaxes);
  bool Init(std::vector<unsigned long long> in_valMaxes);

  unsigned long long GetKey(std::vector<unsigned long long> in_vals);
  std::vector<unsigned long long> InvertKey(unsigned long long);
  void Clean();
  
 private:
  std::string m_handlerName;
  bool m_doDebug;
  std::vector<unsigned long long> m_valMaxes;
  std::vector<unsigned long long> m_multipliers;

  unsigned long long GetNearestTen(unsigned long long inVal);
};

#endif
