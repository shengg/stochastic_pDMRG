#ifndef STOCHASTICPT_H
#define STOCHASTICPT_H

#include <vector>
#include "global.h"
#include "timer.h"
#include "sampling.h"
#include <bitset>

double local_energy(const bitstring& ci, int integralIndex)
{
  std::clock_t startTime = std::clock();

  double energy = coreEnergy[integralIndex];
  int n = ci.size();

  std::vector<int> occupy;
  for(int i=0;i<n;i++)
    if(ci[i])
      occupy.push_back(i);
  for(auto i: occupy)
    energy += v_1[integralIndex](i, i);

  for(int i=0;i<occupy.size();i++)
  for(int j=i+1;j<occupy.size();j++)
  {
    energy += v_2[integralIndex](occupy[i],occupy[j],occupy[i],occupy[j]);
    energy -= v_2[integralIndex](occupy[j],occupy[i],occupy[i],occupy[j]);
  }
  //HdiagonalT += (std::clock() -startTime)/ (double) CLOCKS_PER_SEC;
  return energy;

}

#endif
