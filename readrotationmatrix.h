#ifndef READ_ROTATIONMATRIX_H
#define READ_ROTATIONMATRIX_H

#include <vector>
#include "global.h"
#include "StateInfo.h"
#include "initblocks.h"
#include "ObjectMatrix.h"
#include "Stackwavefunction.h"

namespace SpinAdapted{
  void ReadQuanta(std::vector<SpinQuantum>& quanta, std::vector<int>& quantaStates, std::vector<int>& newquantaStates, int dotsite);
  void ReadRotationMatrix(std::vector<Matrix>& RotationM, std::vector<int>& quantaStates, std::vector<int>& newquantaStates, int dotsite);
  //void MakeRotationMatrix(std::vector<Matrix>& rotationmatrix, StateInfo newStateInfo, int dotsite);
  void BuildFromRotationMatrix(int statea);
  void MakeWavefunction(StateInfo& sysStateInfo, StackWavefunction& wave ,int wavenum);
  //double do_one();
  //void CanonicalizeWavefunction(SweepParams &sweepParams, const bool &forward, int currentstate);
}

#endif
