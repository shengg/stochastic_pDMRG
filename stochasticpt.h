#ifndef STOCHASTICPT_H
#define STOCHASTICPT_H

#include <vector>
#include "global.h"

namespace SpinAdapted{
  class StateQuantaInfo
  {
    struct leftquantainfo
    {
      int first; // The number of left quanta.
      int second;// The number of states in the quanta.
      int third; // After collecting combined quanta from left quanta and dot quanta, the position of this quanta in combined quanta.
      leftquantainfo() :first(0), second(0), third(0){};
      leftquantainfo(int a, int b, int c) :first(a), second(b), third(c){};
    
      friend class boost::serialization::access;
      template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
          ar & first \
            & second\
            & third;
        }
    
    };
    typedef std::pair<int,int> intpair;
    int currentstate;
    std::vector< std::map < intpair, intpair> > quantalink;
    std::vector< std::map < intpair, leftquantainfo> > reversequantalink;
    std::vector< std::vector<Matrix> > RotationMatrices;
    std::vector< std::vector<int> > Spins;
    public:
    double norm;
    public:
    StateQuantaInfo(int state=0): currentstate(state), quantalink(), reversequantalink(), RotationMatrices(){};
    void StoreQuantaInformation(SweepParams &sweepParams, const bool &forward);
    void readRotationandQuanta();
    double getcoeff(const std::vector<int>& ci);
    void preparestate();
    void initstate();
    void get_slater();
    double sampling(std::vector<int>& ci);
    static double local_energy(const std::vector<int>& ci);
  };
}
#endif
