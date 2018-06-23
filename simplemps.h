#ifndef NONSPINMPS_H
#define NONSPINMPS_H
#include <vector>
#include <map>
#include "global.h"
#include "timer.h"
#include "sampling.h"
#include <bitset>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/map.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

class simplemps
{
  public:
    virtual void build(int state)=0 ;//Derized class has different serialization.
    virtual void convertblockfiles()=0;
    virtual double getcoeff(const bitstring& ci)=0;
    virtual double sampling(bitstring& ci)=0;
    virtual void initstate()=0;

    simplemps(){};
    virtual ~simplemps(){};
    virtual const double get_norm() const =0;
    virtual void print_timer() =0;
};

class Abelianmps :public simplemps
{

  friend class boost::serialization::access;
  template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & physicaldim \
         & num_spatial_orbs \
         & currentstate \
         & SiteMatrices \
         & Leftindex \
         & Rightindex \
         & LefttoRight \
         & RighttoLeft \
         & norm;
    }

    int physicaldim;
    int num_spatial_orbs;
    int currentstate;
    double norm;
    double mpscoeffT=0.0;
    double mpssampleT=0.0;
    double HdiagonalT = 0.0;

  public:
    Abelianmps(){;}
    std::vector<std::vector<std::vector<Matrix> > > SiteMatrices;
    std::vector<std::vector<std::vector<int> > > Leftindex;
    std::vector<std::vector<std::vector<int> > > Rightindex;
    std::vector<std::vector<std::vector<int> > > LefttoRight;
    std::vector<std::vector<std::vector<int> > > RighttoLeft;
    void build(int state);
    void convertblockfiles();
    double getcoeff(const bitstring& ci);
    double sampling(bitstring& ci);
    void initstate();
    const double get_norm() const {return norm;}
    void print_timer(){
      cout <<"coeff time: "<<mpscoeffT<<endl;
      cout <<"sample time: "<<mpssampleT<<endl;
    }
    //void initstate();


    //int currentstate;
};

class NonAbelianmps :public simplemps
{

  friend class boost::serialization::access;
  template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & physicaldim \
         & num_spatial_orbs \
         & currentstate \
         & SiteMatrices \
         & Leftindex \
         & Rightindex \
         & LefttoRight \
         & RighttoLeft \
         & spin_vec\
         & norm;
    }

    int physicaldim;
    int num_spatial_orbs;
    int currentstate;
    double norm;
    double mpscoeffT=0.0;
    double mpssampleT=0.0;
    double HdiagonalT = 0.0;

  public:
    NonAbelianmps(){;}
    std::vector<std::vector<std::map<std::pair<int,int>, Matrix> > > SiteMatrices;
    std::vector<std::vector<std::vector<std::vector<int> > > > Leftindex;
    std::vector<std::vector<std::vector<std::vector<int> > > > Rightindex;
    std::vector<std::vector<std::vector<std::vector<int> > > > LefttoRight;
    std::vector<std::vector<std::vector<std::vector<int> > > > RighttoLeft;
    std::vector<std::vector<int> > spin_vec;//It is the 2*s of each quanta num.
    void build(int state);
    void convertblockfiles();
    double getcoeff(const bitstring& ci);
    double sampling(bitstring& ci);
    void initstate();
    const double get_norm() const {return norm;}
    void print_timer(){
      cout <<"coeff time: "<<mpscoeffT<<endl;
      cout <<"sample time: "<<mpssampleT<<endl;
    }
    //void initstate();


    //int currentstate;
};
#endif
