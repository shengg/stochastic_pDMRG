#include "MatrixBLAS.h"
#include "Stackspinblock.h"
#include "initblocks.h"
#include "input.h"
#include "rotationmat.h"
#include "Stackwavefunction.h"
#include "global.h"
#include "sweep.h"
#include <unordered_map>
#include <boost/serialization/map.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <random>
#include "MatrixBLAS.h"
#include "IntegralMatrix.h"
#include <chrono>
#include "stochasticpt.h"
using namespace SpinAdapted;
using namespace SpinAdapted::Sweep;

void ReadInput(char* conf);

void getComplementarySites(std::vector<int> &sites, std::vector<int> &complementarySites) ;


void StateQuantaInfo::StoreQuantaInformation(SweepParams &sweepParams, const bool &forward)
{
    bool restart=false, warmUp = false;
    std::vector<int> sites, spinsites;

    int new_site, wave_site;
    new_site = 0;
    sites.push_back(new_site);
    if (dmrginp.spinAdapted())
      spinsites.push_back(new_site);
    else {
      spinsites.push_back(2*new_site);
      spinsites.push_back(2*new_site+1);
    }
    StateInfo stateinfo; 
    std::vector<Matrix> RotationMatrix;

    makeStateInfo(stateinfo, new_site);
    StateInfo::store(forward, sites, stateinfo, sweepParams.current_root());
    //LoadRotationMatrix (spinsites, RotationMatrix, currentstate);

    for (; sweepParams.get_block_iter() < sweepParams.get_n_iters()+1; ) {
      new_site++;
      sites.push_back(new_site);
      if (dmrginp.spinAdapted())
        spinsites.push_back(new_site);
      else {
        spinsites.push_back(2*new_site);
        spinsites.push_back(2*new_site+1);
      }

      pout <<"Sites: ";
      for (int i : sites)
        pout << i<<" ";
      pout <<endl;

      StateInfo combinedstateinfo, siteState;
      StateInfo leftStateInfo=stateinfo;
      //TensorProduct provides a pointer to left and right StateInfo.
      //Load Stateinfo will overwrite the left StateInfo pointed.

      makeStateInfo(siteState, new_site);
      TensorProduct(leftStateInfo, siteState, combinedstateinfo, NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);

      combinedstateinfo.CollectQuanta();
      StateInfo uncollectedstateinfo = *(combinedstateinfo.unCollectedStateInfo);
      pout << uncollectedstateinfo <<endl;
      pout << combinedstateinfo<<endl;

      LoadRotationMatrix (spinsites, RotationMatrix, currentstate);

      SpinAdapted::StateInfo::transform_state(RotationMatrix, combinedstateinfo, stateinfo);
      StateInfo::store(forward, sites, stateinfo, sweepParams.current_root());
      pout << stateinfo <<endl;

      //typedef std::unordered_map < std::pair<int,int>, std::pair<int,int> > 
      typedef std::pair<int,int> intpair;
      std::map < intpair, intpair> combinedquanta;
      std::map < intpair, leftquantainfo> reversequanta;
      //Print A^{n_i} dimension.
      for (int i=0; i<combinedstateinfo.quanta.size();i++)
      {
        //Some quanta num discarded after renormalization.
        auto iter = std::find(stateinfo.quanta.begin(), stateinfo.quanta.end(), combinedstateinfo.quanta[i]);
        if (iter !=stateinfo.quanta.end() )
        {
          int index = std::distance( stateinfo.quanta.begin(), iter );
          int quantastate =0;

          for(int k=0; k < combinedstateinfo.oldToNewState[i].size(); k++)
          {
            int oldquanta = combinedstateinfo.oldToNewState[i][k];
            std::pair<int,int> combinepair (uncollectedstateinfo.leftUnMapQuanta[oldquanta],uncollectedstateinfo.rightUnMapQuanta[oldquanta]);


            //TODO
            //Discarded quanta is still in Rotation matrix.
            //std::pair<int,int> pointedpair(index,quantastate);
            std::pair<int,int> pointedpair(index,quantastate);
            combinedquanta[combinepair] = pointedpair;
            std::pair<int,int> d_r_quanta (uncollectedstateinfo.rightUnMapQuanta[oldquanta],index);
            int leftquantastates = leftStateInfo.quantaStates[uncollectedstateinfo.leftUnMapQuanta[oldquanta] ];
            leftquantainfo left_index(uncollectedstateinfo.leftUnMapQuanta[oldquanta], leftquantastates, quantastate);
            reversequanta[d_r_quanta] = left_index;
            quantastate += uncollectedstateinfo.quantaStates[oldquanta];
          }
        }
      }
      for (auto it =combinedquanta.begin(); it !=combinedquanta.end();it++)
        pout <<it->first.first <<" "<<it->first.second<<"-> " <<it->second.first<<" "<<it->second.second<<endl;


      {
        std::string file;
        file = str(boost::format("%s%s%d%s%d%s") % dmrginp.load_prefix() % "/quantalink-state" % new_site %"-"%currentstate %".tmp" );
        std::ofstream s(file.c_str());
        boost::archive::text_oarchive oa(s);
        oa << combinedquanta;
        s.close();
      }
      {
        std::string file;
        file = str(boost::format("%s%s%d%s%d%s") % dmrginp.load_prefix() % "/reversequantalink-state" % new_site %"-"%currentstate %".tmp" );
        std::ofstream s(file.c_str());
        boost::archive::text_oarchive oa(s);
        oa << reversequanta;
        s.close();
      }
    

      ++sweepParams.set_block_iter();
    }
    new_site++;

    //TODO
    //The wave function.
    std::vector<int> complementarySites;
    if (dmrginp.spinAdapted())
    {
      complementarySites.assign(1,dmrginp.last_site()-1);
    }
    else
    {
      complementarySites.push_back(dmrginp.last_site()/2-1);

    }
    //wave function and rotation matrix files are named by spin orbitals.
    //StateInfo files are named by spatial orbitals.
    getComplementarySites(spinsites,complementarySites);

    StackWavefunction lastsitewave;
    StateInfo waveinfo;
    //TODO
    lastsitewave.LoadWavefunctionInfo (waveinfo, spinsites, currentstate, true);
    norm = DotProduct(lastsitewave,lastsitewave);

    pout <<"Last site state info " <<endl;
    pout <<stateinfo<<endl;
    pout <<waveinfo<<endl;
    pout <<*(waveinfo.leftStateInfo)<<endl;
    pout <<*(waveinfo.rightStateInfo)<<endl;
    pout <<lastsitewave<<endl;

    //TODO
    //Assuming quanta num waveinfo is not collected.
    //There is one entangle bond between the last site and other sites in each quanta num.
    typedef std::pair<int,int> intpair;
    std::map < intpair, intpair> combinedquanta;
    std::map < intpair, leftquantainfo> reversequanta;

    StateInfo combinedstateinfo, siteState;
    makeStateInfo(siteState, new_site);
    TensorProduct(stateinfo, siteState, combinedstateinfo, PARTICLE_SPIN_NUMBER_CONSTRAINT);


    for (int i=0; i<combinedstateinfo.quanta.size();i++)
    {
      int quantastate =0;

      std::pair<int,int> combinepair (combinedstateinfo.leftUnMapQuanta[i],combinedstateinfo.rightUnMapQuanta[i]);

      std::pair<int,int> pointedpair(i,quantastate);
      combinedquanta[combinepair] = pointedpair;
      std::pair<int,int> d_r_quanta (combinedstateinfo.rightUnMapQuanta[i],i);
      leftquantainfo leftquanta(combinedstateinfo.leftUnMapQuanta[i],stateinfo.quantaStates[combinedstateinfo.leftUnMapQuanta[i]] , quantastate);
      reversequanta[d_r_quanta] = leftquanta;

    }
    {
    std::string file;
    file = str(boost::format("%s%s%d%s%d%s") % dmrginp.load_prefix() % "/quantalink-state" % new_site %"-"%currentstate %".tmp" );
    std::ofstream s(file.c_str());
    boost::archive::text_oarchive oa(s);
    oa << combinedquanta;
    s.close();
    }
    {
    std::string file;
    file = str(boost::format("%s%s%d%s%d%s") % dmrginp.load_prefix() % "/reversequantalink-state" % new_site %"-"%currentstate %".tmp" );
    std::ofstream s(file.c_str());
    boost::archive::text_oarchive oa(s);
    oa << reversequanta;
    s.close();
    }
    for (auto it =combinedquanta.begin(); it !=combinedquanta.end();it++)
      pout <<it->first.first <<" "<<it->first.second<<"-> " <<it->second.first<<" "<<it->second.second<<endl;





    std::vector<Matrix> lastRotationMatrix;
    for (int a=0; a<lastsitewave.nrows(); a++)
      for (int b=0; b<lastsitewave.ncols(); b++)
      {

        if (lastsitewave.allowed(a, b)) {
	        StackMatrix& lM = lastsitewave.operator_element(a, b);
	        Matrix tM = RotationMatrix[a];
          if (RotationMatrix[a].Ncols()==0) continue;
	        Matrix nM(1,1);
          nM(1,1) = 0.0;
	        MatrixMultiply(tM, 't', lM, 'n', nM, 1.0);
          lastRotationMatrix.push_back(nM);
          pout << nM<<endl;
        }
      }

      sites.push_back(new_site);
      if (dmrginp.spinAdapted())
        spinsites.push_back(new_site);
      else {
        spinsites.push_back(2*new_site);
        spinsites.push_back(2*new_site+1);
      }
      SaveRotationMatrix (spinsites, lastRotationMatrix, currentstate);



}

void StateQuantaInfo::readRotationandQuanta()
{
quantalink.resize(0);
reversequantalink.resize(0);
RotationMatrices.resize(0);

std::vector<Matrix> RotationMatrix;
std::vector<Matrix> reducedRotationMatrix;//Remove the discarded quanta num.

int niters;
int new_site = 0;
std::vector<int> spinsites, sites;
sites.push_back(new_site);
if (dmrginp.spinAdapted())
{

  niters = dmrginp.last_site()-1;
  spinsites.push_back(new_site);
}
else {
  niters = dmrginp.last_site()/2-1;
  spinsites.push_back(2*new_site);
  spinsites.push_back(2*new_site+1);
}

StateInfo stateinfo;
StateInfo::restore(true, sites, stateinfo, currentstate);
std::vector<int> spin;
for(SpinQuantum i : stateinfo.quanta)
  spin.push_back(i.totalSpin.getirrep());
Spins.push_back(spin);

if (mpigetrank() == 0) {
  for(int i=0; i<niters; i++)
  {
    new_site++;
    sites.push_back(new_site);
    if (dmrginp.spinAdapted())
    {
      spinsites.push_back(new_site);
    }
    else {
      spinsites.push_back(2*new_site);
      spinsites.push_back(2*new_site+1);
    }
    
    reducedRotationMatrix.clear();
    LoadRotationMatrix (spinsites, RotationMatrix, currentstate);
    pout <<"total size: " <<RotationMatrix.size()<<endl;
    for (Matrix m : RotationMatrix)
    {
      if(m.Ncols()) reducedRotationMatrix.push_back(m);
    }
    RotationMatrices.push_back(reducedRotationMatrix);
  
    StateInfo stateinfo;
    StateInfo::restore(true, sites, stateinfo, currentstate);
    std::vector<int> spin;
    for(SpinQuantum i : stateinfo.quanta)
      spin.push_back(i.totalSpin.getirrep());
    Spins.push_back(spin);

    {
      std::string file;
      file = str(boost::format("%s%s%d%s%d%s") % dmrginp.load_prefix() % "/quantalink-state" % new_site %"-"%currentstate %".tmp" );
      std::ifstream s(file.c_str());
      boost::archive::text_iarchive oa(s);
      std::map < intpair, intpair> combinedquanta;
      oa >> combinedquanta;
      s.close();
      quantalink.push_back(combinedquanta);
    }
    {
      std::string file;
      file = str(boost::format("%s%s%d%s%d%s") % dmrginp.load_prefix() % "/reversequantalink-state" % new_site %"-"%currentstate %".tmp" );
      std::ifstream s(file.c_str());
      boost::archive::text_iarchive oa(s);
      std::map < intpair, leftquantainfo> reversequanta;
      oa >> reversequanta;
      s.close();
      reversequantalink.push_back(reversequanta);
    }
  }
}
#ifndef SERIAL
      mpi::communicator world;
      mpi::broadcast(calc, quantalink, 0);
      mpi::broadcast(calc, reversequantalink, 0);
      mpi::broadcast(calc, RotationMatrices, 0);
#endif
/*
  new_site++;
  if (dmrginp.spinAdapted())
  {
    spinsites.push_back(new_site);
  }
  else {
    spinsites.push_back(2*new_site);
    spinsites.push_back(2*new_site+1);
  }
  LoadRotationMatrix (spinsites, RotationMatrix, currentstate);
  RotationMatrices.push_back(RotationMatrix);
  */

}

void StateQuantaInfo::get_slater()
{
  std::vector<int> ci;
  ci.assign(dmrginp.last_site(),0);
  double norm = 0.0;
  
  for(int n=0;n<pow(2,ci.size());n++)
  {
    for(int i=0;i<ci.size();i++)
      ci[i] = (n & ( 1 << i )) >> i;
    int particle_num = 0;
    for(int i=0;i<ci.size();i++)
      particle_num+= ci[i];
    if (particle_num!=dmrginp.total_particle_number()) continue;

    double coeff = getcoeff(ci);
    norm +=coeff*coeff;
    //if (abs(coeff)>1e-4)
    {
    for(int i=0;i<ci.size();i++)
      pout << ci[i];
    pout <<": ";
    pout << coeff;
    pout <<endl;
    }
  }
  cout <<"END"<<endl;

}

double StateQuantaInfo::getcoeff(const std::vector<int>& ci)
{
    double coeff = 1.0;
    int quantanum;
    int statenum=0;

    int niters;
    if (dmrginp.spinAdapted())
    {
    
      niters = dmrginp.last_site()-2;
    }
    else {
      niters = dmrginp.last_site()/2-2;
    }

    if (dmrginp.spinAdapted())
      quantanum = ci[0];
    else {
      quantanum = 2*ci[0]+ci[1];
    }

    //TODO
    //call DGEMV rather than DGEMM
    Matrix lbasis(1,1);
    lbasis(1,1) = 1.0;
    for(int i=0; i< niters+1;i++)
    {
      int localquantanum = dmrginp.spinAdapted()? ci[i+1]: 2*ci[2*i+2]+ci[2*i+3];
      auto iter = quantalink[i].find(std::make_pair(quantanum,localquantanum));
      if(iter==quantalink[i].end()) 
      {
        //pout << "Quanta: "<< quantanum<<" "<<localquantanum<<" was not found"<<endl;
        //pout << "This slater determinant avoids the wave function symmetry"<<endl;
        coeff = 0.0;
        return coeff;
      }
      auto pos = iter->second;
      //pout <<"Row position" << pos.first<<" " <<pos.second<<endl;

      Matrix rotationmatrix;
      rotationmatrix.ReSize(lbasis.Ncols(),RotationMatrices[i][pos.first].Ncols());
      //pout <<"total size: " <<RotationMatrices[i].size()<<endl;
      //pout <<"Size: " << rotationmatrix.Nrows() << " " << rotationmatrix.Ncols()<<endl;
      //pout <<"Size: " << RotationMatrices[i][pos.first].Nrows() << " " << RotationMatrices[i][pos.first].Ncols()<<endl;
      for(int k=1;k<=rotationmatrix.Nrows();k++)
      for(int l=1;l<=rotationmatrix.Ncols();l++)
      {
        rotationmatrix(k,l) = RotationMatrices[i][pos.first]((pos.second)+k,l);
      }

      Matrix newlbasis;
      newlbasis.ReSize(lbasis.Nrows(), rotationmatrix.Ncols());
      newlbasis= 0.0;
      MatrixMultiply(lbasis,'n', rotationmatrix,'n',newlbasis,1.0);
      //MatrixMultiply(lbasis,'n', rotationmatrix,'n',newlbasis,1.0, 0.0);

      //pout <<newlbasis<<endl;
      quantanum = (pos.first);
      lbasis = newlbasis;
      //pout <<lbasis<<endl;

    }
    assert(lbasis.Ncols()==1);
    assert(lbasis.Nrows()==1);
    //pout <<lbasis<<endl;
    coeff = lbasis(1,1);
    return coeff;


}

void StateQuantaInfo::initstate()
{
    bool direction=false;
    SweepParams sweepParams;
    sweepParams.current_root() = currentstate;
    Sweep::InitializeStateInfo(sweepParams, !direction, currentstate);
    Sweep::InitializeStateInfo(sweepParams, direction, currentstate);
    Sweep::CanonicalizeWavefunction(sweepParams, !direction, currentstate);
    Sweep::CanonicalizeWavefunction(sweepParams, direction, currentstate);
    Sweep::CanonicalizeWavefunction(sweepParams, !direction, currentstate);
    Sweep::InitializeStateInfo(sweepParams, !direction, currentstate);
    Sweep::InitializeStateInfo(sweepParams, direction, currentstate);
}

void StateQuantaInfo::preparestate()
{
  initstate();
  SweepParams sweepParams;
  sweepParams.set_sweep_iter() = 0;
  sweepParams.current_root() = currentstate;
  if (mpigetrank() == 0) {
    StoreQuantaInformation(sweepParams, true);
  }
  readRotationandQuanta();
}

double StateQuantaInfo::sampling(std::vector<int>& ci)
{
    int statenum=0;

    int niters;
    if (dmrginp.spinAdapted())
    {
    
      niters = dmrginp.last_site()-1;
    }
    else {
      niters = dmrginp.last_site()/2-1;
    }
    std::vector<int> reverseci;

    std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> f(0.0,1.0);



    ColumnVector oldwave;
    int rightquanta;
    for(int i=niters-1; i>=0 ;i--)
    {
      int localquanta =0;
      std::vector<double> coeff2(4,0.0);
      std::vector<ColumnVector> newwave(4);
      for(localquanta=0;localquanta<4;localquanta++)
      {
        int leftquanta;

        if(i==niters-1)
        {
          for(auto quantapair:reversequantalink[i])
          {
            if(quantapair.first.first == localquanta)
            {
              oldwave = RotationMatrices[i][quantapair.first.second];
              auto quanta_l = reversequantalink[i].at(quantapair.first);
              coeff2[localquanta] = oldwave(1)*oldwave(1);
              break;
            }
          }
        }
        else
        {
          auto quanta_d_r = std::make_pair(localquanta, rightquanta);
          auto iter = reversequantalink[i].find(quanta_d_r);
          if (iter==reversequantalink[i].end())
          {
            coeff2[localquanta] = 0.0;
          }
          else{

            auto quanta_l = iter->second;
            Matrix& rotateM1 = RotationMatrices[i][quanta_d_r.second];
            Matrix rotateM2;
            rotateM2.ReSize(quanta_l.second, rotateM1.Ncols());
            for(int k=1;k<=rotateM2.Nrows();k++)
            for(int l=1;l<=rotateM2.Ncols();l++)
              rotateM2(k,l) = rotateM1(k+quanta_l.third,l);
            newwave[localquanta].ReSize(rotateM2.Nrows());
            newwave[localquanta] = 0.0;multiplyH_2indexmultiplyH_2indexmultiplyH_2index
            MatrixMultiply(rotateM2,'n', oldwave,'n',newwave[localquanta], 1.0, 0.0);
            coeff2[localquanta] = dotproduct(newwave[localquanta],newwave[localquanta]);
          }
        }

      }
      std::vector<double> sum(4);
      std::partial_sum (coeff2.begin(),coeff2.end() , sum.begin());
      double x = f(generator);
      //cout <<"X: "<<x <<endl;
      //cout <<"coeff: "<<coeff2[0]<< " " <<coeff2[1]<< " "<<coeff2[2]<< " "<<coeff2[3]<< " "<<endl;
      if(x<sum[0]/sum[3])
      {
        localquanta = 0;
      }
      else if(x<sum[1]/sum[3])
      {
        localquanta = 1;
      }
      else if(x<sum[2]/sum[3])
      {
        localquanta = 2;
      }
      else
      {
        localquanta = 3;
      }
      if (i==niters-1)
      {
        for(auto quantapair:reversequantalink[i])
        {
          if(quantapair.first.first == localquanta)
          {
            rightquanta = quantapair.second.first;
            oldwave = RotationMatrices[i][quantapair.first.second];
            break;
          }
        }
        reverseci.push_back(localquanta%2);
        reverseci.push_back(localquanta/2);
      }
      else
      {
        auto quanta_d_r = std::make_pair(localquanta, rightquanta);
        auto quanta_l = reversequantalink[i].at(quanta_d_r);
        reverseci.push_back(localquanta%2);
        reverseci.push_back(localquanta/2);
        if (i==0)
        {
          reverseci.push_back(quanta_l.first%2);
          reverseci.push_back(quanta_l.first/2);
        }
        rightquanta = quanta_l.first;
        oldwave = newwave[localquanta];
      }

    }
    ci.resize(reverseci.size());
    for(int i=0;i<reverseci.size();i++)
      ci[i] = reverseci[reverseci.size()-1-i];
    assert(oldwave.Nrows()==1);
    return oldwave(1);



}

double StateQuantaInfo::sampling_spinadpated(std::vector<int>& ci)
{
    int statenum=0;

    int niters;
    if (dmrginp.spinAdapted())
    {
    
      niters = dmrginp.last_site()-1;
    }
    else {
      niters = dmrginp.last_site()/2-1;
    }
    std::vector<int> reverseci;

    std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> f(0.0,1.0);



    ColumnVector oldwave;
    int rightquanta;
    for(int i=niters-1; i>=0 ;i--)
    {
      int localquanta =0;
      std::vector<double> coeff2(4,0.0);
      std::vector<ColumnVector> newwave(4);
      for(localquanta=0;localquanta<4;localquanta++)
      {
        int leftquanta;

        if(i==niters-1)
        {
          for(auto quantapair:reversequantalink[i])
          {
            if(quantapair.first.first == localquanta)
            {
              oldwave = RotationMatrices[i][quantapair.first.second];
              auto quanta_l = reversequantalink[i].at(quantapair.first);
              coeff2[localquanta] = oldwave(1)*oldwave(1);
              break;
            }
          }
        }
        else
        {
          auto quanta_d_r = std::make_pair(localquanta, rightquanta);
          auto iter = reversequantalink[i].find(quanta_d_r);
          if (iter==reversequantalink[i].end())
          {
            coeff2[localquanta] = 0.0;
          }
          else{

            auto quanta_l = iter->second;
            Matrix& rotateM1 = RotationMatrices[i][quanta_d_r.second];
            Matrix rotateM2;
            rotateM2.ReSize(quanta_l.second, rotateM1.Ncols());
            for(int k=1;k<=rotateM2.Nrows();k++)
            for(int l=1;l<=rotateM2.Ncols();l++)
              rotateM2(k,l) = rotateM1(k+quanta_l.third,l);
            newwave[localquanta].ReSize(rotateM2.Nrows());
            newwave[localquanta] = 0.0;multiplyH_2indexmultiplyH_2indexmultiplyH_2index
            MatrixMultiply(rotateM2,'n', oldwave,'n',newwave[localquanta], 1.0, 0.0);
            coeff2[localquanta] = dotproduct(newwave[localquanta],newwave[localquanta]);
          }
        }

      }
      std::vector<double> sum(4);
      std::partial_sum (coeff2.begin(),coeff2.end() , sum.begin());
      double x = f(generator);
      //cout <<"X: "<<x <<endl;
      //cout <<"coeff: "<<coeff2[0]<< " " <<coeff2[1]<< " "<<coeff2[2]<< " "<<coeff2[3]<< " "<<endl;
      if(x<sum[0]/sum[3])
      {
        localquanta = 0;
      }
      else if(x<sum[1]/sum[3])
      {
        localquanta = 1;
      }
      else if(x<sum[2]/sum[3])
      {
        localquanta = 2;
      }
      else
      {
        localquanta = 3;
      }
      if (i==niters-1)
      {
        for(auto quantapair:reversequantalink[i])
        {
          if(quantapair.first.first == localquanta)
          {
            rightquanta = quantapair.second.first;
            oldwave = RotationMatrices[i][quantapair.first.second];
            break;
          }
        }
        reverseci.push_back(localquanta%2);
        reverseci.push_back(localquanta/2);
      }
      else
      {
        auto quanta_d_r = std::make_pair(localquanta, rightquanta);
        auto quanta_l = reversequantalink[i].at(quanta_d_r);
        reverseci.push_back(localquanta%2);
        reverseci.push_back(localquanta/2);
        if (i==0)
        {
          reverseci.push_back(quanta_l.first%2);
          reverseci.push_back(quanta_l.first/2);
        }
        rightquanta = quanta_l.first;
        oldwave = newwave[localquanta];
      }

    }
    ci.resize(reverseci.size());
    for(int i=0;i<reverseci.size();i++)
      ci[i] = reverseci[reverseci.size()-1-i];
    assert(oldwave.Nrows()==1);
    return oldwave(1);



}

double StateQuantaInfo::local_energy(const std::vector<int>& ci)
{
  int integralIndex =0;
  double energy = 0;
  energy += coreEnergy[integralIndex];
  for(int i=0;i<ci.size();i++)
    energy += ci[i]*v_1[integralIndex](i,i);
  for(int i=0;i<ci.size();i++)
  for(int j=i+1;j<ci.size();j++)
    energy += ci[i]*ci[j]*v_2[integralIndex](i,j,i,j);
  for(int i=0;i<ci.size();i++)
  for(int j=i+1;j<ci.size();j++)
    energy -= ci[i]*ci[j]*v_2[integralIndex](j,i,i,j);
  //for(int i=0;i<ci.size()/2;i++)
  //  energy -= ci[2*i]*ci[2*i+1]*v_2[integralIndex](2*i,2*i+1,2*i,2*i+1);
  return energy;


}

int main(int argc, char* argv[])
{
  //test(argc,argv);
//  for(auto i: argv)
//    cout <<string(i)<<endl;
#ifndef SERIAL
  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world;
  //MPI_Comm Calc;
  if (world.rank() == 0) {
    pout << "Runing with " << world.size() << " processors" << endl;
  }
#endif

//  for(auto i: argv)
//    cout <<string(i)<<endl;

  ReadInput(argv[1]);
  dmrginp.matmultFlops.resize(numthrds, 0.);

  double* stackmemory = new double[dmrginp.getMemory()];
  Stackmem.resize(numthrds);
  Stackmem[0].data = stackmemory;
  Stackmem[0].size = dmrginp.getMemory();

  dmrginp.initCumulTimer();
  StateQuantaInfo baseState(0);
  baseState.preparestate();
  StateQuantaInfo pertuberState(1000);
  pertuberState.preparestate();
  //get_slater();
  std::vector<int> ci;
  const unsigned long num_sample = atoi(argv[2]);
  pout << "Num of samples: "<< num_sample<<endl;
  //const unsigned long num_sample = 20000;
  double base_H0=0.0;
  double base_H01=0.0;
  double test = 0.0;
  for (int i=0;i<num_sample;i++)
  {
    double coeff = baseState.sampling(ci);
    if (abs(StateQuantaInfo::local_energy(ci)) < 1e-13) cout <<"ERROR" <<"Too small E_{ci}"<<endl;
    base_H0+= 1/(StateQuantaInfo::local_energy(ci)*num_sample);
    test += StateQuantaInfo::local_energy(ci)/num_sample;

    base_H01 += pertuberState.getcoeff(ci)/(StateQuantaInfo::local_energy(ci)*num_sample*coeff);
    //for(auto i:ci)
    //  cout <<i;
    //cout <<":" << coeff<<endl;
  }

  //cout <<world.rank()<<" "<<"base_H0: " << base_H0<<endl;


  double perturber_H0=0.0;
  double perturber_H01=0.0;
  for (int i=0;i<num_sample;i++)
  {
    double coeff = pertuberState.sampling(ci);
    perturber_H0+= 1/(StateQuantaInfo::local_energy(ci)*num_sample);
    perturber_H01 += baseState.getcoeff(ci)/(StateQuantaInfo::local_energy(ci)*num_sample*coeff);
  }


#ifndef SERIAL
  double global_base_H0;
  double global_base_H01;
  double global_perturber_H0;
  double global_perturber_H01;
  reduce(world, base_H0, global_base_H0, std::plus<double>(),0);
  reduce(world, base_H01, global_base_H01, std::plus<double>(),0);
  reduce(world, perturber_H0, global_perturber_H0, std::plus<double>(),0);
  reduce(world, perturber_H01, global_perturber_H01, std::plus<double>(),0);
  base_H0       =  global_base_H0;
  base_H01      =  global_base_H01;
  perturber_H0  =  global_perturber_H0;
  perturber_H01 =  global_perturber_H01;
#endif
//  if(world.rank())
//  {
//  MPI_Reduce(&base_H0      , &test, 1, MPI_DOUBLE, MPI_SUM,0, Calc);
//  MPI_Reduce(&base_H01     , &base_H01, 1, MPI_DOUBLE, MPI_SUM, 0, Calc);
//  MPI_Reduce(&perturber_H0 , &perturber_H0, 1, MPI_DOUBLE, MPI_SUM, 0, Calc);
//  MPI_Reduce(&perturber_H01, &perturber_H01, 1, MPI_DOUBLE, MPI_SUM, 0, Calc);
//  }
//  else
//  {
//  MPI_Reduce(&base_H0      , &test, 1, MPI_DOUBLE, MPI_SUM,0, Calc);
//  //MPI_Allreduce(MPI_IN_PLACE, &base_H0, 1, MPI_DOUBLE, MPI_SUM, Calc);
//  //MPI_Reduce(MPI_IN_PLACE, &base_H0, 1, MPI_DOUBLE, MPI_SUM, 0, Calc);
//  MPI_Reduce(MPI_IN_PLACE, &base_H01, 1, MPI_DOUBLE, MPI_SUM, 0, Calc);
//  MPI_Reduce(MPI_IN_PLACE, &perturber_H0, 1, MPI_DOUBLE, MPI_SUM, 0, Calc);
//  MPI_Reduce(MPI_IN_PLACE, &perturber_H01, 1, MPI_DOUBLE, MPI_SUM, 0, Calc);
//  cout <<test<<"TEST"<<endl;
//
//  }
//  cout <<world.rank()<<" "<<"base_H0: " << base_H0<<endl;
  if(!world.rank())
  {
    base_H0 /= world.size();
    base_H01 /= world.size();
    perturber_H0 /= world.size();
    perturber_H01 /= world.size();
    perturber_H01 *=pertuberState.norm;
    perturber_H0 *=pertuberState.norm;
    cout <<"base_H0: " << base_H0<<endl;
    cout <<"base_H01: " << base_H01<<endl;
    cout <<"perturber_norm: " << pertuberState.norm<<endl;
    cout <<"TEST: " <<test <<endl;
    cout <<"perturber_H0: " << perturber_H0<<endl;
    cout <<"perturber_H01: " << perturber_H01<<endl;

    double pt2_energy = -1.0*perturber_H0+perturber_H01*perturber_H01/(base_H0);
    cout <<"perturbation energy: " << pt2_energy<<endl;
    pt2_energy = -1.0*perturber_H0+base_H01*base_H01/(base_H0);
    cout <<"perturbation energy: " << pt2_energy<<endl;
  }



  return 0;
}

