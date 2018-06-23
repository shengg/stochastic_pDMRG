/*
 * =====================================================================================
 *
 *       Filename:  readrotatationmatrix.C
 *
 *    Description: Read rotation matrix from external DMRG code to build spinblocks.
 *
 *        Version:  1.0
 *        Created:  02/13/2017 12:29:47 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Sheng Guo 
 *
 * =====================================================================================
 */


#include "MatrixBLAS.h"
#include "Stackspinblock.h"
#include "initblocks.h"
#include "input.h"
#include "rotationmat.h"
#include "global.h"
#ifndef SERIAL
#include "mpi.h"
#include <boost/mpi.hpp>
#include "readrotationmatrix.h"
#include "sweep.h"
#include "Stackdensity.h"
#include "pario.h"
#include "Stackwavefunction.h"
#include "stackguess_wavefunction.h"
#include "operatorfunctions.h"
#endif

void ReadInput(char* conf);
void restart(double sweep_tol, bool reset_iter);

using namespace SpinAdapted;
using namespace SpinAdapted::Sweep;

void getComplementarySites(std::vector<int> &sites, std::vector<int> &complementarySites) ;

int test(int argc, char* argv[])
{
  std::vector<SpinQuantum> quanta;
  std::vector<int> quantaStates;
  std::vector<int> newquantaStates;
  std::vector<double> oldmatrix;
  for(int dotsite = 0; dotsite < 10; dotsite++){
    cout << "Site : " << dotsite <<endl;
    ReadQuanta(quanta, quantaStates, newquantaStates, dotsite+1);
    cout << "Quanta " <<endl;
    int states = 0;
    for(int i=0; i< quanta.size(); ++i)
    {
      states += quantaStates[i];

      cout << quanta[i] <<" " <<quantaStates[i] << " " <<newquantaStates[i]<<endl;
    }

    std::vector<Matrix> RotationM;
    ReadRotationMatrix(RotationM, quantaStates, newquantaStates, dotsite);
    cout << " rotation matrix " <<endl;
    for(int i=0; i< RotationM.size(); ++i)
    {

      cout << RotationM[i]<<endl;
    }
  }
  return 0;

}

int main(int argc, char* argv[])
{
  //test(argc,argv);
#ifndef SERIAL
  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world;
  if (world.rank() == 0) {
    pout << "Runing with " << world.size() << " processors" << endl;
  }
#endif

  ReadInput(argv[1]);
  dmrginp.matmultFlops.resize(numthrds, 0.);

  double* stackmemory = new double[dmrginp.getMemory()];
  Stackmem.resize(numthrds);
  Stackmem[0].data = stackmemory;
  Stackmem[0].size = dmrginp.getMemory();

  dmrginp.initCumulTimer();


  //BuildBlock(-1);
  BuildFromRotationMatrix(0);
  BuildFromRotationMatrix(-1);
  SweepParams sweepParams;
  sweepParams.set_sweep_iter() = 0;
  sweepParams.current_root() = -1;
	Sweep::CanonicalizeWavefunction(sweepParams, false, 0);
	Sweep::CanonicalizeWavefunction(sweepParams, false, -1);
	Sweep::CanonicalizeWavefunction(sweepParams, true, 0);
	Sweep::CanonicalizeWavefunction(sweepParams, true, -1);
	Sweep::CanonicalizeWavefunction(sweepParams, false, 0);
	Sweep::CanonicalizeWavefunction(sweepParams, false, -1);
  sweepParams.savestate(true,1);
  //do_one();
	//restart(1e-7, false);
  return 0;
}

  void SpinAdapted::BuildFromRotationMatrix(int statea) {

    SweepParams sweepParams;
    sweepParams.set_sweep_parameters();
    sweepParams.set_block_iter() = 0;


    bool forward = true, restart=false, warmUp = false;
    std::vector<int> sites;

    int new_site, wave_site;
    new_site = 0;
    sites.push_back(new_site);
    StateInfo stateInfo1; 
    makeStateInfo(stateInfo1, new_site);


//{
//  std::vector< std::vector<Csf> > ladders;
//  std::vector< Csf > dets; 
//  dets = CSFUTIL::spinfockstrings(sites, ladders);
//  StateInfo s = StateInfo(dets);
//  cout <<"Determinants" <<endl;
//  for(auto i: dets)
//    cout <<i <<endl;
//
//  if (dmrginp.add_noninteracting_orbs() && dmrginp.molecule_quantum().get_s().getirrep() != 0 && dmrginp.spinAdapted() && new_site == 0)
//  {
//    SpinQuantum sq = dmrginp.molecule_quantum();
//    sq = SpinQuantum(sq.get_s().getirrep(), sq.get_s(), IrrepSpace(0));
//    int qs = 1, ns = 1;
//    StateInfo addstate(ns, &sq, &qs), newstate; 
//
//    TensorProduct(addstate, s, stateInfo1, NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);
//  }
//  else
//  {
//    stateInfo1 = s;
//  }
//
//}
//






    
    {

    std::vector<SpinQuantum> quanta;
    std::vector<int> quantaStates, newquantaStates, renormquantaStates;
    renormquantaStates.resize(stateInfo1.quantaStates.size());
    ReadQuanta(quanta, quantaStates, newquantaStates, new_site+1);

      for(int i=0; i< quanta.size(); ++i)
      {
        auto it = find(stateInfo1.quanta.begin(), stateInfo1.quanta.end(), quanta[i]);
        assert(it!=stateInfo1.quanta.end());
        assert(quantaStates[i] ==stateInfo1.quantaStates[it-stateInfo1.quanta.begin()] );
      }

      for(int i=0; i< stateInfo1.quanta.size(); ++i)
      {
        auto it = find(quanta.begin(), quanta.end(), stateInfo1.quanta[i]);
        if(it==quanta.end()){
          renormquantaStates[i] = 0;
        }
        else{
          int j = it-quanta.begin();
          renormquantaStates[i] = newquantaStates[j];
        }
      }
      std::vector<Matrix> rotation1;
      ReadRotationMatrix(rotation1, stateInfo1.quantaStates, renormquantaStates, 0);
      StateInfo renormState1;
      SpinAdapted::StateInfo::transform_state(rotation1, stateInfo1, renormState1);
      stateInfo1 = renormState1;


    }




    StateInfo::store(true, sites, stateInfo1, statea);



    for (; sweepParams.get_block_iter() < sweepParams.get_n_iters()-1; ) {
      new_site++;
      wave_site = new_site+1;
      std::vector<int> complementarySites, spindotsites(1, new_site), oldsites = sites, oldcomplement;
      sites.push_back(new_site);

      StateInfo siteState, newState1, newStateInfo, rotatednewStateInfo; 
      makeStateInfo(siteState, new_site);
      TensorProduct(stateInfo1, siteState, newState1, NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);
      newState1.CollectQuanta();


      std::vector<SpinQuantum> quanta;
      std::vector<int> quantaStates, newquantaStates, renormquantaStates;
      renormquantaStates.resize(newState1.quantaStates.size());
      ReadQuanta(quanta, quantaStates, newquantaStates, new_site+1);
      for(int i=0; i< quanta.size(); ++i)
      {
        auto it = find(newState1.quanta.begin(), newState1.quanta.end(), quanta[i]);
        assert(it!=newState1.quanta.end());
        assert(quantaStates[i] ==newState1.quantaStates[it-newState1.quanta.begin()] );
      }

      for(int i=0; i< newState1.quanta.size(); ++i)
      {
        auto it = find(quanta.begin(), quanta.end(), newState1.quanta[i]);
        if(it==quanta.end()){
          renormquantaStates[i] = 0;
        }
        else{
          int j = it-quanta.begin();
          renormquantaStates[i] = newquantaStates[j];
        }
      }

      std::vector<Matrix> rotation1;
      ReadRotationMatrix(rotation1, newState1.quantaStates, renormquantaStates, new_site);
      SaveRotationMatrix (sites, rotation1, statea);

      StateInfo renormState1;
      SpinAdapted::StateInfo::transform_state(rotation1, newState1, renormState1);
      //FIXME
      //assert(newquantaStates==renormState1.quantaStates);

      StateInfo::store(forward, sites, renormState1, statea);
      stateInfo1 = renormState1;
      ++sweepParams.set_block_iter();
    }
    StackWavefunction wave;
    MakeWavefunction(stateInfo1, wave, statea);

  }

  void SpinAdapted::MakeWavefunction(StateInfo& leftleftStateInfo, StackWavefunction& wave ,int wavenum){

    //Build Wavefunction at end of a sweep.
    //Always using reverse sweep, because of the singlet embedding at the first site.
    StateInfo lastsiteState, secondsiteState, sysStateInfo, reducedsysStateInfo, renormStateInfo;
    makeStateInfo(secondsiteState, dmrginp.last_site()-2);
    TensorProduct(leftleftStateInfo, secondsiteState, sysStateInfo, NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);
    sysStateInfo.CollectQuanta();

    std::vector<int> quantaStates, newquantaStates, renormquantaStates;
    renormquantaStates.resize(sysStateInfo.quantaStates.size());
    std::vector<SpinQuantum> quanta;
    ReadQuanta(quanta, quantaStates, newquantaStates, dmrginp.last_site()-1);
    std::vector<Matrix> LeftRotation;


      for(int i=0; i< quanta.size(); ++i)
      {
        auto it = find(sysStateInfo.quanta.begin(), sysStateInfo.quanta.end(), quanta[i]);
        assert(it != sysStateInfo.quanta.end());
        assert(quantaStates[i] ==sysStateInfo.quantaStates[it-sysStateInfo.quanta.begin()] );
      }
      for(int i=0; i< sysStateInfo.quanta.size(); ++i)
      {
        auto it = find(quanta.begin(), quanta.end(), sysStateInfo.quanta[i]);
        if(it==quanta.end()){
          renormquantaStates[i] = 0;
        }
        else{
          int j = it-quanta.begin();
          renormquantaStates[i] = newquantaStates[j];
        }
      }

    ReadRotationMatrix(LeftRotation, sysStateInfo.quantaStates, renormquantaStates, dmrginp.last_site()-2);
    SpinAdapted::StateInfo::transform_state(LeftRotation, sysStateInfo, renormStateInfo);


    makeStateInfo(lastsiteState, dmrginp.last_site()-1);

    StateInfo oldcombinedState;//The sys and dot (the last site) combined together.
    TensorProduct(renormStateInfo, lastsiteState , oldcombinedState, PARTICLE_SPIN_NUMBER_CONSTRAINT);

    std::vector<int> lastquantaStates(oldcombinedState.quanta.size(),1);
    //std::vector<int> lastquantaStates(renormquantaStates.size(),1);
    std::vector<Matrix> wavedata;
    ///ReadQuanta(lastquanta, lastquantaStates, lastnewquantaStates, dmrginp.last_site());
    // only one quanta in last rotation matrix. It is just a vector (a wave fucntion).
    //assert(lastquanta.size()==1); 
    //assert(lastnewquantaStates[0] ==1);
    ReadRotationMatrix(wavedata, oldcombinedState.quantaStates, lastquantaStates, dmrginp.last_site()-1);




    //StackWavefunction wave;
    wave.initialise(dmrginp.effective_molecule_quantum_vec(), sysStateInfo, lastsiteState, true);  
    wave.Clear();

    for(int i=0; i< sysStateInfo.quanta.size();++i)
    {
      for (int j = 0; j < lastsiteState.quanta.size(); ++j) {
        if(wave.allowed(i, j) && LeftRotation[i].Ncols())
        {
          for(int k =0; k < oldcombinedState.quanta.size(); ++k)
            if(sysStateInfo.quanta[i] == renormStateInfo.quanta[oldcombinedState.leftUnMapQuanta[k]])
            {
	            StackMatrix& a = wave.operator_element(i,j);
	            Matrix& b = wavedata[k];
	            MatrixMultiply (LeftRotation[i], 'n', b, 'n', a, 1.); 
              continue;
            }
        }
      }
    }

    std::vector<int> sites;
    for(int i=0; i<dmrginp.last_site()-1;++i)
      sites.push_back(i);
    StateInfo bigState;//The sys , dot and dot (the last site) combined together.
    TensorProduct(sysStateInfo, lastsiteState , bigState, PARTICLE_SPIN_NUMBER_CONSTRAINT);

    wave.SaveWavefunctionInfo(bigState, sites, wavenum);
  }

  void SpinAdapted::ReadQuanta(std::vector<SpinQuantum>& quanta, std::vector<int>& quantaStates, std::vector<int>& newquantaStates, int dotsite){
    //Read quanta from external DMRG code
    //The quanta was stored in the form of Block SpinQuantum
    std::vector<double> data; //The quanta is in numpy ndarray form. [N, S, totalStates]
    char file [500];
    sprintf (file, "%s%d", "quanta", dotsite);
    std::ifstream ifs(file, std::ios::binary| std::ios::ate);
    long size = ifs.tellg();
    assert(sizeof(double)==8);
    assert(size%32==0);
    data.resize(size/8);
    ifs.seekg (0, std::ios::beg);
    ifs.read((char *)(data.data()), size);
    ifs.close();
    quanta.resize(data.size()/4);
    quantaStates.resize(data.size()/4);
    newquantaStates.resize(data.size()/4);
    for(int i=0; i < quanta.size(); ++i)
    {
      quanta[i] = SpinQuantum((int)(data[4*i]+0.5), SpinSpace((int)(2*data[4*i+1]+0.5)), IrrepSpace());
      quantaStates[i] = (int) (data[4*i+2] + 0.5);
      newquantaStates[i] = (int) (data[4*i+3] + 0.5);

    }
  }

  void SpinAdapted::ReadRotationMatrix(std::vector<Matrix>& RotationM, std::vector<int>& quantaStates, std::vector<int>& newquantaStates, int dotsite){
    char file [500];
    RotationM.resize(quantaStates.size());
    std::vector<double> data;
    data.resize(10000);

    sprintf (file, "%s%d", "rotL", dotsite);
    std::ifstream ifs(file, std::ios::binary);
    assert(sizeof(double)==8);
    for(int i=0; i< RotationM.size(); ++i)
    {
      if((newquantaStates[i] == 0) || (quantaStates[i] == 0)) continue;
      RotationM[i].ReSize(quantaStates[i], newquantaStates[i]);
      ifs.read((char *)(RotationM[i].data()), sizeof(double)*quantaStates[i]*newquantaStates[i]);
    }

    ifs.close();
  }

