#include "MatrixBLAS.h"
#include "Stackspinblock.h"
#include "initblocks.h"
#include "input.h"
#include "timer.h"
#include <ctime>
#include "rotationmat.h"
#include "Stackwavefunction.h"
#include "global.h"
#include "sweep.h"
#include <unordered_map>
#include <unordered_set>
#include <boost/serialization/map.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <random>
#include "IntegralMatrix.h"
#include <chrono>
#include "stochasticpt.h"
#include "sampling.h"
#include "heatbath.h"
#include "simplemps.h"
using namespace SpinAdapted;
using namespace SpinAdapted::Sweep;
int bitstring::n_orb;
void ReadInput(char* conf);

void getComplementarySites(std::vector<int> &sites, std::vector<int> &complementarySites) ;



void compressedMPS(const unsigned long num_sample)
{
  //std::shared_ptr<NonAbelianmps> zeromps, sampledmps;
  std::shared_ptr<simplemps> zeromps, sampledmps;
  if(dmrginp.spinAdapted())
  {
    std::shared_ptr<NonAbelianmps> mps1 = std::make_shared<NonAbelianmps>();
    std::shared_ptr<NonAbelianmps> mps2 = std::make_shared<NonAbelianmps>();
    //Use their own shared_ptr for serialization in boost broadcast;
    mps1->build(0);
    mps2->build(1000);
    zeromps = mps1;
    sampledmps = mps2;
  }
  else{
    std::shared_ptr<Abelianmps> mps1 = std::make_shared<Abelianmps>();
    std::shared_ptr<Abelianmps> mps2 = std::make_shared<Abelianmps>();
    //Use their own shared_ptr for serialization in boost broadcast;
    mps1->build(0);
    mps2->build(1000);
    zeromps = mps1;
    sampledmps = mps2;
  }
  pout <<"QV|0> norm"<<sampledmps->get_norm()<<endl;
  /*
  {
  heatbath baseheatbath;
  baseheatbath.precompute();
    //sampledmps.build(1000);
    double norm = 0.0;
    double largest = 0.0;
    double h00 = 0.0;
    double h11 = 0.0;
    double h10 = 0.0;
    for(long n=0;n<pow(2,dmrginp.last_site()*2);n++)
    {
      bitstring ci(n);
    int particle_num = 0;
    for(int i=0;i<ci.size();i++)
      particle_num+= ci[i];
    if (particle_num!=dmrginp.total_particle_number()) continue;

    double coeff = zeromps->getcoeff(ci);
    norm +=coeff*coeff;
    largest = max(coeff*coeff,largest);
    h00 += coeff*coeff/local_energy(ci, 0);

      std::unordered_map<bitstring, double> sd_hashtable1;
      baseheatbath.allexcite(ci, 1.0,sd_hashtable1, 1e-13);
      sd_hashtable1[ci] += local_energy(ci, 1);
      double overlap=0.0;
      for(auto iter1=sd_hashtable1.begin();iter1!=sd_hashtable1.end();iter1++)
      {
        overlap += zeromps->getcoeff(iter1->first)*iter1->second;
      }
      h11 += overlap*overlap/local_energy(ci, 0);
      h10 += overlap*zeromps->getcoeff(ci)/local_energy(ci, 0);
      sd_hashtable1.clear();

    }
    cout <<"Norm: " <<norm<<endl;
    cout <<"Largest: " <<largest<<endl;
    cout <<"h00"<<h00<<endl;
    cout <<"h11"<<h11<<endl;
    cout <<"h10"<<h10<<endl;
    abort();
    cout <<"Begin sample"<<endl;
    for(int i=0;i<10;i++)
    {

      bitstring determinant;
      double coeff = zeromps->sampling(determinant);
      cout <<"ul"<<determinant<<endl;
      coeff = zeromps->getcoeff(determinant);
      cout <<"coeff*coeff"<< coeff*coeff<<endl;

    }
    abort();

  }
  */

  //get_slater();
  std::vector<int> ci;
  //const unsigned long num_sample = 20000;
  heatbath baseheatbath;
  //baseheatbath.set_tol() = dmrginp.stochasticpt_tol();
  baseheatbath.precompute();
  double H11=0.0;
  double H11_2=0.0;
  double H11_new=0.0;
  double H11_2_new=0.0;
  double H00=0.0;
  double H00_2=0.0;
  double H01=0.0;
  double H01_2=0.0;
  double expectation1=0.0;
  double test0=0.0;
  double test0_2=0.0;
  double test1=0.0;
  double test1_2=0.0;
  double test2=0.0;
  double test2_2=0.0;
  double PT_onesample = 0.0;
  double PT_onesample_2 = 0.0;
  double PT_twosample = 0.0;
  double PT_twosample_2 = 0.0;
  int n_H_samples = dmrginp.stochasticpt_Hsamples();
  std::unordered_map<bitstring, double> sd_hashtable0;
  std::unordered_map<bitstring, double> sd_hashtable0_p;
  std::unordered_map<bitstring, double> sd_hashtable1;
  std::unordered_map<bitstring, double> sd_hashtable1sum;
  std::unordered_map<bitstring, double> det_energy;
  std::unordered_map<bitstring, double> det_coeff;

  double minele = 1000.0;
  double maxele = 0.0;
  double timer1=0.0;
  double timer2=0.0;
  double timer3=0.0;
  double tol = dmrginp.stochasticpt_tol();
  std::clock_t startTime = std::clock();
  std::clock_t begintime;
  //sd_hashtable.reserve(num_sample);
  {
    begintime = std::clock();
    double pt2=0.0;
    double pt2_one=0.0;
    double pt2_two=0.0;
    double norm= 0.0;
    double H01_temp=0.0;
    pout <<"Number of samples: "<< mpigetsize() <<" X "<< num_sample<<endl;
    cout <<"Process "<<mpigetrank()<<": began."<<endl;
    for (int i=0;i<num_sample;i++)
    {
      bitstring determinant;
      double coeff = zeromps->sampling(determinant);
      //cout <<"coeff"<<coeff<<endl;
      //cout <<"coeff^2"<<coeff*coeff<<endl;
      //double coeff = baseState.sampling(determinant);
      double local_energy0;
      if(det_energy.find(determinant)==det_energy.end())
      {
        local_energy0= baseheatbath.local_energy(determinant, 0);
        det_energy[determinant] = local_energy0;
      }
      else{
        local_energy0= det_energy[determinant];

      }
      H00 += 1/(local_energy0*num_sample);
      H00_2 += 1/(local_energy0*local_energy0*num_sample);
      //norm += 1/(coeff*coeff);
    }
    timer1 += (std::clock() - begintime)/ (double) CLOCKS_PER_SEC;
    begintime = std::clock();
    long num_coeff=0;
    for (int i=0;i<num_sample;i++)
    {

      bitstring determinant;
      //double prob;
      double coeff = sampledmps->sampling(determinant);

      //cout <<"coeff"<<coeff<<endl;
      //cout <<"coeff^2"<<coeff*coeff/sampledmps->get_norm()<<endl;
      double e = baseheatbath.local_energy(determinant, 0);

      double overlap = baseheatbath.Expectation(determinant, zeromps.get(), fabs(tol));
      //baseheatbath.allexcite(determinant, 1.0,sd_hashtable1, fabs(tol*100));
      ////baseheatbath.allexcite(iter->first, 1.0,sd_hashtable1, fabs(tol/coeff));
      //sd_hashtable1[determinant] += local_energy(determinant, 1);
      //num_coeff += sd_hashtable1.size();

      //double overlap=0.0;
      //for(auto iter1=sd_hashtable1.begin();iter1!=sd_hashtable1.end();iter1++)
      //{
      //  overlap +=  zeromps->getcoeff(iter1->first)*iter1->second;
      //  //testmps.getcoeff(iter1->first);
      //}
      double temp =sampledmps->get_norm()/(e);
      test1 += temp/num_sample;
      test1_2 += temp*temp/num_sample;
      //TODO 
      //coeff can be very small number, however, it is close to overlap.
      //Allways use overlap/coeff. It could reduce variance a lot.
      #pragma optimize( "", off )
      double ratio = overlap/coeff;
      temp *= ratio*ratio;
      #pragma optimize( "", on ) 
      H11 += temp/num_sample;
      H11_2 += temp*temp/num_sample;

      temp = sampledmps->get_norm()*overlap*zeromps->getcoeff(determinant)/(e*coeff*coeff);
      H01 += temp/num_sample;
      H01_2 += temp*temp/num_sample;


      sd_hashtable1.clear();

      if(i%(max(num_sample/10,(const unsigned long)1))==0)
      {
        //pout <<"Finished " <<i <<" samples"<<endl;
        cout <<"Process "<<mpigetrank()<<": Finished " <<i <<" samples"<<endl;
      }

      //sd_hashtable0[determinant] += 1.0;
    }
    timer2 += (std::clock() - begintime)/ (double) CLOCKS_PER_SEC;
    begintime = std::clock();

    pout <<"Num of coeff to compute: "<< num_coeff<<endl;
    timer3 += (std::clock() - begintime)/ (double) CLOCKS_PER_SEC;
  }

  cout <<"Process "<<mpigetrank()<<": Finished."<<endl;
  cout <<"Process "<<mpigetrank()<<": Approx sampling |c| time: " <<(std::clock() - startTime)/ (double) CLOCKS_PER_SEC<<endl;


  {
  boost::mpi::communicator world;
  PT_onesample = all_reduce(world, PT_onesample, std::plus<double>())/world.size();
  PT_onesample_2 = all_reduce(world, PT_onesample_2, std::plus<double>())/world.size();
  PT_twosample = all_reduce(world, PT_twosample, std::plus<double>())/world.size();
  PT_twosample_2 = all_reduce(world, PT_twosample_2, std::plus<double>())/world.size();
  H11 = all_reduce(world, H11, std::plus<double>())/world.size();
  H11_2 = all_reduce(world, H11_2, std::plus<double>())/world.size();
  H00 = all_reduce(world, H00, std::plus<double>())/world.size();
  H00_2 = all_reduce(world, H00_2, std::plus<double>())/world.size();
  H01 = all_reduce(world, H01, std::plus<double>())/world.size();
  H01_2 = all_reduce(world, H01_2, std::plus<double>())/world.size();
  test0 = all_reduce(world, test0, std::plus<double>())/world.size();
  test0_2 = all_reduce(world, test0_2, std::plus<double>())/world.size();
  test1 = all_reduce(world, test1, std::plus<double>())/world.size();
  test1_2 = all_reduce(world, test1_2, std::plus<double>())/world.size();
  test2 = all_reduce(world, test2, std::plus<double>())/world.size();
  test2_2 = all_reduce(world, test2_2, std::plus<double>())/world.size();
  pout <<"TEST0: " <<test0<<" + "<<sqrt((test0_2-test0*test0)/(num_sample*world.size()))<<endl;
  pout <<"TEST1: " <<test1<<" + "<<sqrt((test1_2-test1*test1)/(num_sample*world.size()))<<endl;
  pout <<"TEST2: " <<test2<<" + "<<sqrt((test2_2-test2*test2)/(num_sample*world.size()))<<endl;
  //pout <<"PT one sample" << PT_onesample<<" + "<<sqrt((PT_onesample_2-PT_onesample*PT_onesample)/(num_sample*world.size()))<<endl;
  //pout <<"PT two sample" << PT_twosample<<" + "<<sqrt((PT_twosample_2-PT_twosample*PT_twosample)/(num_sample*world.size()))<<endl;
  pout <<"<0|VQ(H_0-E_0)^(-1)QV|0>" << H11<<" + "<<sqrt((H11_2-H11*H11)/(num_sample*world.size()))<<endl;
  pout <<"<0|(H_0-E_0)^(-1)|0>" << H00<<" + "<<sqrt((H00_2-H00*H00)/(num_sample*world.size()))<<endl;
  pout <<"<0|(H_0-E_0)^(-1)QV|0>" << H01<<" + "<<sqrt((H01_2-H01*H01)/(num_sample*world.size()))<<endl;
  pout <<"PT energy: " << -H11+H01*H01/H00 << " + " <<sqrt((H11_2-H11*H11)/(num_sample*world.size()))<<endl;
  }

  if(mpigetrank() ==0)
  {
    cout <<"Timer for head node"<<endl;
    zeromps->print_timer();
    sampledmps->print_timer();
    cout <<"Single Excitation: " << baseheatbath.singleT<<endl;
    cout <<"Double Excitation: " << baseheatbath.doubleT<<endl;
    cout <<"Excitation: " << baseheatbath.excitation_T<<endl;
    cout <<"Energy: " << baseheatbath.energy_T<<endl;
    cout <<"factor: " << baseheatbath.factor_T<<endl;
    cout <<"Phase1: " << timer1<<endl;
    cout <<"Phase2: " << timer2<<endl;
    cout <<"Phase3: " << timer3<<endl;
  }

  return;


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
  //TODO
  //region = boost::interprocess::mapped_region{segment, boost::interprocess::read_only};

  //int m_norbs = dmrginp.spinAdapted()?dmrginp.last_site()*2:dmrginp.last_site();
  //long oneIntegralMem = (m_norbs/2*(m_norbs/2+1))/2; 
  //long twoedim =v_2[0].matDim;
  //long twoIntegralMem = twoedim*(twoedim+1)/2;
  //for(int i=0;i<v_1.size();i++)
  //{
  //  v_1[i].set_data() = static_cast<double*>(region.get_address()) + (oneIntegralMem+twoIntegralMem)*i;
  //  v_2[i].set_data() = static_cast<double*>(region.get_address()) + oneIntegralMem + (oneIntegralMem+twoIntegralMem)*i;

  //}
  //dmrginp.matmultFlops.resize(1, 0.);
  dmrginp.initCumulTimer();
  bitstring::n_orb= dmrginp.spinAdapted()?dmrginp.last_site()*2:dmrginp.last_site();



  const int numthrd=1;
  dmrginp.matmultFlops.resize(numthrds, 0.);

  int size = 1, rank =0;

  MAX_THRD = dmrginp.thrds_per_node()[mpigetrank()];
  int mkl_thrd = dmrginp.mkl_thrds();

  #ifdef _OPENMP
    omp_set_num_threads(MAX_THRD);
  #ifdef _HAS_INTEL_MKL 
    mkl_set_num_threads(mkl_thrd);
    mkl_set_dynamic(0);
  #endif
  omp_set_nested(1);
  #endif
  cout.precision (12);
  cout << std::fixed;
  double* stackmemory = new double[dmrginp.getMemory()];
  Stackmem.resize(numthrds);
  Stackmem[0].data = stackmemory;
  Stackmem[0].size = dmrginp.getMemory();
  //dmrginp.initCumulTimer();



  long num_sample = dmrginp.stochasticpt_nsamples();
  //check_heatbath(num_sample);
  //check_sampling_approx(num_sample);
  //exactpt();
  //printoutall_combined(num_sample);
  //printoutall(num_sample);
  //printoutall_twobatches(num_sample);
  //printoutall_twoH(num_sample);
  //splitspace(num_sample);
  //hci(num_sample);
  //check_sampling_approx_combined(num_sample);
  //check_sampling_approx(num_sample);
  compressedMPS(num_sample);
  //approxsampling_twoH(num_sample);

  //check_overlap(num_sample);
  //double cputime = globaltimer.totalcputime();
  //double walltime = globaltimer.totalwalltime();
  //pout << setprecision(3) <<"\n\n\t\t\t BLOCK CPU  Time (seconds): " << cputime << endl;
  //pout << setprecision(3) <<"\t\t\t BLOCK Wall Time (seconds): " << walltime << endl;
  delete [] stackmemory;
  return 0;
}

