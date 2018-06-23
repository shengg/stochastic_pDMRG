#include <vector>
#include "global.h"
#include "input.h"
#include "heatbath.h"
#include <ctime>
#include <random>
#include "IntegralMatrix.h"
#include <chrono>

using namespace SpinAdapted;
using namespace std;

void heatbath::precompute(){

  int num_spatial_orbs = dmrginp.spinAdapted()?dmrginp.last_site():dmrginp.last_site()/2;

  Direct.ReSize(num_spatial_orbs, num_spatial_orbs);
  Exchange.ReSize(num_spatial_orbs, num_spatial_orbs);
  //Direct =Matrix(num_spatial_orbs, num_spatial_orbs);
  //Exchange =Matrix(num_spatial_orbs, num_spatial_orbs);
  for(int p=0;p<num_spatial_orbs;p++)
    for(int q=0;q<num_spatial_orbs;q++)
    {
      Direct(p+1,q+1) = v_2[0](2*p,2*q,2*p,2*q);
      Exchange(p+1,q+1) = v_2[0](2*p,2*q,2*q,2*p);
    }

  /*
  {

  doubleH.resize(num_spatial_orbs*2);
  for(int p=0;p<num_spatial_orbs*2;p++)
  {
    doubleH[p].resize(num_spatial_orbs*2);
    for(int q=0;q<p;q++)
    {
      //doubleH.resize(num_spatial_orbs*2*num_spatial_orbs*2);
      for(int r=0;r<num_spatial_orbs*2;r++)
      for(int s=0;s<r;s++)
      {
        if(s!=p && s!=q && r!=p && r!=q )
        //if(s!=r &&s!=p && s!=q && r!=p && r!=q )
        {
          double h = v_2[integralIndex](r,s,p,q) -v_2[integralIndex](r,s,q,p)+v_2[integralIndex](s,r,q,p) -v_2[integralIndex](s,r,p,q) ;
        //doubleH[p][q].push_back(std::make_tuple(v_2[integralIndex](r,s,p,q), r, s));
          if(fabs(h)>tol)
            doubleH[p][q].push_back(std::make_tuple(h, r, s));

        }

      }
      std::sort(doubleH[p][q].begin(),doubleH[p][q].end(),[](std::tuple<double, int, int> a, std::tuple<double, int, int> b){return fabs(std::get<0>(a)) > abs(std::get<0>(b));});
    }
  }
  }
  */

  /*
  if(mpigetrank()!=0)
  {
    for(int i=0;i <v_2.size();i++)
      v_2[i].ReSize(0);
  }
  else{

      segment.truncate((v_1[0].GetRepresentation().Storage()+v_2[0].GetRepresentation().Storage())*m_num_Integrals*sizeof(double)); 
      region = boost::interprocess::mapped_region{segment, boost::interprocess::read_write};
      memset(region.get_address(), 0., (oneIntegralMem+twoIntegralMem)*m_num_Integrals*sizeof(double));
    v1.set_data() = static_cast<double*>(region.get_address()) + (oneIntegralMem+twoIntegralMem)*integralIndex;
    v2.set_data() = static_cast<double*>(region.get_address()) + oneIntegralMem + (oneIntegralMem+twoIntegralMem)*integralIndex;
  }

  #ifndef SERIAL
  boost::mpi::communicator world;
  broadcast(world, doubleH, 0);
  #endif
*/

}


    void heatbath::allexcite(const std::vector<int>& sd_in, double coeff,  std::unordered_map<longbitarray, double>& sd_table, double tol){
      std::clock_t startTime = std::clock();
      double sample_p = 1.0;
      std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
      longbitarray sdbits(sd_in);
      std::vector<int> occupy;
      std::vector<int> unoccupy;
      for(int i=0;i<sd_in.size();i++)
      {
        if(sd_in[i]==1) occupy.push_back(i);
        else unoccupy.push_back(i);
      }
      for(int p: occupy)
        for(int q: occupy)
          if(p!=q)
          {
            for(int i=0;i<doubleH[p][q].size();i++)
            {
              if(fabs(std::get<0>(doubleH[p][q][i])) < abs(tol/coeff))
                break;
              if(sdbits.getocc(std::get<1>(doubleH[p][q][i]))) continue;
              if(sdbits.getocc(std::get<2>(doubleH[p][q][i]))) continue;
              longbitarray newsdbits = sdbits;
              newsdbits.unset(p);
              newsdbits.unset(q);
              newsdbits.set(std::get<1>(doubleH[p][q][i]));
              newsdbits.set(std::get<2>(doubleH[p][q][i]));
              int factor = permutefactor(sd_in, std::get<1>(doubleH[p][q][i]), std::get<2>(doubleH[p][q][i]), q, p);
              sd_table[newsdbits] += 0.5*coeff*std::get<0>(doubleH[p][q][i])*factor;
            }
          }
      for(int p: occupy)
      {
        for(int r: unoccupy)
        {

          double h = v_1[integralIndex](r,p);
          for(int q: occupy)
          {
            h += q==p? 0.0: v_2[integralIndex](r,q, p, q);
            h -= p==q? 0.0: v_2[integralIndex](r,q, q, p);
          }
          if(fabs(h)> abs(tol/coeff))
          {
            longbitarray newsdbits = sdbits;
            newsdbits.unset(p);
            newsdbits.set(r);
            int factor = permutefactor(sd_in, r, p);
            sd_table[newsdbits] += coeff*h*factor;
          }
        }
      }
      excitation_T+= (std::clock() -startTime)/(double) CLOCKS_PER_SEC;

    }

    void heatbath::allexcite(const bitstring& sd_in, double coeff,  std::unordered_map<bitstring, double>& sd_table, double tol){
      std::clock_t startTime = std::clock();
      double sample_p = 1.0;
      //std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
      std::vector<int> occupy;
      std::vector<int> unoccupy;
      for(int i=0;i<sd_in.size();i++)
      {
        if(sd_in[i]==1) occupy.push_back(i);
        else unoccupy.push_back(i);
      }
      for(int p=0;p<occupy.size();p++)
        for(int q=0;q<p;q++)
          for(int r=0;r<unoccupy.size();r++)
            for(int s=0;s<r;s++)
          {
            int i=occupy[p], j=occupy[q], k=unoccupy[r], l=unoccupy[s];

            double h = v_2[integralIndex](k,l,i,j) - v_2[integralIndex](l,k,i,j)+v_2[integralIndex](l,k,j,i) - v_2[integralIndex](k,l,j,i);
              if(fabs(h) < fabs(tol/coeff))
                continue;
              bitstring newsdbits = sd_in;
              int factor = newsdbits.excitation(k,l,j,i);
              if(sd_table.find(newsdbits)!=sd_table.end())
              {
                cout <<"determinant Used"<<endl;
              }
              sd_table[newsdbits] += 0.5*coeff*h*factor;
          }
      for(int p: occupy)
      {
        for(int r: unoccupy)
        {

          double h = v_1[integralIndex](r,p);
          for(int q: occupy)
          {
            h += q==p? 0.0: v_2[integralIndex](r,q, p, q);
            h -= p==q? 0.0: v_2[integralIndex](r,q, q, p);
          }
          if(fabs(h)> abs(tol/coeff))
          {
            bitstring newsdbits = sd_in;
            int factor = newsdbits.excitation(r, p);
            sd_table[newsdbits] += coeff*h*factor;
          }
        }
      }

      excitation_T+= (std::clock() -startTime)/(double) CLOCKS_PER_SEC;
    }

    double heatbath::Expectation(const bitstring& sd_in, simplemps* zeromps, double tol){
      std::clock_t startTime = std::clock();
      double o = local_energy(sd_in,1)*zeromps->getcoeff(sd_in);
      double sample_p = 1.0;
      //std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
      std::vector<int> occupy;
      std::vector<int> unoccupy;
      for(int i=0;i<sd_in.size();i++)
      {
        if(sd_in[i]==1) occupy.push_back(i);
        else unoccupy.push_back(i);
      }
      for(int p=0;p<occupy.size();p++)
        for(int q=0;q<p;q++)
          for(int r=0;r<unoccupy.size();r++)
            for(int s=0;s<r;s++)
          {
            int i=occupy[p], j=occupy[q], k=unoccupy[r], l=unoccupy[s];

            double h = v_2[integralIndex](k,l,i,j) - v_2[integralIndex](l,k,i,j)+v_2[integralIndex](l,k,j,i) - v_2[integralIndex](k,l,j,i);
              if(fabs(h) < fabs(tol))
                continue;
              bitstring newsdbits = sd_in;
              int factor = newsdbits.excitation(k,l,j,i);

              o += 0.5*h*factor*zeromps->getcoeff(newsdbits);
          }
      for(int p: occupy)
      {
        for(int r: unoccupy)
        {

          double h = v_1[integralIndex](r,p);
          for(int q: occupy)
          {
            h += q==p? 0.0: v_2[integralIndex](r,q, p, q);
            h -= p==q? 0.0: v_2[integralIndex](r,q, q, p);
          }
          if(fabs(h)> abs(tol))
          {
            bitstring newsdbits = sd_in;
            int factor = newsdbits.excitation(r, p);
            o += h*factor*zeromps->getcoeff(newsdbits);
          }
        }
      }

      excitation_T+= (std::clock() -startTime)/(double) CLOCKS_PER_SEC;
      return o;
    }

double heatbath::EnergyAfterExcitation(vector<int>& occupy, int i, int j, int a, int b, double Energyd) {
  /*!
     Calculates the new energy of a determinant after double excitation.

     .. note:: Assumes that i!=j and a!=b

     :Arguments:

       vector<int>& closed:
           Occupied orbitals in a vector.
       int i:
           Orbital index for destruction operator.
       int j:
           Orbital index for destruction operator.
       int a:
           Orbital index for creation operator.
       int b:
           Orbital index for creation operator.
       double Energyd:
           Old determinant energy.

     :Returns:

        double E:
            Energy after excitation.
   */

  std::clock_t startTime = std::clock();
  assert(i%2==a%2);
  assert(j%2==b%2);
  double E = Energyd;
  E += v_1[integralIndex](a, a)+v_1[integralIndex](b, b) -v_1[integralIndex](i, i)-v_1[integralIndex](j, j)   ;

/*
  for (int I: occupy) {
    if (I == i|| I ==j) continue;
    E +=  v_2[integralIndex](I,a,I,a)-v_2[integralIndex](I,i,I,i)+v_2[integralIndex](I,b,I,b)-v_2[integralIndex](I,j,I,j);
    E +=  -v_2[integralIndex](I,a, a, I)+v_2[integralIndex](I,i,i, I)-v_2[integralIndex](I,b, b, I)+v_2[integralIndex](I,j,j, I);
}
  E += v_2[integralIndex](a, b, a, b) - v_2[integralIndex](a, b, b, a) -v_2[integralIndex](i, j, i, j) +v_2[integralIndex](i, j, j, i);
*/
  for(int I: occupy)
  {
    if (I==i || I==j) continue;
    E += Direct(a/2+1, I/2+1) - Direct(i/2+1, I/2+1);
    if (i%2 == I%2)
      E += -Exchange(a/2+1, I/2+1) + Exchange(i/2+1, I/2+1);
  }
  for(int I: occupy)
  {
    if (I==i || I==j) continue;
    E += Direct(b/2+1, I/2+1) - Direct(j/2+1, I/2+1);
    if (j%2 == I%2)
      E += -Exchange(b/2+1, I/2+1) + Exchange(j/2+1, I/2+1);
  }
  E += Direct(a/2+1, b/2+1) - Direct(i/2+1, j/2+1);
  if (i%2==j%2)
  E += -Exchange(a/2+1, b/2+1) + Exchange(i/2+1, j/2+1);

  energy_T += (std::clock() -startTime)/ (double) CLOCKS_PER_SEC;
  return E;
}

double heatbath::EnergyAfterExcitation(vector<int>& occupy, int i, int a, double Energyd) {
  /*!
     Calculates the new energy of a determinant after single excitation.

     .. note:: Assumes that the spin of i and a orbitals is the same

     :Arguments:

       vector<int>& closed:
           Occupied orbitals in a vector.
       int i:
           Orbital index for destruction operator.
       int a:
           Orbital index for creation operator.
       double Energyd:
           Old determinant energy.

     :Returns:

        double E:
            Energy after excitation.
   */

  std::clock_t startTime = std::clock();
  double E = Energyd;
  E += v_1[integralIndex](a, a)-v_1[integralIndex](i, i);

  for (int I :occupy) {
    if (I == i) continue;
    //E +=  v_2[integralIndex](I,a,I,a)-v_2[integralIndex](I,i,I,i);
    E +=  Direct(I/2+1, a/2+1) - Direct(I/2+1, i/2+1);
    if ( (I%2) == (i%2) )
    //E +=  -v_2[integralIndex](I,a, a, I)+v_2[integralIndex](I,i,i, I);
    E +=  -Exchange(I/2+1, a/2+1) + Exchange(I/2+1, i/2+1);
  }
  energy_T += (std::clock() -startTime)/ (double) CLOCKS_PER_SEC;
  return E;
}

double heatbath::local_energy(const bitstring& ci, int integralIndex)
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
        //energy += Direct(occupy[i]/2+1,occupy[j]/2+1);
        //if(i%2==j%2)
        //  energy -= Exchange(occupy[j]/2+1,occupy[i]/2+1);
      }
      //HdiagonalT += (std::clock() -startTime)/ (double) CLOCKS_PER_SEC;
      energy_T += (std::clock() -startTime)/(double) CLOCKS_PER_SEC;
      return energy;
    
    }


