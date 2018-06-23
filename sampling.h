#ifndef SAMPLING_H
#define SAMPLING_H

#include <vector>
#include <bitset>
#include "global.h"

template <typename T>
struct counter
{
    static int objects_created;
    static int objects_alive;

    counter()
    {
        ++objects_created;
        ++objects_alive;
    }
    
    counter(const counter&)
    {
        ++objects_created;
        ++objects_alive;
    }
  protected:
    ~counter() // objects should never be removed through pointers of this type
    {
        --objects_alive;
    }
};
template <typename T> int counter<T>::objects_created( 0 );
template <typename T> int counter<T>::objects_alive( 0 );

class longbitarray 
{
  public:
  static const int size=6; //num of 64-bit integers used for ci string.
  static const uint64_t one;
  std::array<uint64_t, size> data;
  longbitarray(){unset();}
  longbitarray(const std::vector<int>& sd){
    unset();
    for(int i=0;i<sd.size();i++)
    {
      if(sd[i])
        set(i);
      else  
        unset(i);
    }
  }
  bool operator==(const longbitarray& b) const
  {
    for(int i=0;i<size;i++)
      if (data[i]!=b.data[i])
        return false;
    return true;
  }
  inline int getocc(int i) const{ return (data[i/64] >> (i%64)) & (uint64_t)1;}
  //inline void flip(int i){ data[i/64] = data[i/64] ^ ((uint64_t)1 <<(i%64)) ;}
  //inline void set(int i){ data[i/64] = data[i/64] | ((uint64_t)1 <<(i%64)) ;}
  //inline void unset(int i){ set(i); flip(i);}
  inline void flip(int i){ data[i/64] ^= ((uint64_t)1 <<(i%64));}
  inline void set(int i){ data[i/64] |= ((uint64_t)1 <<(i%64));}
  inline void unset(int i){ data[i/64] &= ~((uint64_t)1 <<(i%64));}
  inline void unset() {for(int i=0;i<size;i++) data[i] = 0;}
  inline std::vector<int> to_vector(int n) const {
    std::vector<int> sd(n);
    for(int i=0;i<n;i++)
    {
      sd[i] = getocc(i);
    }
    return sd;
  }
  inline void from_vector(const std::vector<int>& sd){
    unset();
    for(int i=0;i<sd.size();i++)
    {
      if(sd[i])
        set(i);
      else  
        unset(i);
    }
  }
};

namespace std
{
  template <>
  struct hash<longbitarray>
  {
    size_t operator()(const longbitarray& k) const
      {
        return hash<uint64_t>()(k.data[0]);
        //return hash<uint64_t>()(k.data.back());
      }
  };
}

const int bitlength=64*6;
//typedef std::bitset<bitlength> bitstring;
class bitstring : public std::bitset<bitlength>
{
  public:
    //friend class std::hash<bitset<bitlength> >;
    static int n_orb;
    size_t length() const {return n_orb;}
    size_t size() const {return n_orb;}
    bitstring(long n) : std::bitset<bitlength>(n){}
    bitstring() : std::bitset<bitlength>(0){}
    const int Sz() const{
      int sz=0;
      for(int i=0;i<n_orb/2;i++)
        sz += this->operator[](2*i) -this->operator[](2*i+1);
      return sz;
    }
    int excitation(int q, int p)
    {
      //Return permutation factor for excitation operators.
      assert(q<n_orb && p <n_orb);
      assert(this->operator[](p));
      assert(!(this->operator[](q)));
      this->reset(p);
      this->set(q);
      int nelec = 0;
      auto newbits = *this;
      newbits >>= min(p,q)+1;
      newbits <<= min(p,q) +(bitlength-max(p,q)) +1;
      return newbits.count()%2? -1:1;
      //for(int i=min(p,q)+1; i<max(p,q);i++)
      //  nelec += this->operator[](i);
      //return nelec%2? -1:1;
    }
    int excitation(int r, int s, int q, int p)
    {
      //Return permutation factor for excitation operators.
      assert(q<n_orb && p <n_orb && r<n_orb && s<n_orb);
      assert(this->operator[](p));
      assert(this->operator[](q));
      assert(!(this->operator[](r)));
      assert(!(this->operator[](s)));
      int nelec = 0;

      auto newbits = *this;
      newbits >>= min(p,q)+1;
      newbits <<= min(p,q) +(bitlength-max(p,q)) +1;
      nelec += newbits.count();

      //for(int i=min(p,q)+1; i<max(p,q);i++)
      //  nelec += this->operator[](i);
      this->reset(p);
      this->reset(q);

      newbits = *this;
      newbits >>= min(r,s)+1;
      newbits <<= min(r,s) +(bitlength-max(r,s)) +1;
      nelec += newbits.count();

      //for(int i=min(r,s)+1; i<max(r,s);i++)
      //  nelec += this->operator[](i);
      this->set(r);
      this->set(s);
      int factor = nelec%2? -1:1;
      if(r>s) factor *=-1;
      if(p>q) factor *=-1;
      return factor;
    }
};

namespace std
{
  template <>
  struct hash<bitstring>
  {
    size_t operator()(const bitstring& k) const
      {
        return hash<bitset<bitlength> >()(k);
        //return hash<uint64_t>()(k.data.back());
      }
  };
}


#endif
