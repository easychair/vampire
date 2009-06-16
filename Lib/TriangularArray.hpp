/**
 * @file TriangularArray.hpp
 * Defines class TriangularArray.
 */


#ifndef __TriangularArray__
#define __TriangularArray__

#include "Allocator.hpp"

namespace Lib {

template<typename T>
class TriangularArray
{
private:
  //private and undefined operator= and copy constructor to avoid implicitly generated ones
  TriangularArray(const TriangularArray&);
  TriangularArray& operator=(const TriangularArray&);
public:
  TriangularArray(size_t side)
  : _2sideMinus1(2*side-1), _data(0)
  {
    ASS_G(side,0);
    _capacity=dataSize();
    void * mem=ALLOC_KNOWN(_capacity*sizeof(T), "Lib::TriangularArray");
    _data=new(mem) T[_capacity];
  }
  ~TriangularArray()
  {
    T* p=_data+_capacity;
    while(p!=_data) {
      (--p)->~T();
    }
    DEALLOC_KNOWN(_data,_capacity*sizeof(T), "Lib::TriangularArray");
  }

  void setSide(size_t side)
  {
    ASS_G(side,0);
    _2sideMinus1=2*side-1;
    if(dataSize()>_capacity) {
      size_t newCapacity=max(_capacity*2,dataSize());
      void * mem=ALLOC_KNOWN(newCapacity*sizeof(T), "Lib::TriangularArray");

      T* p=_data+_capacity;
      while(p!=_data) {
        (--p)->~T();
      }
      DEALLOC_KNOWN(_data,_capacity*sizeof(T), "Lib::TriangularArray");

      _capacity=newCapacity;
      _data=new(mem) T[_capacity];
    }
  }

  inline
  size_t dataSize() const
  {
    return ((_2sideMinus1+1)*(_2sideMinus1+3)) / 8;
  }

  T& get(size_t x, size_t y)
  {
    ASS_GE(x,y);
    ASS_L(x,(_2sideMinus1+1)/2);
    ASS_L(x + (y*(_2sideMinus1-y))/2, _capacity);
    return _data[x + (y*(_2sideMinus1-y))/2];
  }

  void set(size_t x, size_t y, T val)
  {
    ASS_GE(x,y);
    ASS_L(x,(_2sideMinus1+1)/2);
    ASS_L(x + (y*(_2sideMinus1-y))/2, _capacity);
    _data[x + (y*(_2sideMinus1-y))/2]=val;
  }



private:
  /**
   * It's more convenient to have the side times two minus 1,
   * than simply side, as we save two operations in the
   * access operation.
   */
  size_t _2sideMinus1;
  size_t _capacity;
  T* _data;
};

};

#endif /* __TriangularArray__ */
