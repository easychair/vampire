/*
 * This file is part of the source code of the software program
 * Vampire. It is protected by applicable
 * copyright laws.
 *
 * This source code is distributed under the licence found here
 * https://vprover.github.io/license.html
 * and in the source directory
 */
/**
 * @file RuntimeStatistics.cpp
 * Implements class RuntimeStatistics.
 */

#include "Lib/Comparison.hpp"
#include "Lib/Int.hpp"
#include "Lib/Sort.hpp"
#include "Lib/Stack.hpp"

#include "RuntimeStatistics.hpp"

#if RUNTIME_STATS

namespace Debug
{

using namespace std;

void RSMultiCounter::print(std::ostream& out)
{
  out << "% " << name() << ":"<< endl;
  for(size_t i=0;i<_counters.size();i++) {
    if(_counters[i]) {
      out << "  " << i << ": " << _counters[i] <<endl;
    }
  }
}

RSMultiStatistic::~RSMultiStatistic()
{
  for(size_t i=0;i<_values.size();i++) {
    ValList::destroy(_values[i]);
  }
}

void RSMultiStatistic::print(std::ostream& out)
{
  out << "% " << name() << ":"<< endl;
  for(size_t i=0;i<_values.size();i++) {
    if(_values[i]) {

      ValList* vals=_values[i];

      ValList::Iterator vit(vals);
      size_t cnt=0;
      long long int sum=0;
      int min=vals->head();
      int max=vals->head();
      while(vit.hasNext()) {
	int val=vit.next();
	cnt++;
	sum+=val;
	if(val<min) {
	  min=val;
	}
	if(val>max) {
	  max=val;
	}
      }
      out << "  " << i << ": " <<
              "cnt: "+Int::toString(cnt)+
              ", avg: "+Int::toString(static_cast<float>(sum)/cnt)+
              ", min: "+Int::toString(min)+
              ", max: "+Int::toString(max) << endl;
    }
  }
}


RuntimeStatistics* RuntimeStatistics::instance()
{
  static RuntimeStatistics inst;

  return &inst;
}

RuntimeStatistics::RuntimeStatistics()
{
}

RuntimeStatistics::~RuntimeStatistics()
{
  while (_objs.isNonEmpty()) {
    delete _objs.pop();
  }
}

void RuntimeStatistics::print(std::ostream& out)
{
  out<<"% ----  Runtime statistics ----"<<endl;

  ObjSkipList::Iterator it(_objs);
  while(it.hasNext()) {
    it.next()->print(out);
  }

  out<<"% -----------------------------"<<endl;
}


}

#endif
