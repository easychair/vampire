/*
 * This file is part of the source code of the software program
 * Vampire. It is protected by applicable
 * copyright laws.
 *
 * This source code is distributed under the licence found here
 * https://vprover.github.io/license.html
 * and in the source directory
 */

inline DArray<int> predLevels() {
  DArray<int> out(env.signature->predicates());
  out.init(out.size(), 1);
  return out;
}

using namespace Kernel;

template<class SigTraits>
inline KboWeightMap<SigTraits> toWeightMap(unsigned introducedSymbolWeight, KboSpecialWeights<SigTraits> ws, const Map<unsigned, KboWeight>& xs, unsigned sz) 
{
  auto df = KboWeightMap<SigTraits>::dflt();
  df._specialWeights = ws;

  DArray<KboWeight> out(sz);
  for (unsigned i = 0; i < sz; i++) {
    auto w = xs.getPtr(i);
    out[i] = w == NULL ? df.symbolWeight(i) : *w;
  }
  return  {
    ._weights = out.clone(),
    ._introducedSymbolWeight = introducedSymbolWeight,
    ._specialWeights         = ws,
  };
}

inline void __weights(Map<unsigned, KboWeight>& ws) {
}

template<class A, class... As>
inline void __weights(Map<unsigned, KboWeight>& ws, std::pair<A, KboWeight> a, std::pair<As, KboWeight>... as) {
  ws.insert(std::get<0>(a).functor(), std::get<1>(a));
  __weights(ws, as...);
}

template<class... As>
inline Map<unsigned, KboWeight> weights(std::pair<As, KboWeight>... as) {
  Map<unsigned, KboWeight> out;
  __weights(out, as...);
  return out;
}
