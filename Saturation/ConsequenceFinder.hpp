/**
 * @file ConsequenceFinder.hpp
 * Defines class ConsequenceFinder for splitting with backtracking.
 */

#ifndef __ConsequenceFinder__
#define __ConsequenceFinder__

#include "../Forwards.hpp"

#include "../Lib/Allocator.hpp"
#include "../Lib/Array.hpp"
#include "../Lib/Event.hpp"
#include "../Lib/Int.hpp"

#include "../Kernel/Clause.hpp"

#include "../Inferences/TautologyDeletionISE.hpp"

namespace Saturation {

using namespace Kernel;
using namespace Inferences;

class ConsequenceFinder {
public: 
  ~ConsequenceFinder();

  void init(SaturationAlgorithm* sa);

  bool isRedundant(Clause* cl);

  void onNewPropositionalClause(Clause* cl);

private:

  void onClauseInserted(Clause* cl);
  void onClauseRemoved(Clause* cl);

  void indexClause(unsigned indexNum, Clause* cl, bool add);

  SaturationAlgorithm* _sa;

  typedef SkipList<Clause*,Int> ClauseSL;

  ZIArray<ClauseSL*> _index;
  ZIArray<bool> _redundant;

  TautologyDeletionISE _td;

  /** SubscriptionData for the @b onClauseInserted method */
  SubscriptionData _sdInsertion;
  /** SubscriptionData for the @b onClauseRemoved method */
  SubscriptionData _sdRemoval;
};

}

#endif // __ConsequenceFinder__
