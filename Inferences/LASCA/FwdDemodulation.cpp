#include "Inferences/LASCA/FwdDemodulation.hpp"
#include "Saturation/SaturationAlgorithm.hpp"

#define DEBUG(...)  // DBG(__VA_ARGS__)
using Demod = Inferences::LASCA::Demodulation;

////////////////////////////////////////////////////////////////////////////////////////////////////
// INDEXING
////////////////////////////////////////////////////////////////////////////////////////////////////

namespace Inferences {
namespace LASCA {

#if VDEBUG
void FwdDemodulation::setTestIndices(Stack<Indexing::Index*> const& indices) 
{
  _index = (decltype(_index)) indices[0]; 
  _index->setShared(_shared);
}
#endif

void FwdDemodulation::attach(SaturationAlgorithm* salg)
{  
  ASS(!_index);

  this->ForwardSimplificationEngine::attach(salg);
  _index=static_cast<decltype(_index)> (
	  _salg->getIndexManager()->request(LASCA_FWD_DEMODULATION_SUBST_TREE));
  _index->setShared(_shared);
}

void FwdDemodulation::detach()
{

  CALL("Superposition::detach");
  ASS(_salg);

  _index=0;
  _salg->getIndexManager()->release(LASCA_FWD_DEMODULATION_SUBST_TREE);
  this->ForwardSimplificationEngine::detach();
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// RULE APPLICATION
////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * Perform forward simplification on @b cl
 *
 * Return true if the simplification is applicable on @b cl,
 * set @b replacement to a replacing clause if there is one (otherwise keep @b replacement = nullptr)
 *
 * @b premises will contain clauses that justify the simplification
 * performed.
 */
bool FwdDemodulation::perform(Clause* toSimplify, Clause*& replacement, ClauseIterator& premises)
{
  ASS_EQ(replacement, NULL)
  Stack<Literal*> simplified;
  for (auto rhs : Rhs::iter(*_shared, toSimplify)) {
    // DEBUG("simplifyable position: ", pos.term, " in ", *pos.lit)
    for (auto lhs : _index->generalizations(rhs.term)) {
      auto simplified = Demodulation::apply(*_shared, *lhs.data, rhs);
      if (simplified.isSome()) {
        replacement = simplified.unwrap();
        premises    = pvi(getSingletonIterator(lhs.data->clause()));
        return true;
      }
    }
  }
  premises = ClauseIterator::getEmpty();
  return false;
}

} // namespace LASCA
} // namespace Inferences
