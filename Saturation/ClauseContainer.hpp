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
 * @file ClauseContainer.hpp
 * Defines class ClauseContainer
 *
 */

#ifndef __ClauseContainer__
#define __ClauseContainer__

#include "Forwards.hpp"

#include "Lib/Event.hpp"
#include "Lib/VirtualIterator.hpp"
#include "Lib/Deque.hpp"
#include "Lib/Stack.hpp"
#include "Kernel/Clause.hpp"
#include "Lib/Set.hpp"

#include "Lib/Allocator.hpp"

#define OUTPUT_LRS_DETAILS 0

namespace Saturation
{

using namespace Lib;
using namespace Kernel;
using namespace Shell;

class ClauseContainer
{
public:
  virtual ~ClauseContainer() {}
  ClauseEvent addedEvent;
  /**
   * This event fires when a clause is removed from the
   * container because it is no longer needed, e.g. it was
   * backward-simplified, or the container is destroyed.
   * It does not fire for clauses that are removed from the
   * container because they are selected to be further
   * processed by the saturation algorithm (e.g. activated).
   */
  ClauseEvent removedEvent;
  /**
   * This event fires when a clause is removed from the
   * container to be further processed by the saturation
   * algorithm (e.g. activated).
   */
  ClauseEvent selectedEvent;
  virtual void add(Clause* c) = 0;
  void addClauses(ClauseIterator cit);
};

class RandomAccessClauseContainer
: public ClauseContainer
{
public:
  virtual void attach(SaturationAlgorithm* salg);
  virtual void detach();

  virtual unsigned sizeEstimate() const = 0;
  virtual void remove(Clause* c) = 0;
  void removeClauses(ClauseIterator cit);

protected:
  RandomAccessClauseContainer() :_salg(0) {}
  SaturationAlgorithm* getSaturationAlgorithm() { return _salg; }

  virtual void onLimitsUpdated() {}
private:
  SaturationAlgorithm* _salg;
  SubscriptionData _limitChangeSData;
};

class PlainClauseContainer : public ClauseContainer {
public:
  void add(Clause* c) override
  {
    addedEvent.fire(c);
  }
};


class UnprocessedClauseContainer
: public ClauseContainer
{
public:
  virtual ~UnprocessedClauseContainer();
  UnprocessedClauseContainer() : _data(64) {}
  void add(Clause* c) override;
  Clause* pop();
  bool isEmpty() const
  { return _data.isEmpty(); }
private:
  Deque<Clause*> _data;
};

typedef PlainEvent LimitsChangeEvent;

class PassiveClauseContainer
: public RandomAccessClauseContainer
{
public:
  PassiveClauseContainer(bool isOutermost, const Shell::Options& opt, std::string name = "") : _isOutermost(isOutermost), _opt(opt), _name(name) {}
  virtual ~PassiveClauseContainer(){};

  LimitsChangeEvent changedEvent;

  virtual bool isEmpty() const = 0;
  virtual Clause* popSelected() = 0;

  virtual unsigned sizeEstimate() const = 0;

  /*
   * LRS specific methods for computation of Limits
   */
  void updateLimits(long long estReachableCnt);

  virtual void simulationInit() = 0;
  virtual bool simulationHasNext() = 0;
  virtual void simulationPopSelected() = 0;

  // returns whether at least one of the limits was tightened
  virtual bool setLimitsToMax() = 0;
  // returns whether at least one of the limits was tightened
  virtual bool setLimitsFromSimulation() = 0;

  virtual void onLimitsUpdated() = 0;

  /*
   * LRS specific methods and fields for usage of limits
   */
  virtual bool ageLimited() const = 0;
  virtual bool weightLimited() const = 0;

  virtual bool fulfilsAgeLimit(Clause* c) const = 0;
  // note: w here denotes the weight as returned by weight().
  // age is to be recovered from inference
  // this method internally takes care of computing the corresponding weightForClauseSelection.
  virtual bool fulfilsAgeLimit(unsigned w, unsigned numPositiveLiterals, const Inference& inference) const = 0;

  virtual bool fulfilsWeightLimit(Clause* cl) const = 0;
  // note: w here denotes the weight as returned by weight().
  // age is to be recovered from inference
  // this method internally takes care of computing the corresponding weightForClauseSelection.
  virtual bool fulfilsWeightLimit(unsigned w, unsigned numPositiveLiterals, const Inference& inference) const = 0;
  
  virtual bool childrenPotentiallyFulfilLimits(Clause* cl, unsigned upperBoundNumSelLits) const = 0;

protected:
  bool _isOutermost;
  const Shell::Options& _opt;

public:
  std::string _name;
};

class ActiveClauseContainer
: public RandomAccessClauseContainer
{
public:
  ActiveClauseContainer() {}

  void add(Clause* c) override;
  void remove(Clause* c) override;

  unsigned sizeEstimate() const override { return _clauses.size(); }
  ClauseIterator clauses() const { return pvi(_clauses.iter()); }

protected:
  void onLimitsUpdated() override;
private:
  Set<Clause*> _clauses;
  // const Shell::Options& _opt;
};


};

#endif /*__ClauseContainer__*/
