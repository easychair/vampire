/*
 * This file is part of the source code of the software program
 * Vampire. It is protected by applicable
 * copyright laws.
 *
 * This source code is distributed under the licence found here
 * https://vprover.github.io/license.html
 * and in the source directory
 */

namespace AdditionGeneralizationImpl {

using namespace std;

/**
 * Rule' 2)
 *   generalize multiplication
 *   C[X + d] 
 *   ------------ 
 *   C[X]
 *   where 
 *   - all occurences of X are in terms of the form `X + d`
 *   - sound due to substitution X -> X - d
 *   - generalization since obviously
 */

template<class NumTraits> class MonomSet;

using GenMap = Map<Variable, AnyNumber<MonomSet>, StlHash>;

template<class NumTraits>
class MonomSet 
{
  using Monom = Kernel::Monom<NumTraits>;
  using Const = typename NumTraits::ConstantType;
  using MonomFactors = Kernel::MonomFactors<NumTraits>;

  Stack<Monom> _cancellable;

  MonomSet(decltype(_cancellable) cancel) : _cancellable(cancel) {}

public:
  using Lattice = MonomSet;
  MonomSet& operator=(MonomSet&&) = default;
  MonomSet(MonomSet&&) = default;

  static MonomSet bot() 
  { return MonomSet(decltype(_cancellable){}); }

  MonomSet(Variable var, Perfect<Polynom<NumTraits>> poly) : MonomSet(decltype(_cancellable)()) 
  {
    _cancellable.reserve(poly->nSummands() - 1);
    for (auto const& monom : poly->iterSummands()) {
      if (monom.tryVar() != some(var)) {
        _cancellable.push(monom);
      }
    }
  }

  MonomSet intersect(MonomSet&& rhs) && {
    auto& lhs = *this;
    return MonomSet(intersectSortedStack(std::move(lhs._cancellable), std::move(rhs._cancellable)));
  }

  Stack<Monom> const& summands() const 
  { return _cancellable; }

  bool isBot() const 
  { return _cancellable.isEmpty(); }

  friend std::ostream& operator<<(std::ostream& out, MonomSet const& self)
  { return out << self._cancellable; }
};


struct Preprocess 
{
  GenMap& map;

  template<class NumTraits> 
  void operator()(Perfect<Polynom<NumTraits>> poly)
  {
    // a variable might occur twice within one sum.
    Set<Variable, StlHash> didOccur;
    for (auto monom : poly->iterSummands()) {
      auto var = monom.tryVar();

      if (var.isSome() && !didOccur.contains(var.unwrap())) {
        auto v = var.unwrap();
        didOccur.insert(v);
        auto gen = MonomSet<NumTraits>(v, poly);
        map.updateOrInit(v,
            [&](AnyNumber<MonomSet> old_) 
            { 
              auto old = old_.downcast<NumTraits>().unwrap();
              auto result = std::move(old).intersect(std::move(gen));
              return AnyNumber<MonomSet>(std::move(result));
            },
            [&]() { return AnyNumber<MonomSet>(std::move(gen)); });
      } else {
        for (auto factor : monom.factors->iter()) {
           if (factor.term.template is<Variable>()) {
             auto v = factor.term.template unwrap<Variable>();
             map.replaceOrInsert(v, MonomSet<NumTraits>::bot());
           }
        }
      }

    }
  }
};


/** 
 * A polymorphic closure to bottom-up evaluate clause bottom-up that replaces all occurences of the factors in the field `toRem`
 */
struct Generalize 
{
  Variable var;
  AnyNumber<MonomSet>& gen;
  bool doOrderingCheck;

  template<class NumTraits>
  Perfect<Polynom<NumTraits>> operator()(Perfect<Polynom<NumTraits>> poly, PolyNf* generalizedArgs)  
  {
    using Monom = Kernel::Monom<NumTraits>;

    auto found = poly->iterSummands()
      .find([&](Monom p) 
          { return p.tryVar() == some(var); });
    if (found.isNone()) {
      return perfect(poly->replaceTerms(generalizedArgs));
    }

    auto& toCancel = gen.downcast<NumTraits>().unwrap().summands();


    Stack<Monom> out(poly->nSummands() - toCancel.size());

    unsigned p = 0;
    unsigned genOffs = 0;

    auto pushGeneralized = [&]()  
    { 
      auto factors = perfect(poly->summandAt(p).factors->replaceTerms(&generalizedArgs[genOffs]));
      auto coeff = poly->summandAt(p).numeral;

      genOffs += factors->nFactors();
      p++;

      return out.push(Monom(coeff, factors));
    };

    auto skipGeneralized = [&]() 
    {
      genOffs += poly->summandAt(p).factors->nFactors();
      p++;
    };

    unsigned c = 0; 
    while (c < toCancel.size() && poly->summandAt(p) < toCancel[c]  ) {
      pushGeneralized();
    }
    while (p < poly->nSummands() && c < toCancel.size()) {
      if (toCancel[c] == poly->summandAt(p)) {
        skipGeneralized();
        c++;
      } else {
        ASS_L(poly->summandAt(p), toCancel[c]);
        pushGeneralized();
      }
    }
    while (p < poly->nSummands()) {
      pushGeneralized();
    }

    return perfect(Polynom<NumTraits>(std::move(out)));
      //
    // using Monom = Monom<NumTraits>;
    // auto numeral = num.downcast<NumTraits>().unwrap();
    // auto newFactors = perfect(monom.factors->replaceTerms(evaluatedArgs));
    // for (auto f : monom.factors->iter()) {
    //   if (f.tryVar() == Option<Variable>(var)) {
    //     ASS_EQ(numeral, monom.numeral)
    //     return Monom(decltype(numeral)(1), newFactors);
    //   }
    // }
    // return Monom(monom.numeral, newFactors);
  }
};


struct IsBot 
{
  template<class C>
  bool operator()(C const& lattice)
  { return lattice.isBot(); }
};

/** applies the rule */ 
SimplifyingGeneratingInference1::Result applyRule(Clause* cl, bool doOrderingCheck) 
{
  DEBUG("input clause: ", *cl);

  GenMap map;

  for (auto poly : iterPolynoms(cl)) {
    poly.apply(Preprocess {map});
  }

  Option<typename GenMap::Entry &> selected;
  for (auto& e : iterTraits(map.iter()) ) {
    if (!e.value().apply(IsBot{})) {
      /* we use the entry with the least variable, in order to prevent non-determinism */
      if (selected.isNone() || e.key() < selected.unwrap().key()) {
        selected = decltype(selected)(e);
      }
    }
  }

  if (selected.isNone()) {
    DEBUG("not applicable")
    return SimplifyingGeneratingInference1::Result::nop(cl);
  } else {
    auto& e = selected.unwrap();
    DEBUG("selected generalization: ", e.key(), " ", e.value());
    Generalize gen { e.key(), e.value(), doOrderingCheck };
    return generalizeBottomUp(cl, EvaluatePolynom<Generalize> {gen});
  }

}

} // namespace AdditionGeneralizationImpl
