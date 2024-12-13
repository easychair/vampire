/*
 * This file is part of the source code of the software program
 * Vampire. It is protected by applicable
 * copyright laws.
 *
 * This source code is distributed under the licence found here
 * https://vprover.github.io/license.html
 * and in the source directory
 */

namespace VariableMultiplicationGeneralizationImpl {

using namespace std;

/**
*  Rule 3)
*    generalize variable multiplication
*    C[X0 ⋅ X1 ⋅ ... ⋅ Xn] 
*    ------------ 
*    C[X0]
*    where 
*    - all occurences of Xi are in terms of the form `X0 ⋅ X1 ⋅ ... ⋅ Xn`
*    - X0 /= Xi (for i /= 0)
*    - all Xi are distinct variables
*    - sound due to substitution { X1 -> 1 ..., XN -> 1 }
*    - obviously a generalization 
*/

/** 
 * Type for associating objects with integer ids. It is mainly only used in order to use IntUnionFind with other types than int.
 */
template<class A>
class IntMap
{
  Map<A, unsigned, StlHash> _map;
  Stack<A> _stack;
public:
  unsigned insert(A obj) 
  {
    return _map.getOrInit(std::move(obj), [&](){ 
      auto idx = _stack.size();
      _stack.push(obj);
      return idx; 
    });
  }

  unsigned size() const 
  { return _stack.size(); }

  unsigned toInt(A const& obj) const
  { return _map.get(obj); }

  A const& fromInt(int obj)  const
  { return _stack[obj]; }
};



/**
 * Represents the region of variables `X0 ⋅ X1 ⋅ ... ⋅ Xn` from the rule
 *
 * Two regions can be intersected. 
 */
class VariableRegion 
{
  Coproduct<
    Stack<AnyNumber<MonomFactor>>, /* <- always sorted; contains only factors of the form `variable ^ n` */
    Top /* <- unininitialized */
      > _self;

public:
  VariableRegion() : _self(Top{}) {}
  VariableRegion(Stack<AnyNumber<MonomFactor>>&& self) : _self(self) {}
  VariableRegion(VariableRegion &&) = default;
  VariableRegion& operator=(VariableRegion &&) = default;

  /* intersects two regions */
  VariableRegion meet(VariableRegion rhs) 
  {
    auto& lhs = *this;
    if (lhs.isUninit()) return VariableRegion(std::move(rhs));
    if (rhs.isUninit()) return VariableRegion(std::move(lhs));
    return VariableRegion(intersectSortedStack(std::move(lhs.unwrap()), std::move(rhs.unwrap())));
  }

  Stack<AnyNumber<MonomFactor>> const& unwrap() const 
  { return _self.template unwrap<Stack<AnyNumber<MonomFactor>>>(); }

  Stack<AnyNumber<MonomFactor>> & unwrap() 
  { return _self.template unwrap<Stack<AnyNumber<MonomFactor>>>(); }

  friend std::ostream& operator<<(std::ostream& out, VariableRegion const& self) 
  {
    return self.isUninit() ? out << "Top"
                        : out << self.unwrap();
  }

  bool isUninit() const 
  { return _self.template is<Top>(); }

  bool isInit() const 
  { return !isUninit(); }
};


/**
 * Polymorphic closure that processes each subterm of the clause. 
 *
 * All variables that occur together in a product (i.e. a monom), are being associated with the same "connected component" 
 * using the IntUnionFind components.
 *
 * For each of these components a minimal VariableRegion is kept stored in the field varRegions.
 */
struct Preprocess 
{
  IntUnionFind& components;
  IntMap<Variable> &varMap;
  Stack<VariableRegion> &varRegions;

  VariableRegion& varSet(int v)
  { return varRegions[v]; }

  int root(Variable v) const
  { return components.root(varMap.toInt(v)); }

  template<class NumTraits> 
  void operator()(Perfect<Polynom<NumTraits>> p) 
  {
    for (auto summand : p->iterSummands()) {

      auto varIter = summand.factors->iter()
            .filter([](MonomFactor<NumTraits> factor) { return factor.term.template is<Variable>(); });

      auto varIter2 = varIter;
      auto varStack = VariableRegion(
          varIter2
            .map([](MonomFactor<NumTraits> factor) { return AnyNumber<MonomFactor>(factor); })
            .template collect<Stack>());

      auto vars = varIter.map([](MonomFactor<NumTraits> factor) { return factor.term.template unwrap<Variable>(); });

      if (vars.hasNext())  {
        auto cur = root(vars.next());

        varSet(cur) = std::move(varSet(cur)).meet(std::move(varStack));

        for (auto var : vars) {
          cur = joinRegions(cur, root(var));
        }

      }
    }
  }

  void dbgState() const {
    DEBUG("---------------------");
    for (unsigned i = 0; i < varMap.size(); i++) {
      DEBUG(varMap.fromInt(i), " -> ", varMap.fromInt(components.root(i)), " -> ", varRegions[components.root(i)]);
    }
    DEBUG("---------------------");
  }

  int joinRegions(int v, int w)
  {
    if (v == w) return v;

    components.doUnion(v,w);
    auto r = components.root(v);
    varSet(r) = std::move(varSet(v)).meet(std::move(varSet(w)));

    return r;
  }

};

/** 
 * A polymorphic closure to bottom-up evaluate clause bottom-up that replaces all occurences of the factors in the field `toRem`
 */
struct Generalize 
{
  Stack<AnyNumber<MonomFactor>> const& toRem; /* <- always expected to be sorted */
  bool doOrderingCheck;

  template<class NumTraits>
  Monom<NumTraits> operator()(Monom<NumTraits> p, PolyNf* evaluatedArgs)  
  {
    using Pair = Monom<NumTraits>;
    return Pair(p.numeral, perfect(MonomFactors<NumTraits>(filter(p.factors, evaluatedArgs))));
  }

  template<class NumTraits>
  Stack<MonomFactor<NumTraits>> filter(Perfect<MonomFactors<NumTraits>> const& factors, PolyNf* evaluatedArgs)
  {
    Stack<MonomFactor<NumTraits>> out;
    unsigned rmI = 0;
    unsigned m = 0;

    auto skip = [&]() { rmI++; m++; };
    auto push = [&]() { out.push(MonomFactor<NumTraits>(evaluatedArgs[m], factors->factorAt(m).power)); m++; };


    while (m < factors->nFactors() && rmI < toRem.size()) {
      auto factor = factors->factorAt(m);
      auto rm = toRem[rmI].template downcast<NumTraits>();
      if (rm.isNone()) {
        push();
      } else if (factor == rm.unwrap()) {
        skip();
      } else if (factor < rm.unwrap()) {
        push();
      } else {
        ASS_L(rm.unwrap(), factor)
        rmI++;
      }
    }
    while (m < factors->nFactors()) {
      push();
    }

    std::sort(out.begin(), out.end());
    return out;
  }
};

/** 
 * applies the rule
 */ 
SimplifyingGeneratingInference1::Result applyRule(Clause* cl, bool doOrderingCheck) 
{
  DEBUG("input clause: ", *cl);
  IntMap<Variable> varMap;

  /* initialization */
  for (auto var : iterVars(cl)) {
    varMap.insert(var);
  }
  if (varMap.size() == 0) {
    DEBUG("no variables. generalization not applicable");
    return SimplifyingGeneratingInference1::Result::nop(cl);
  }

  IntUnionFind components(varMap.size());
  Stack<VariableRegion> varRegions(varMap.size());;
  for (unsigned i = 0; i < varMap.size(); i++)  {
    varRegions.push(VariableRegion());
  }

  /* preprocessing. finds all products `X0 ⋅ X1 ⋅ ... ⋅ Xn` such that the rule is applicable */
  for (auto poly : iterPolynoms(cl)) {
    poly.apply(Preprocess {components, varMap, varRegions});
  }


  /* create a stack of all variables that shall be removed in the final step */

  Stack<AnyNumber<MonomFactor>> remove;

  components.evalComponents();
  for (auto comp : iterTraits(IntUnionFind::ComponentIterator(components))) {
    auto& maybeRegion = varRegions[components.root(comp.next())];
    if (maybeRegion.isInit()) {
      auto& region = maybeRegion.unwrap();

      /* one variable with power one needs to be kept, per varible region */
      auto var = iterTraits(region.iter())
        .map([](auto& x) { return x; })
        .filter([](auto& p) { return p.apply([](auto& t){ return t.tryVar(); }).isSome(); })
        .tryNext();

      if (var.isSome()) {
        for (auto varPower : region) {
          if (varPower != var.unwrap()) {
            remove.push(varPower);
          }
        }
      }
    }
  }

  /* apply the substitution `X0 ⋅ X1 ⋅ ... ⋅ Xn ==> X0`  */
  DEBUG("removing variables: ", remove)
  if (remove.isEmpty()) {
    return SimplifyingGeneratingInference1::Result::nop(cl);
  } else {
    std::sort(remove.begin(), remove.end());
    Generalize gen { remove, doOrderingCheck };
    return generalizeBottomUp(cl, EvaluateMonom<Generalize> {gen});
  }
}

} // namespace VariableMultiplicationGeneralizationImpl 


