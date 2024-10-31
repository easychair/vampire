/*
 * This file is part of the source code of the software program
 * Vampire. It is protected by applicable
 * copyright laws.
 *
 * This source code is distributed under the licence found here
 * https://vprover.github.io/license.html
 * and in the source directory
 */
#ifndef SUBSTITUTIONTHEORY_HPP
#define SUBSTITUTIONTHEORY_HPP

#include <algorithm>
#include <type_traits>

#include "./vector_map.hpp"
#include "./types.hpp"
#include "./log.hpp"

#include "Kernel/Term.hpp"
#include "Lib/DHMap.hpp"

namespace subsat {


/// Domain of the substitution: Vampire's FOL variables
using VampireVar = unsigned;


/// Since vampire variables aren't necessarily contiguous,
/// we map them to contiguous integers ("position") first.
using VampireVarPos = std::uint32_t;


/// Range of the substitution: Vampire's FOL terms
using VampireTerm = ::Kernel::TermList;
#define SHOWVAR(var) (VampireTerm{var, false})
#define SHOWTERM(term) (term.toString())

static_assert(sizeof(VampireTerm) == 8, "unexpected term size");


class BindingsManager final
{
public:
  // empty substitution theory
  BindingsManager()
  {
    m_bindings.reserve(32);
    m_bindings_storage.reserve(128);
  }

  // prevent accidental copies
  BindingsManager(BindingsManager& other) = delete;
  BindingsManager& operator=(BindingsManager& other) = delete;

  // but we allow moves
  BindingsManager(BindingsManager&& other) = default;
  BindingsManager& operator=(BindingsManager&& other) = default;

public:
  using BindingsEntry = std::pair<VampireVar, VampireTerm>;
  using BindingsStorage = vector<BindingsEntry>;

  class Binder final
  {
    enum class State {
      Active,
      Committed,
      Reset,
    };
  public:
    Binder(BindingsStorage& bindings_storage) noexcept
        : m_bindings_storage{bindings_storage}
        , m_bindings_start{bindings_storage.size()}
    { }

    Binder(Binder&) = delete;
    Binder& operator=(Binder&) = delete;
    Binder(Binder&&) = default;
    Binder& operator=(Binder&&) = delete;

    ~Binder() noexcept
    {
      if (m_state == State::Active) {
        reset();
      }
    }

    bool bind(VampireVar var, VampireTerm term)
    {
      ASS(m_state == State::Active);
      for (auto it = m_bindings_storage.cbegin() + m_bindings_start; it != m_bindings_storage.cend(); ++it) {
        if (it->first == var) {
          // 'var' is already bound => binding succeeds iff we bind to the same term again
          return (it->second == term);
        }
      }
      // 'var' is not bound yet => store new binding
      m_bindings_storage.emplace_back(var, term);
      return true;
    }

    void specVar(VampireVar /* var */, VampireTerm /* term */)
    {
      ASSERTION_VIOLATION;
    }

    void reset() noexcept
    {
      ASS(m_state == State::Active);
      ASS_GE(m_bindings_storage.size(), m_bindings_start);
      m_bindings_storage.resize(m_bindings_start);  // will not throw because we shrink the size
      m_state = State::Reset;
    }

    size_t size() const noexcept
    {
      ASS_GE(m_bindings_storage.size(), m_bindings_start);
      return m_bindings_storage.size() - m_bindings_start;
    }

  private:
    size_t index() const noexcept
    {
      return m_bindings_start;
    }

    /// Users should call BindingsManager::commit_bindings() instead.
    void commit() noexcept
    {
      ASS(m_state == State::Active);
      m_state = State::Committed;
    }

    friend class BindingsManager;

  private:
    BindingsStorage& m_bindings_storage;
    /// First index of the current bindings in m_bindings_storage
    size_t m_bindings_start;
    /// Current state of the binder
    State m_state = State::Active;
  };  // class Binder

  /// Only one binder may be active at any time.
  Binder start_binder() noexcept
  {
    return {m_bindings_storage};
  }

  /// Commit bindings to storage.
  void commit_bindings(Binder& binder, subsat::Var b)
  {
    binder.commit();
    while (b.index() >= m_bindings.size()) {
      m_bindings.emplace_back();
    }
    BindingsRef& bindings = m_bindings[b];
    ASS_REP(!bindings.is_valid(), "Bindings for b are already set");
    bindings.index = binder.index();
    bindings.size = binder.size();
    ASS(bindings.is_valid());
  }

  struct BindingsRef final
  {
    uint32_t index = std::numeric_limits<uint32_t>::max();
    uint32_t size;

    constexpr bool is_valid() const noexcept
    {
      return index < std::numeric_limits<uint32_t>::max();
    }

    /// last index + 1
    constexpr uint32_t end() const noexcept
    {
      return index + size;
    }
  };

  BindingsRef const& get_bindings(subsat::Var b) const
  {
    return m_bindings[b];
  }

  size_t size() const
  {
    return m_bindings.size();
  }

public:
  void clear() noexcept
  {
    m_bindings_storage.clear();
    m_bindings.clear();
  }

  bool empty() const noexcept
  {
    bool is_empty = m_bindings.empty();
    ASS(!is_empty || m_bindings_storage.empty());
    return is_empty;
  }

  void reserve(std::uint32_t bool_var_count, std::uint32_t bindings_per_bool_var, std::uint32_t vampire_var_count)
  {
    m_bindings_storage.reserve(bool_var_count * bindings_per_bool_var);
    m_bindings.reserve(bool_var_count);
  }

  std::pair<VampireVar, VampireTerm> get_single_binding(uint32_t index) {
    return m_bindings_storage[index];
  }

  friend class SubstitutionTheory;

private:
  // TODO: bindings in the array could be stored in a heap structure (i.e., one heap per binder).
  // Then we don't need to do linear search.
  // Furthermore, one could use the heap structure only for terms with a large number of variables.
  BindingsStorage m_bindings_storage;

  vector_map<subsat::Var, BindingsRef> m_bindings;
};






class SubstitutionTheory final
{
public:
  using BindingsManager = typename ::subsat::BindingsManager;
private:
  using BindingsEntry = typename BindingsManager::BindingsEntry;
  using BindingsRef = typename BindingsManager::BindingsRef;

public:
  // empty substitution theory
  SubstitutionTheory()
  {
  }

  // prevent accidental copies
  SubstitutionTheory(SubstitutionTheory& other) = delete;
  SubstitutionTheory& operator=(SubstitutionTheory& other) = delete;

  // but we allow moves
  SubstitutionTheory(SubstitutionTheory&& other) = default;
  SubstitutionTheory& operator=(SubstitutionTheory&& other) = default;

  void setBindings(BindingsManager* bm)
  {
    m_bm = bm;
    m_prepared = false;
  }

  void prepare_for_solving()
  {
    if (!m_bm || m_prepared) {
      return;
    }
    for (uint32_t b_idx = 0; b_idx < m_bm->m_bindings.size(); ++b_idx) {
      subsat::Var b{b_idx};
      BindingsRef const& bindings = m_bm->m_bindings[b];
      if (!bindings.is_valid()) {
        continue;
      }
      for (uint32_t i = bindings.index; i < bindings.end(); ++i) {
        BindingsEntry const& entry = m_bm->m_bindings_storage[i];
        VampireVar var = entry.first;
        VampireTerm term = entry.second;
        LOG_TRACE(SHOWVAR(var) << " -> " << SHOWTERM(term));
        VampireVarPos pos = m_var_pos.findOrInsert(var, m_var_pos.size());
        while (pos >= m_bindings_by_pos.size()) {
          m_bindings_by_pos.emplace_back();
        }
        m_bindings_by_pos[pos].push_back({b, term});
      }
    }
    m_prepared = true;
  }

public:
  void clear() noexcept
  {
    m_bm = nullptr;
    m_prepared = false;
    for (VampireVarPos pos = 0; pos < m_var_pos.size(); ++pos) {
      m_bindings_by_pos[pos].clear();
    }
    m_var_pos.reset();
  }

  bool empty() const noexcept
  {
    bool is_empty = !m_bm;
    if (is_empty) {
      ASS(m_var_pos.isEmpty());
#if VDEBUG
      for (auto const& vec : m_bindings_by_pos) {
        ASS(vec.empty());
      }
#endif
    }
    return is_empty;
  }

  void reserve(std::uint32_t bool_var_count, std::uint32_t bindings_per_bool_var, std::uint32_t vampire_var_count)
  {
    m_bindings_by_pos.reserve(vampire_var_count);
  }

public:
    /// Call this when a boolean variable has been set to true.
    ///
    /// The callback 'propagate' should have the following type:
    ///
    ///   PropagateCallback ~ bool propagate(subsat::Lit propagated, Lit reason)
    ///
    /// The callback will be called whenever a theory propagation is possible.
    /// It should assign the given literal and return true, if this is possible.
    /// It should return false if the assignment is conflicting.
    ///
    /// 'enable' will return true if all propagations succeeded, and false otherwise.
    template < typename PropagateCallback >
    bool enable(subsat::Var b, PropagateCallback propagate)
    {
      if (!m_bm || !m_bm->m_bindings.contains(b)) {
        return true;
      }

      // Retrieve the bindings corresponding to the given boolean variable
      BindingsRef const bindings = m_bm->m_bindings[b];

      // Exhaustively propagate the conflicting atoms
      for (std::uint32_t i = bindings.index; i < bindings.end(); ++i) {
        // Substitution constraint (X0 -> t) for b
        BindingsEntry const entry = m_bm->m_bindings_storage[i];
        VampireVar const var = entry.first;
        VampireVarPos const pos = m_var_pos.get(var);
        // Go through possible other mappings for X0
        for (auto const& q : m_bindings_by_pos[pos]) {
          if (q.second != entry.second) {
            // conflicting substitution constraints
            subsat::Var const conflicting_var = q.first;  // conflicting variable
            subsat::Lit const reason = ~b;
            if (!propagate(~conflicting_var, reason)) {
              return false;
            }
          }
        }
      }
      return true;
    }  // enable(...)

private:
  BindingsManager const* m_bm = nullptr;

  /// @brief Set to true to keep the same bindings after a call to clear().
  bool m_prepared = false;

  // Unfortunately vampire variables don't need to be contiguous,
  // so we map them to contiguous integers ("position") first.
  // TODO: replace bindings_by_pos by a single contiguous vector, and instead of VampireVarPos; use an index/length into that vector.
  //       (probably only makes sense if we construct that data lazily? in the current version we'd need an extra loop to determine counts for each variable first.)
  //       OTOH the current version isn't that bad, since we do map the vampire variables into a contiguous range first. so the vectors in m_bindings_by_pos are being reused anyway.
  Lib::DHMap<VampireVar, VampireVarPos, Kernel::IdentityHash> m_var_pos;
  vector_map<VampireVarPos, vector<std::pair<subsat::Var, VampireTerm>>> m_bindings_by_pos;
};


}  // namespace SATSubsumption

#endif /* !SUBSTITUTIONTHEORY_HPP */
