#ifndef SUBSAT_HPP
#define SUBSAT_HPP

#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <new>
#include <vector>

// Ensure NDEBUG and VDEBUG are synchronized
#ifdef NDEBUG
static_assert(!VDEBUG, "");
#else
static_assert(VDEBUG, "");
#endif

// TODO:
// Once this works, make a separate version 'matchsat',
// which keeps an array of matches as well.
// (see my notes on SAT+CSP)

namespace SMTSubsumption {

using std::uint32_t;

enum class Value : signed char {
  False = -1,
  Unknown = 0,
  True = 1,
};

class Lit;

/// Boolean variable represented by its integer index.
/// Use consecutive indices starting at 0.
class Var final {
  uint32_t m_index;

public:
  explicit constexpr Var(uint32_t index) noexcept
      : m_index{index}
  {
    // assert(m_index <= Var::max_index());  // TODO: how to assert in constexpr constructor?
  }

  [[nodiscard]] constexpr uint32_t index() const noexcept
  {
    return m_index;
  }

  [[nodiscard]] static constexpr uint32_t max_index() noexcept
  {
    return (1u << 31) - 2;
  }

  [[nodiscard]] static constexpr Var invalid() noexcept
  {
    return Var{std::numeric_limits<uint32_t>::max()};
  }

  // [[nodiscard]] constexpr Lit operator~() const noexcept
  // {
  //   return Lit{*this, false};
  // }

  // [[nodiscard]] constexpr operator Lit() const noexcept
  // {
  //   return Lit{*this, true};
  // }
}; // Var

static_assert(Var::max_index() == static_cast<uint32_t>(INT32_MAX - 1), "unexpected max variable index");
static_assert(Var::max_index() < Var::invalid().index(), "valid variable indices overlap with invalid sentinel value");

/// Boolean literals represented by integer index.
/// The least significant bit indicates the sign.
///
/// Mapping from variable indices to literal indices:
///    Lit{0} ... 0
///   ~Lit{0} ... 1
///    Lit{1} ... 2
///   ~Lit{1} ... 3
///      :
///      :
class Lit final {
  uint32_t m_index;

private:
  friend class Clause;
  /// Uninitialized value (for clause constructor)
  Lit() noexcept
  // : m_index{Lit::invalid().index()}
  {
  }

  explicit constexpr Lit(uint32_t index) noexcept
      : m_index{index}
  {
    // assert(m_index <= Lit::max_index()); // TODO: how to assert in constexpr constructor?
  }

public:
  explicit constexpr Lit(Var var, bool positive) noexcept
      : Lit{2 * var.index() + static_cast<uint32_t>(!positive)}
  {
  }

  [[nodiscard]] static constexpr Lit pos(Var var) noexcept
  {
    return Lit{var, true};
  }

  [[nodiscard]] static constexpr Lit neg(Var var) noexcept
  {
    return Lit{var, false};
  }

  [[nodiscard]] constexpr uint32_t index() const noexcept
  {
    return m_index;
  }

  [[nodiscard]] static constexpr uint32_t max_index() noexcept
  {
    static_assert(Var::max_index() < (std::numeric_limits<uint32_t>::max() - 1) / 2, "cannot represent all literals");
    return 2 * Var::max_index() + 1;
  }

  [[nodiscard]] static constexpr Lit invalid() noexcept
  {
    return Lit{std::numeric_limits<uint32_t>::max()};
  }

  [[nodiscard]] constexpr bool is_positive() const noexcept
  {
    return (m_index & 1) == 0;
  }

  [[nodiscard]] constexpr bool is_negative() const noexcept
  {
    return !is_positive();
  }

  [[nodiscard]] constexpr Lit operator~() const noexcept
  {
    return Lit{m_index ^ 1};
  }
}; // Lit

static_assert(Lit::max_index() < Lit::invalid().index(), "valid literal indices overlap with invalid sentinel value");

inline void* subsat_alloc(std::size_t size)
{
#ifdef SUBSAT_STANDALONE
  void* p = std::malloc(size);
#else
  void* p = ALLOC_UNKNOWN(size, "subsat");
#endif
  if (!p && size > 0) {
    throw std::bad_alloc();
  }
  return p;
}

inline void subsat_dealloc(void* p)
{
#ifdef SUBSAT_STANDALONE
  std::free(p);
#else
  DEALLOC_UNKNOWN(p, "subsat");
#endif
}






class Clause final {
  uint32_t m_size;   // number of literals
  Lit m_literals[2]; // actual size is m_size, but C++ does not officially support flexible array members (as opposed to C)

public:
  using iterator = Lit const*;

  [[nodiscard]] iterator begin() const noexcept
  {
    return &m_literals[0];
  }

  [[nodiscard]] iterator end() const noexcept
  {
    return begin() + m_size;
  }

  Lit& operator[](size_t idx) noexcept
  {
    assert(idx < m_size);
    return m_literals[idx];
  }

  Lit const& operator[](size_t idx) const noexcept
  {
    assert(idx < m_size);
    return m_literals[idx];
  }

  /// Number of bytes required by a clause containing 'size' literals.
  static size_t bytes(uint32_t size) noexcept
  {
    size_t const embedded_literals = std::extent_v<decltype(m_literals)>;
    size_t const additional_literals = (size >= embedded_literals) ? (size - embedded_literals) : 0;
    size_t const total_bytes = sizeof(Clause) + sizeof(Lit) * additional_literals;
    return total_bytes;
  }

  /// Allocate a clause with enough space for 'size' literals.
  static Clause* create(uint32_t size)
  {
    void* p = subsat_alloc(bytes(size));
    return new (p) Clause{size};
  }

  // static void* operator new(size_t, uint32_t num_literals)
  // {
  //   size_t const contained_literals = std::extent_v<decltype(m_literals)>;
  //   size_t const additional_literals = std::max(0, static_cast<size_t>(num_literals) - contained_literals);
  //   size_t const total_bytes = sizeof Clause + sizeof Lit * additional_literals;
  //   return ::operator new(total_bytes);
  // }

private:
  // NOTE: do not use this constructor directly
  // because it does not allocate enough memory for the literals
  Clause(uint32_t size) noexcept
      : m_size{size}
  {
  }

  // cannot copy/move because of flexible array member
  Clause(Clause&) = delete;
  Clause(Clause&&) = delete;
  Clause& operator=(Clause&) = delete;
  Clause& operator=(Clause&&) = delete;
}; // Clause

template <typename Key>
class IndexMember {
  typename std::invoke_result_t<typename Key::index> operator()(Key key) const
  {
    return key.index();
  }
};

template <typename Integer>
class IndexIdentity {
  Integer operator()(Integer key) const noexcept
  {
    return key;
  }
};

template <typename Key> struct DefaultIndex;
template <> struct DefaultIndex<Var> {
  using type = IndexMember<Var>;
};
template <> struct DefaultIndex<Lit> {
  using type = IndexMember<Lit>;
};
template <> struct DefaultIndex<uint32_t> {
  using type = IndexIdentity<uint32_t>;
};

// template <typename Key, typename T, typename Index = IndexMember<Key>>
template <typename Key, typename T, typename Index = DefaultIndex<Key>>
class ivector {
public:
  using key_type = Key;
  using value_type = T;
  using reference = value_type&;
  using const_reference = value_type const&;
  using size_type = typename std::vector<T>::size_type;

  reference operator[](key_type key) { return m_data[index(key)]; }
  const_reference operator[](key_type key) const { return m_data[index(key)]; }

  void reserve(size_type new_cap) { m_data.reserve(new_cap); }

private:
  size_type index(Key key) const
  {
    Index index;
    return index(key);
  }

private:
  std::vector<T> m_data;
};





class Solver {
public:
  using ClauseRef = uint32_t;
  using Level = uint32_t;

  /// Ensure space for a new variable and return it.
  [[nodiscard]] Var new_variable()
  {
    reserve_variables(1);
    return Var{m_used_variables++};
  }

  /// Reserve space for n additional variables.
  void reserve_variables(uint32_t count)
  {
    if (m_used_variables + count <= m_allocated_variables) {
      return;
    }
    m_allocated_variables = m_used_variables + count;
    m_values.reserve(2 * m_allocated_variables);
    // TODO
  }

  void add_empty_clause()
  {
  }

  void add_unit_clause(Lit* lit)
  {
  }

  void add_binary_clause(Lit* lit1, Lit* lit2)
  {
  }

  void add_clause(Clause* clause)
  {
  }

private:
  void assign(Lit lit, Value value)
  {
  }

  void propagate_units()
  {
  }

  void propagate()
  {
  }

  void analyze()
  {
  }

private:
  bool m_inconsistent;
  uint32_t m_used_variables;
  uint32_t m_allocated_variables;
  std::vector<Lit> m_units;
  ivector<Lit, ClauseRef> m_reasons;
  ivector<Lit, Value> m_values; ///< Current assignment of literals
  ivector<ClauseRef, Clause*> m_clauses;
  // std::vector<Clause*> clauses;
}; // Solver



} // namespace SMTSubsumption

#endif /* !SUBSAT_HPP */
