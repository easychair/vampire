/*
 * This file is part of the source code of the software program
 * Vampire. It is protected by applicable
 * copyright laws.
 *
 * This source code is distributed under the licence found here
 * https://vprover.github.io/license.html
 * and in the source directory
 */

#include "Lib/Coproduct.hpp"
#include "Test/UnitTesting.hpp"

using namespace std;
using namespace Kernel;
using namespace Lib;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// USAGE EXAMPLES
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TEST_FUN(examples__is_variant_01) {
  // a Coproduct can be initialized with different variants
  auto x = Coproduct<int, int, float>::variant<2>(1.0f);
 
  // one can check which variant we have with the template function is
  ASS(!x.is<0>())
  ASS(!x.is<1>())
  ASS(x.is<2>())
  // x.is<4>() // <- does not compile
 
  // we can get the value using unwrap
  ASS_EQ(1.0f, x.unwrap<2>());

}


TEST_FUN(examples__is_variant_02) {
  // for coproducts where all variants are of distinct types we can leave away the variand index for construction
  auto x = Coproduct<int, float>(1.0f);
 
  // and in this case we can also use the type for is and unrap
  ASS(x.is<float>())
  ASS(!x.is<int>())
  ASS_EQ(x.unwrap<float>(), 1.0f);

  // x.is<char*>() // <- does not compile
}


TEST_FUN(examples__is_variant_03) {
  auto const x = Coproduct<int, float>(1.0f);
 
  // `as` is the combined version of `is` and `unwrap`:
  ASS(x.as<float>().isSome())
  ASS(x.as<int>().isNone())

  ASS_EQ(x.as<float>().toOwned(), Option<float>(1.0f));
  ASS_EQ(x.as<int>()  .toOwned(), Option<int>());

  // x.is<char*>() // <- does not compile
}

TEST_FUN(examples__equal_01) {
  // two corproduct variants can have the same type and value, but a distinct tag.
  auto x = Coproduct<int, int, float>::variant<0>(0);
  auto y = Coproduct<int, int, float>::variant<1>(0);
  ASS(x != y)
}

TEST_FUN(examples__equal_02) {
  // they are only equal if their tag and their content matches
  auto x = Coproduct<int, int, float>::variant<0>(0);
  auto y = Coproduct<int, int, float>::variant<0>(0);
  ASS(x == y)
}

TEST_FUN(examples__equal_03) {
  auto x = Coproduct<int, int, float>::variant<0>(0);
  auto y = Coproduct<int, int, float>::variant<0>(1);
  ASS(x != y)
}

TEST_FUN(examples__match_01) {
  // we can also transform the content using the match method
  auto x = Coproduct<int, float>(1);
  auto isGreaterThanZero = x.match(
      [](int   i) { return i > 0;  },
      [](float f) { return f > 0.0f; }
  );
  ASS(isGreaterThanZero)
}

TEST_FUN(examples__match_02) {
  // Further we can create polymorphic function structs if each match branch does the same thing
  auto x = Coproduct<int, float>(1);

  std::string str = x.apply([](auto const& c) {
    std::stringstream out;
    out << c;
    return out.str(); 
  });
  ASS_EQ(str, "1")
}

TEST_FUN(examples__compare) {
  // Coproducts are orderd first by tag, then by value.
  using Co = Coproduct<int, double>;
  ASS(Co(1) < Co(1.0))
  ASS(Co(2) < Co(1.0))
  ASS(Co(2) < Co(3))
  ASS(Co(1.0) < Co(2.0))
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MISC TESTS
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TEST_FUN(unwrap_01) {
  auto x = Coproduct<int, int, float>::variant<0>(0);
  auto y = Coproduct<int, int, float>::variant<1>(0);
  ASS(x.unwrap<0>() == y.unwrap<1>());
}


TEST_FUN(unwrap_02) {
  auto x = Coproduct<int, int, float>::variant<0>(0);
  auto y = Coproduct<int, int, float>::variant<1>(1);
  ASS(x.unwrap<0>() != y.unwrap<1>());
}


struct NonCopy {
  bool content;
  bool wasMoved;
  NonCopy(bool content) : content(content), wasMoved(false) {}
  NonCopy(NonCopy&& other) 
    : content(other.content)
    , wasMoved(false)  {
    other.wasMoved = true;
  }
  NonCopy& operator=(NonCopy&& other) {
    content = other.content;
    other.wasMoved = true;
    return *this;
  }

  bool operator==(const NonCopy& other) const {
    return content == other.content;
  }
  friend std::ostream& operator<<(std::ostream& out, const NonCopy& x)  {
    return out << "NonCopy(" << x.content << ")";
  }
};

TEST_FUN(move_01) {

  auto y = Coproduct<int, NonCopy>::variant<1>(NonCopy( true ));

  ASS((y == Coproduct<int,NonCopy>::variant<1>(NonCopy( true ))));
  ASS((y != Coproduct<int,NonCopy>::variant<1>(NonCopy( false ))));
  ASS((y != Coproduct<int,NonCopy>::variant<0>(1)));

  y = Coproduct<int, NonCopy>::variant<0>(1);

  ASS((y == Coproduct<int,NonCopy>::variant<0>(1)));
  ASS((y != Coproduct<int,NonCopy>::variant<1>(NonCopy( false ))));
  ASS((y != Coproduct<int,NonCopy>::variant<0>(0)));
}

TEST_FUN(apply01) {

  struct Bar {
    int f() { return 42; }
  } b;

  struct Foo {
    int f() { return -42; }
  } f;


  using C = Coproduct<Foo,Bar>;
  auto foo = C(f);
  auto& fooRef = foo;
  auto bar = C(b);

  auto applyF = [](auto& x) { return x.f(); };
  auto& applyFref = applyF;
  static_assert( std::is_same<decltype(foo.template unwrap<0>()), Foo& >::value, "");
  static_assert( std::is_same<decltype(foo.template unwrap<1>()), Bar& >::value, "");
  static_assert( std::is_same<decltype(applyF((foo.template unwrap<1>()))), int >::value, "");
  static_assert( std::is_same<decltype(applyF((foo.template unwrap<0>()))), int >::value, "");
  static_assert( std::is_same<decltype(applyFref((fooRef.template unwrap<0>()))), int >::value, "");

  auto doApply = [](auto f) {

    auto foo = C(Foo{});
    auto doApplication = [&]() {
      return f(foo.template unwrap<0>());
    };
    return doApplication();
  };
  ASS_EQ(doApply(applyF), -42)

  ASS_EQ( 42, bar.apply([](auto& x) { 
  static_assert( std::is_same<decltype(x), Foo& >::value
              || std::is_same<decltype(x), Bar& >::value, "");
        return x.f(); }))
  ASS_EQ(-42, foo.apply([](auto& x) { return x.f(); }))

}

// TEST_FUN(move_02) {
//
//   auto x = Coproduct<int, NonCopy>::variant<1>(NonCopy( true ));
//
//   auto y = x;
//
//   ASS(!x.unwrap<NonCopy>().wasMoved)
//   ASS(!y.unwrap<NonCopy>().wasMoved)
//
// }
