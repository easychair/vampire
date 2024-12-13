/*
 * This file is part of the source code of the software program
 * Vampire. It is protected by applicable
 * copyright laws.
 *
 * This source code is distributed under the licence found here
 * https://vprover.github.io/license.html
 * and in the source directory
 */

#ifndef STL_HPP
#define STL_HPP

#include <algorithm>
#include <memory>
#include <type_traits>
#include <utility>

namespace Lib {

template<class T>
T* move_to_heap(T t) { return new T(std::move(t)); }
                
template< typename T, 
          std::enable_if_t<!std::is_pointer<T>::value, bool> = true  
        >
std::shared_ptr<T> make_shared(T t)
{ return std::shared_ptr<T>(move_to_heap(std::move(t))); }

}  // namespace Lib

#endif /* !STL_HPP */
