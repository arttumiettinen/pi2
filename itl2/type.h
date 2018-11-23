#pragma once
/*

Copyright (c) 2018, ITHare.com
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef ithare_util_type_type_h_included
#define ithare_util_type_type_h_included

#include <cstdint>
#include <type_traits>
#ifdef ITHARE_UTIL_INTERNAL_DEBUG
#include <assert.h> internal debug only
#define ITHARE_UTIL_ASSERT(x) assert(x)
#else
#define ITHARE_UTIL_ASSERT(x)
#endif

//namespace ithare { namespace util { namespace type {
 
namespace intuitive {
  using MAX_EFFICIENT_INT = intptr_t;//platform-dependent, but intptr_t is not a bad starting point
									 //  if it is suboptimal for some platform - open an issue, we'll add an #ifdef for that platform 

//DUPLICATED CODE (1/2). It is obviously possible to get rid of duplicates, 
//                       but it is currently unclear how it may affect performance of the compiled code, so for now we prefer it this way 
template<class TA, class TB>
inline constexpr bool lt(TA a, TB b) {
  static_assert(std::is_integral<TA>::value);
  static_assert(std::is_integral<TB>::value);
  constexpr bool aSigned = std::is_signed<TA>::value;
  constexpr bool bSigned = std::is_signed<TB>::value;
  if constexpr(aSigned == bSigned)
    return a < b;//both signed or both unsigned - no casts required, C promotions will do just fine
  else {//different is_signed, let's make TSIGNED always-signed, and TUNSIGNED - always-unsigned
    using TSIGNED = typename std::conditional<aSigned,TA,TB>::type;
    using TUNSIGNED = typename std::conditional<aSigned,TB,TA>::type;
 
    static_assert(sizeof(TSIGNED)+sizeof(TUNSIGNED)==sizeof(TA)+sizeof(TB));//self-check
    if constexpr(sizeof(TSIGNED)>sizeof(TUNSIGNED))
      return a < b;//again, no casts required, C promotions will do just fine (promoting b to TA which is signed)
    if constexpr(sizeof(TUNSIGNED)<sizeof(int)) {
      return a < b;//again, no casts required, C promotion-to-int will do just fine (promoting both to int which is signed) 
    } 
     
    //at this point, we have sizeof(TUNSIGNED) >= sizeof(TSIGNED) => no-cast will be counterintuitive
    if constexpr(sizeof(TUNSIGNED)<sizeof(MAX_EFFICIENT_INT)) {
      //we can cast unsigned to MAX_EFFICIENT_INT
      if constexpr(aSigned) {
        ITHARE_UTIL_ASSERT(!bSigned);
        return a < MAX_EFFICIENT_INT(b);
      }
      else {
        ITHARE_UTIL_ASSERT(bSigned);
        return MAX_EFFICIENT_INT(a) < b; 
      } 
     
    }
    else { //last resort: expression 
	  ITHARE_UTIL_ASSERT(sizeof(TUNSIGNED)>=sizeof(TSIGNED));
      if constexpr(aSigned) {
        //return a<0 ? true : TUNSIGNED(a) < b;
        return (a<0) | (TUNSIGNED(a) < b);//sic - single '|', seems to perform better under GCC (avoids branch)
	  }
      else {
        ITHARE_UTIL_ASSERT(bSigned);
        //return b<0 ? false : a < TUNSIGNED(b);
        return (b>=0) & (a < TUNSIGNED(b));//sic - single '&', seems to perform better under GCC (avoids branch)
      }
    }
  }
}

//DUPLICATED CODE (2/2)
template<class TA, class TB>
inline constexpr bool eq(TA a, TB b) {
  static_assert(std::is_integral<TA>::value);
  static_assert(std::is_integral<TB>::value);
  constexpr bool aSigned = std::is_signed<TA>::value;
  constexpr bool bSigned = std::is_signed<TB>::value;
  if constexpr(aSigned == bSigned)
    return a == b;//both signed or both unsigned - no casts required, C promotions will do just fine
  else {//different is_signed, let's make TSIGNED always-signed, and TUNSIGNED - always-unsigned
    using TSIGNED = typename std::conditional<aSigned,TA,TB>::type;
    using TUNSIGNED = typename std::conditional<aSigned,TB,TA>::type;
 
    static_assert(sizeof(TSIGNED)+sizeof(TUNSIGNED)==sizeof(TA)+sizeof(TB));//self-check
    if constexpr(sizeof(TSIGNED)>sizeof(TUNSIGNED))
      return a == b;//again, no casts required, C promotions will do just fine (promoting b to TA which is signed)
    if constexpr(sizeof(TUNSIGNED)<sizeof(int)) {
      return a == b;//again, no casts required, C promotion-to-int will do just fine (promoting both to int which is signed) 
    } 
     
    //at this point, we have sizeof(TUNSIGNED) >= sizeof(TSIGNED) => no-cast will be counterintuitive
    if constexpr(sizeof(TUNSIGNED)<sizeof(MAX_EFFICIENT_INT)) {
      //we can cast unsigned to MAX_EFFICIENT_INT
      if constexpr(aSigned) {
        ITHARE_UTIL_ASSERT(!bSigned);
        return a == MAX_EFFICIENT_INT(b);
      }
      else {
        ITHARE_UTIL_ASSERT(bSigned);
        return MAX_EFFICIENT_INT(a) == b; 
      } 
    }
    else { //last resort: expression
	  ITHARE_UTIL_ASSERT(sizeof(TUNSIGNED)>=sizeof(TSIGNED));
      if constexpr(aSigned)
        return a<0 ? false : TUNSIGNED(a) == b;
      else {
        ITHARE_UTIL_ASSERT(bSigned);
        return b<0 ? false : a == TUNSIGNED(b);
      }
    }
  }
}

template<class TA, class TB>
inline constexpr bool gt(TA a, TB b) {
	return lt(b,a);
}

template<class TA, class TB>
inline constexpr bool le(TA a, TB b) {
	return !lt(b,a);
}

template<class TA, class TB>
inline constexpr bool ge(TA a, TB b) {
	return !lt(a,b);
}

template<class TA, class TB>
inline constexpr bool ne(TA a, TB b) {
	return !eq(a,b);
}

}//namespace intuitive

//}}} //namespace ithare::util::type


#endif //ithare_util_type_type_h_included
