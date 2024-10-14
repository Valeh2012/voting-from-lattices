/*
 * params.hpp
 *
 *  Created on: Sep 5, 2022
 *      Author: valeh
 */

#ifndef PARAMS_HPP_
#define PARAMS_HPP_

#include <gmpxx.h>
#include <nfl.hpp>

namespace parameters{
	using poly_t = nfl::poly_from_modulus<uint64_t, 4096, nfl::params<uint64_t>::kModulusBitsize>;
	const size_t MU=1, LAMBDA=1;
}

using poly_t = nfl::poly_from_modulus<uint64_t, 4096, nfl::params<uint64_t>::kModulusBitsize>;

#define plaintextModulus 10000009  // set p next_prime after NV*LC
#define NV 100                    // number of voters
#define NC 100                     // number of candidates
#define LC 32                      // maximum vote that can be given to a single candidate
#define LS 100                     // total amount of votes a single voter can give - vote budget
#define delta_1 (1L<<32)           // rejection sampling parameter
#define beta_1 1024                // rejection sampling and MSIS parameter
#define k 2                        // soundness boost and automorphism order parameter

#define SYMBYTES 32                // output of hash functions in bytes

#endif /* PARAMS_HPP_ */
