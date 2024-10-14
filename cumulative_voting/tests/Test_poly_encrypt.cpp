/**
 *
 * This file is part of BGV-NFLlib
 *
 * Copyright (C) 2016  CryptoExperts
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
 *
 */

#include <cstddef>

#include <gmpxx.h>
#include <chrono>
#include <iostream>
#include <nfl.hpp>
#include <thread>
#include <vector>
#include <chrono>

/// include the BGV homomorphic encryption library
namespace BGV {
namespace params {
using poly_t = nfl::poly_from_modulus<uint32_t, 128, nfl::params<uint32_t>::kModulusBitsize>;
poly_t::value_type plaintextModulus = 90001;
long plaintextModulusBitsize = 17;
}
}  // namespace BGV::params
#include <BGV.hpp>


#define NTESTS 10000

int main() {
  // Seed (for deterministic values)
  srand(0);

  using P = BGV::params::poly_t;
  std::chrono::steady_clock::time_point begin,end;

  P::value_type phi = nfl::params<P::value_type>::primitive_roots[0];
  for (unsigned int i = 0 ;
	  i < nfl::static_log2<nfl::params<P::value_type>::kMaxPolyDegree>::value - nfl::static_log2<P::degree>::value ; i++)
  {
	phi = nfl::ops::mulmod<P::value_type, nfl::simd::serial>{}(phi, phi, 0);
  }
  std::cout << phi << std::endl;

  // Keygen
  BGV::sk_t secret_key;
  BGV::pk_t public_key(secret_key);

  // Polynomials
  P m[2], t{nfl::uniform()};

  m[0] = {1, 2, 3, 4, 5, 6, 7, 8};
  m[1] = {1,1};
  m[1].ntt_pow_phi();
  std::cout << m[1] << std::endl;

  std::cout << P::get_modulus(0) << std::endl;
  BGV::util::reduce(t, BGV::params::plaintextModulus, P::get_modulus(0), public_key.qDivBy2);
//  std::cout << t << std::endl;
  t.ntt_pow_phi();

  // Encrypt
  BGV::ciphertext_t c;
  auto total=0;
  begin = std::chrono::steady_clock::now();
  for(int i=0;i<NTESTS;i++){
	  BGV::encrypt_poly(c, public_key, t);
  }
  end = std::chrono::steady_clock::now();
  total = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
  std::cout << "Avg enc time: " << total / NTESTS << "[μs]" <<std::endl;

  t.invntt_pow_invphi();
  // Initialize polym
  P polym0;


  // decrypt to the polym
  total=0;
  begin = std::chrono::steady_clock::now();
  for(int i=0;i<NTESTS;i++){
	  BGV::decrypt_poly(polym0, secret_key, public_key, c);
  }
  end = std::chrono::steady_clock::now();
  total = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
  std::cout << "Avg dec time: " << total / NTESTS << "[μs]" <<std::endl;

  bool f = true;
  for(size_t i=0; i< P::degree;i++){
	  f = f & (polym0.data()[i] == t.data()[i]);
	  if(!f) break;
  }
  std::cout << (f ? "Success" : "Failure" )<< std::endl;
//  std::cout << polym0 << std::endl;


  return 0;
}
