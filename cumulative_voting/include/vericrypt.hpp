/*
 * vericrypt.hpp
 *
 *  Created on: Dec 13, 2022
 *      Author: valeh
 */

#ifndef INCLUDE_VERICRYPT_HPP_
#define INCLUDE_VERICRYPT_HPP_

#include <iostream>
#include <string>
#include <memory>
#include <random>
#include <chrono>
#include <vector>
#include <numeric>

//using namespace std;

#include "params.hpp"
#include <BGV.hpp>
#include <BDLOP.hpp>

typedef BGV::ciphertext_t ctx_t;
typedef BGV::pk_t pk_t;
typedef BGV::sk_t sk_t;

using P = parameters::poly_t;

template<size_t KAPPA>
struct Proof{
	BDLOP::comm_t<KAPPA> t;
	array<BDLOP::comm_rnd_t<KAPPA>, k > z;
	uint8_t chash[SYMBYTES];
	ctx_t ct;
	P h;
};


#define kappa (LC+8)
using proof_t = Proof<kappa>;
typedef BDLOP::comm_pk_t<kappa> comm_key_t;

#define nmod(a, m) ((a) < 0 ? (m - (-1*(a))%m) : (a)%m)

void vericrypt(proof_t &p, const pk_t &pubkey,const comm_key_t& commkey, const P &message);

bool verify(proof_t &p, const comm_key_t &commkey, const pk_t &pubkey);

#endif /* INCLUDE_VERICRYPT_HPP_ */
