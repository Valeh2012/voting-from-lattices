#pragma once

#include <cstddef>

#include <algorithm>
#include <iostream>

extern "C" {
#include "fips202.h"
}


namespace utils{

using P = parameters::poly_t;
using V = typename P::value_type;

P sample_ternary(uint8_t mode, uint8_t* hash){
	P p;
	const V pm = P::get_modulus(0) - 1u;
	uint8_t rnd[P::degree];
	V tmp[P::degree];
	if(hash == nullptr) nfl::fastrandombytes(rnd, P::degree); // TODO: remove randomness for null pointer
	else shake256(rnd,sizeof(rnd),hash,32);  // TODO: remove XOF after hashing
	for(std::size_t i=0;i<P::degree;i++){
		tmp[i] = rnd[i] <= mode ? pm + (rnd[i] & 2)  : 0u;
	}
	p.set(tmp, tmp+P::degree, true);
	return p;
}

void sample_uniform(P *out, uint8_t* hash, int j){

	size_t outlen = j*sizeof(P);
	uint8_t *rnd = (uint8_t*) malloc(outlen);
	//	uint8_t rnd[(j+1)*param_N*2];
	memset(rnd, 0, outlen); //TODO: change to incremental api.
	shake256(rnd,outlen,hash, 32);
	for(int i=0;i<j;i++){
		V * ptr = (V*) (rnd+i*sizeof(P));
		out[i].set(ptr, ptr+ P::degree, true);
	}
	free(rnd);
}

void sample_uniform_scalar_ntt(P *out, uint8_t* hash,int j){
	size_t outlen = j*sizeof(V);
	uint8_t *rnd = (uint8_t*) malloc(outlen);
	//	uint8_t rnd[(j+1)*param_N*2];
	memset(rnd, 0, outlen); //TODO: change to incremental api.
	shake256(rnd,outlen,hash, 32);
	for(int i=0;i<j;i++){
		V* ptr = (V*) (rnd + i*sizeof(V));
		fill(out[i].begin(), out[i].end(), *ptr);
	}
	free(rnd);
}


void poly_transpose(P &in){
	V tmp;
	V q = P::get_modulus(0);
	auto l = P::degree;
	in.data()[l/2] = q - in.data()[l/2];
	for(size_t i=1;i< l / 2;i++){
		tmp = in.data()[i];
		in.data()[i] = q - in.data()[l-i];
		in.data()[l-i] = q - tmp;
	}
}


void poly_inner_product_Zq(mpz_t &res, P &a, P &b){
	mpz_t tmp1, tmp2;
	mpz_init(tmp1);
	mpz_init(tmp2);
	for(size_t i=0;i<P::degree;i++){
		mpz_set_ui(tmp1, a.data()[i]);
		mpz_set_ui(tmp2, b.data()[i]);
		mpz_addmul(res, tmp1, tmp2);
	}
}



void galois_transform(P &out, P &p, size_t i){
	// re-think about implementation

	// check gcd(i,2*degree) == 1

	size_t l  = P::degree;

	if( gcd(l, i) != 1 ){
		cerr << "no such transformation" << endl;
		exit(-1);
		return;
	}

	V modulus = P::get_modulus(0);
	P res{0}, tmp;
	int q,r;
	for(size_t j=0;j<l;j++){
		if(p.data()[j] == 0) continue;
		tmp.set(0);
		q = j*i / l;
		r = (j*i)%l;
		if(q & 1) tmp.data()[r] = modulus - p.data()[j];
		else tmp.data()[r] = p.data()[j];
		res = res + tmp;
	}


	out = res;
}

template<typename _Mn>
_Mn modInverse(_Mn a, _Mn m){
	_Mn m0 = m;
	_Mn y = 0, x = 1;

	if (m == 1)
		return 0;

	while (a > 1)
	{
		// q is quotient
		_Mn q = a / m;
		_Mn t = m;

		// m is remainder now, process same as
		// Euclid's algo
		m = a % m, a = t;
		t = y;

		// Update y and x
		y = x - q * y;
		x = t;
	}

	// Make x positive
	if (x < 0)
		x += m0;

	return x;
}

void inv_galois_transform(P& out, P &p, int i){
	// re-think about implementation

	int l = P::degree;

	if( !gcd(l, i) ){
		cerr << "no such transformation" << endl;
		exit(-1);
		return;
	}

	//	cout << i << endl;
	// find multiplicative inverse mod 2l
	i = modInverse(i, 2*l);
	//	cout << i << endl;
	//	poly_type res;
	p.invntt_pow_invphi();
	galois_transform(out,p,i);
	out.ntt_pow_phi();
	//	memcpy(&out, &res, sizeof(out));
}


V _abs(V a){
	P::signed_value_type tmp1;
	BGV::util::center(tmp1, a, P::get_modulus(0), P::get_modulus(0) >> 1);
	return (V) abs(tmp1);
}

template<std::size_t LEN>
V inf_norm(std::array<P, LEN> &polyvec){
	V max_el = 0;
	V tmp = 0;
	for(auto &it : polyvec){
		for(auto &el : it){
			tmp =  _abs(el);
			max_el = max(max_el, tmp);
		}
	}

	return max_el;
}

P coeffs_sum(P const& poly){
	V sum = 0;
	for(auto &c : poly){
		sum += c;
		sum = sum;
	}
	return P{sum};
}

}

