#include <string.h>
#include "bgv.h"
#include "randombytes.h"
#include "symmetric.h"
#include "math.h"
#include "reduce.h"

/*************************************************
* Name:        bgv_genkey
*
* Description: Generate a BGV public and private keypair
*
* Arguments:   - poly pk[2]: public key polynomials
*              - poly sk: secret key polynomial
**************************************************/
void bgv_genkey(poly pk[2], poly * sk)
{
    uint8_t seed[PQMX_SYMBYTES];
    unsigned int nonce=0;

    randombytes(seed, PQMX_SYMBYTES);
    hash_g(seed, seed, PQMX_SYMBYTES);  // Because kyber does so

    poly p = { 0 };
    p.coeffs[0] = PQMX_P;
    poly_ntt(&p);

    poly_uniform(&pk[0], seed, nonce); 
    poly_ntt(&pk[0]);                  // pk[0] = a = uniform poly

    poly e;
    poly_nonuniform(&e, 0x80, seed, nonce++);
    poly_ntt(&e);
    poly_basemul_montgomery(&e, &e, &p);  
    poly_tomont(&e);                         // e = pe

    poly_nonuniform(sk, 0x80, seed, nonce++);
    poly_ntt(sk);

    poly_basemul_montgomery(&pk[1], &pk[0], sk);
    poly_tomont(&pk[1]);
    poly_add(&pk[1], &pk[1], &e);
    poly_reduce(&pk[1]);                     // pk[1]=b=as+pe
    
}

/*************************************************
* Name:        bgv_enc
*
* Description: bgv Encryption of a polynomial
*
* Arguments:   - bgvcp *cp: pointer to output ciphertext
*              - const poly *msg: pointer to message polynomial
*              - const poly pk[2]: public key polynomials
*              - bgvrnd *rnd: pointer to output bgv randomnesses
**************************************************/
void bgv_enc(bgvcp* cp, const poly *msg, const poly pk[2], bgvrnd *rnd)
{
    uint8_t seed[PQMX_SYMBYTES];
    unsigned int nonce=0;

    randombytes(seed, PQMX_SYMBYTES);

    poly p = { 0 };
    p.coeffs[0] = PQMX_P;
    poly_ntt(&p);

    poly r,e1,e2;
    poly_nonuniform(&r, 0x80, seed, nonce++);
    poly_nonuniform(&e1, 0x80, seed, nonce++);
    poly_nonuniform(&e2, 0x80, seed, nonce++);

    rnd->r1 = r;
    rnd->r2 = e1;
    rnd->r3 = e2;

    poly_ntt(&r);
    poly_ntt(&e1);
    poly_ntt(&e2);

    poly_basemul_montgomery(&e1, &e1, &p);  
    poly_tomont(&e1);

    poly_basemul_montgomery(&e2, &e2, &p);  
    poly_tomont(&e2);  

    poly_basemul_montgomery(&cp->u, &pk[0], &r);
    poly_tomont(&cp->u);
    poly_add(&cp->u, &cp->u, &e1);
    poly_reduce(&cp->u);
    
    poly m = *msg;
    poly_ntt(&m);

    poly_basemul_montgomery(&cp->v, &pk[1], &r);
    poly_tomont(&cp->v);
    poly_add(&cp->v, &cp->v, &e2);
    poly_reduce(&cp->v);
    poly_add(&cp->v, &cp->v, &m);
    poly_reduce(&cp->v);
}

/*************************************************
* Name:        bgv_dec
*
* Description: bgv Decryption of a polynomial
*
* Arguments:   - poly *msg: pointer to output message polynomial
*              - const *poly sk: pointer to private key polynomial
*              - bgvcp *cp: pointer to input ciphertext
**************************************************/
void bgv_dec(poly *msg, const poly * sk, const bgvcp* cp)
{
    unsigned int j;
    poly res; 
    poly_basemul_montgomery(&res, sk, &cp->u);
    poly_tomont(&res);
    poly_sub(&res, &cp->v, &res);
    poly_reduce(&res);
    poly_invntt_tomont(&res);
    poly_reduce_mont(&res);
    poly_reduce(&res);

    int64_t t;
    for(j=0;j<PQMX_N;j++){
        t = res.coeffs[j] % PQMX_P; 
		t =  ( t >= 0 ? t : (t+PQMX_P));
        msg->coeffs[j] = t;
    }
    
}