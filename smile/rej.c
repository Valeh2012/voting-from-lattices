#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include "rej.h"
#include "randombytes.h"

int Rej0(poly z[PQMX_ELL], poly v[PQMX_ELL]){
    long seed;
    randombytes((uint8_t*) &seed, sizeof(long));
    gmp_randstate_t rndstate;
    gmp_randinit_mt(rndstate);
    gmp_randseed_ui(rndstate, seed);
    mpfr_set_default_prec(200);

    mpfr_t u;
    mpfr_init(u);
    mpfr_urandomb(u, rndstate);

    double M = sqrt(1.5);
    mpfr_t m;
    mpfr_init(m);
    mpfr_set_d(m, M, MPFR_RNDD);

    mpfr_t zv, v2, tmp, res, s;
    mpfr_inits(zv, v2, tmp, res, s,(mpfr_ptr) 0);
    mpfr_set_d(res, 0.0, MPFR_RNDD);
    mpfr_set_d(v2, 0.0, MPFR_RNDD);
    mpfr_set_d(zv, 0.0, MPFR_RNDD);

    for(int i=0; i< PQMX_ELL; i++){
        for(int j=0; j<PQMX_N;j++){
            mpfr_set_si(tmp, v[i].coeffs[j], MPFR_RNDD);
            mpfr_sqr(tmp, tmp, MPFR_RNDD);
            mpfr_add(v2, v2, tmp, MPFR_RNDD);  

            mpfr_set_si(tmp, z[i].coeffs[j], MPFR_RNDD);
            mpfr_mul_si(tmp, tmp, v[i].coeffs[j], MPFR_RNDD);

            mpfr_add(zv, zv, tmp, MPFR_RNDD);
        }
    } 


    mpfr_sub(res, v2, zv, MPFR_RNDD);
    mpfr_sub(res, res, zv, MPFR_RNDD);    // res= -<z,v>+|v|^2

    mpfr_set_d(s, SIGMA0, MPFR_RNDD);
    mpfr_sqr(s, s, MPFR_RNDD);
    mpfr_mul_si(s, s, 2, MPFR_RNDD);  // 2s^2
    mpfr_div(res, res, s, MPFR_RNDD); // res = (-<z,v>+|v|^2) / 2s^2
    mpfr_exp(res, res, MPFR_RNDD);    // res = exp(res)
    mpfr_div(res, res, m, MPFR_RNDD);  // res = res / m

    int ans = mpfr_cmp(u, res);
    

    mpfr_clear(u);
    mpfr_clear(m);
    mpfr_clear(tmp);
    mpfr_clear(res);
    mpfr_clear(zv);
    mpfr_clear(v2);
    mpfr_clear(s);
    gmp_randclear(rndstate);
    return (ans == 1 ? 1 : 0);
}
int Rej1(commrnd *z, commrnd* v){
    long seed;
    randombytes((uint8_t*) &seed, sizeof(long));
    gmp_randstate_t rndstate;
    gmp_randinit_mt(rndstate);
    gmp_randseed_ui(rndstate, seed);
    mpfr_set_default_prec(200);

    mpfr_t u;
    mpfr_init(u);
    mpfr_urandomb(u, rndstate);

    double M = sqrt(1.5);
    mpfr_t m;
    mpfr_init(m);
    mpfr_set_d(m, M, MPFR_RNDD);

    mpfr_t zv, v2, tmp, res, s;
    mpfr_inits(zv, v2, tmp, res, s, (mpfr_ptr) 0);
    mpfr_set_d(res, 0.0, MPFR_RNDD);
    mpfr_set_d(v2, 0.0, MPFR_RNDD);
    mpfr_set_d(zv, 0.0, MPFR_RNDD);

    poly zz[PQMX_KAPPA], vv[PQMX_KAPPA];
    for(int i=0; i<PQMX_LAMBDA; i++){
        zz[i] = z->s[i];
        vv[i] = v->s[i];
    }
    for(int i=0; i<PQMX_MU; i++){
        zz[PQMX_LAMBDA + i] = z->e[i];
        vv[PQMX_LAMBDA + i] = v->e[i];
    }
    for(int i=0; i<PQMX_M; i++){
        zz[PQMX_LAMBDA + PQMX_MU + i] = z->em[i];
        vv[PQMX_LAMBDA + PQMX_MU + i] = v->em[i];
    }
    polyvec_sub(vv, zz, vv, PQMX_KAPPA);
    polyvec_reduce(vv, PQMX_KAPPA);

    for(int i=0; i< PQMX_KAPPA; i++){
        for(int j=0; j<PQMX_N;j++){
            mpfr_set_si(tmp, vv[i].coeffs[j], MPFR_RNDD);
            mpfr_sqr(tmp, tmp, MPFR_RNDD);
            mpfr_add(v2, v2, tmp, MPFR_RNDD);  

            mpfr_set_si(tmp, zz[i].coeffs[j], MPFR_RNDD);
            mpfr_mul_si(tmp, tmp, vv[i].coeffs[j], MPFR_RNDD);

            mpfr_add(zv, zv, tmp, MPFR_RNDD);
        }
    } 

    int ans = mpfr_cmp_si(zv, 0);
    if(ans < 0){
        mpfr_clear(u);
        mpfr_clear(m);
        mpfr_clear(tmp);
        mpfr_clear(res);
        mpfr_clear(zv);
        mpfr_clear(v2);
        mpfr_clear(s);
        gmp_randclear(rndstate);
        return 1;
    }

    mpfr_sub(res, v2, zv, MPFR_RNDD);
    mpfr_sub(res, res, zv, MPFR_RNDD); // res= -<z,v>+|v|^2

    mpfr_set_d(s, SIGMA1, MPFR_RNDD);
    mpfr_sqr(s, s, MPFR_RNDD);
    mpfr_mul_si(s, s, 2, MPFR_RNDD);  // 2s^2
    mpfr_div(res, res, s, MPFR_RNDD); // res = (-<z,v>+|v|^2) / 2s^2
    mpfr_exp(res, res, MPFR_RNDD);    // res = exp(res)
    mpfr_div(res, res, m, MPFR_RNDD);  // res = res / m

    ans = mpfr_cmp(u, res);

    mpfr_clear(u);
    mpfr_clear(m);
    mpfr_clear(tmp);
    mpfr_clear(res);
    mpfr_clear(zv);
    mpfr_clear(v2);
    mpfr_clear(s);
    gmp_randclear(rndstate);
    return (ans == 1 ? 1 : 0);
    return 0;
}