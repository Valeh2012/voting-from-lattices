#ifndef BGV_H
#define BGV_H

#include <stdlib.h>
#include "params.h"
#include "poly.h"

typedef struct{
  poly r1;
  poly r2;
  poly r3;
}bgvrnd;

typedef struct{
    poly u;
    poly v;
}bgvcp;

#define bgv_genkey PQMX_NAMESPACE(bgv_keygen)
void bgv_genkey(poly pk[2], poly * sk);

#define bgv_enc PQMX_NAMESPACE(bgv_enc)
void bgv_enc(bgvcp* cp, const poly *msg, const poly pk[2], bgvrnd *rnd);

#define bgv_dec PQMX_NAMESPACE(bgv_dec)
void bgv_dec(poly *msg, const poly * sk, const bgvcp* cp);

#endif