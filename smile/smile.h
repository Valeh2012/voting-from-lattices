#ifndef SMILE_H
#define SMILE_H

#include <stdint.h>
#include "params.h"
#include "bgv.h"
#include "comm.h"

typedef struct {
  uint8_t hash[PQMX_SYMBYTES];
  uint8_t hash2[PQMX_SYMBYTES];
  poly z1[PQMX_ELL];
  commrnd z2;
  comm t;
  poly h;
} smile_proof;

void smile_prove(smile_proof *p, const uint8_t rho[2*PQMX_SYMBYTES], const poly bgvpk[2], poly candidates[PQMX_L], bgvrnd *rnd, const bgvcp *cipher, int iota);
int smile_proof_verify(smile_proof *p, const uint8_t rho[2*PQMX_SYMBYTES], const poly bgvpk[2], poly candidates[PQMX_L], const bgvcp *cipher);


#endif