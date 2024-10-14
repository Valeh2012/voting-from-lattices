#ifndef UTIL_H
#define UTIL_H

#include <stdint.h>
#include <stddef.h>
#include "params.h"
#include "poly.h"
#include "polyvec.h"

#define maxab(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define poly_norm PQMX_NAMESPACE(poly_norm)
double poly_norm(const poly *a, int l);

#define polyvec_norm PQMX_NAMESPACE(polyvec_norm)
double polyvec_norm(const poly *a, int l, int len);

#ifdef GAUSSIAN
#define poly_gaussian PQMX_NAMESPACE(poly_gaussian)
void poly_gaussian(poly *r, const double sigma);
#endif
int secure_random(int j);
void shuffle_array(int64_t *arr, int len);

void padding(uint8_t out[PQMX_INDCPA_MSGBYTES], const uint8_t *buf, size_t buflen);
void remove_padding(uint8_t *out, const uint8_t in[PQMX_INDCPA_MSGBYTES]);

#endif