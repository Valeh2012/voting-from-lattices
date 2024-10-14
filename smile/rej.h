#ifndef REJ_H
#define REJ_H

#include "poly.h"
#include "comm.h"

#define Rej0 PQMX_NAMESPACE(Rej0)
int Rej0(poly z[PQMX_ELL], poly v[PQMX_ELL]);

#define Rej1 PQMX_NAMESPACE(Rej1)
int Rej1(commrnd *z, commrnd* v);

#endif