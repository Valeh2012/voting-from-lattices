#include <string.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <time.h>
#include "randombytes.h"
#include "poly.h"
#include "bgv.h"

#define NTESTS 1000

int main(void)
{
  int i;
  uint8_t seed[PQMX_SYMBYTES];

  randombytes(seed, PQMX_SYMBYTES);

  poly bgvpk[2], bgvsk;
  bgv_genkey(bgvpk, &bgvsk);
  bgvrnd rnd;
  bgvcp cp;
  
  uint8_t msg[PQMX_INDCPA_MSGBYTES] = "Hello world!";
  printf("Plaintext: %s\n", (const char *)msg);

  poly m;
  poly_frommsg(&m, msg);

  clock_t tic = clock();
  for(i=0;i<NTESTS;i++)
    bgv_enc(&cp, &m, bgvpk, &rnd);
  clock_t toc = clock();
  printf("bgv encryption time in %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC / NTESTS);
  
  // verify
  tic = clock();
  for(i=0;i<NTESTS;i++)
    bgv_dec(&m, &bgvsk, &cp);
  toc = clock();
  printf("bgv decryption time in %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC / NTESTS);

  uint8_t msg2[PQMX_INDCPA_MSGBYTES] = {0};

  poly_tomsg(msg2, &m);
  printf("Decrypted ciphertext: %s\n",(const char *)msg2);

  poly_frommsg(&m, msg);
  poly m2;

  int failures = 0;
  int l1_norms[NTESTS];
  for(i=0;i<NTESTS;i++){
    bgv_enc(&cp, &m, bgvpk, &rnd);
    bgv_dec(&m2, &bgvsk, &cp);
    poly_tomsg(msg2, &m2);
    failures += memcmp(msg, msg2, PQMX_INDCPA_MSGBYTES) != 0;
  }
  printf("number of failures: %d\n", failures);

  return 0;
}
