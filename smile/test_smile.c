#include <string.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <time.h>
#include "randombytes.h"
#include "poly.h"
#include "comm.h"
#include "util.h"
#include "bgv.h"
#include "smile.h"


int main(void)
{

  uint8_t seed[2*PQMX_SYMBYTES];

  randombytes(seed, 2*PQMX_SYMBYTES);

  poly bgvpk[2], bgvsk;
  bgv_genkey(bgvpk, &bgvsk);


  bgvcp cp;   //ballot box
  bgvrnd bgv_rnd; //bgv randomnesses

  // Randomly generate a ballot and encrpyt
  poly *candidates = (poly *) aligned_alloc(32, PQMX_L*sizeof(poly));
  poly m;
  uint8_t paddedballot[PQMX_INDCPA_MSGBYTES];
  char* choice = (char *) malloc(32);
  for(int i=0; i< PQMX_L; i++){
    memset(choice, 0, 32);
    sprintf(choice,"cand_%d", i);
    padding(paddedballot, (uint8_t*) choice, strlen(choice));
    poly_frommsg(&candidates[i], paddedballot);
  }
  
  int i = secure_random(PQMX_L);
  bgv_enc(&cp, &candidates[i], bgvpk, &bgv_rnd);
  printf("Encrypted ballot: cand_%d\n", i);
  
//   polyvec_invntt_tomont(bgvpk,2);
//   polyvec_reduce_mont(bgvpk,2);

  // generate zero-knowledge proof 
  smile_proof p;
  clock_t tic = clock();
  smile_prove(&p, seed, bgvpk, candidates, &bgv_rnd, &cp, i);
  clock_t toc = clock();
  printf("proof generated in %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);

  // verify the proof
  tic = clock();
  int res = smile_proof_verify(&p, seed, bgvpk, candidates, &cp);
  toc = clock();
  printf("proof %s in %f seconds\n", (res ? "rejected" : "verified"), (double)(toc - tic) / CLOCKS_PER_SEC);

  // Start tallying
  uint8_t ballot[32], ballotOld[32];
  bgv_dec(&m, &bgvsk, &cp);
  poly_tomsg(paddedballot, &m);
  remove_padding(ballot, paddedballot);
  poly_tomsg(paddedballot, &candidates[i]);
  remove_padding(ballotOld, paddedballot);
  printf("Decrypted ballot: ");
  puts((const char*)ballot);


  printf("\n");
  if(strcmp((const char*)ballot, (const char*)ballotOld)){
    printf("Decrypted ciphertext does not match actual plaintext!\n");
  }

  printf("proof size: %ld\n", sizeof(smile_proof));


  free(p.t.t0);
  free(p.t.tm);
  free(p.z2.e);
  free(p.z2.em);
  free(p.z2.s);
  
  return 0;
}

