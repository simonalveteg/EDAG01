#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

#include "error.h"
#include "poly.h"


typedef struct poly_t {
  /*int *exponents;*/ // better but more annoying.
  int *coefficients;
  int degree; // to know how large the arrays are.
} poly_t;

poly_t* empty_poly_of_degree(int degree) {
  poly_t* poly = calloc(1, sizeof(poly_t));
  poly->coefficients = calloc(degree+1, sizeof(int));
  poly->degree = degree;

  return poly;
}

int next_number(const char* str, int* start, int length) {
  int i = *start;
  int cint;
  int num = 0; 
  while (isdigit(str[i])) {
    cint = str[i] - '0';
    if (num != 0) 
      cint += 10*num; // Add previous number.

    num = cint;
    i++;
  }
  *start = i; // Update index in calling function.
  return num;
}

poly_t* new_poly_from_string(const char* str) {
  
}

void free_poly(poly_t* poly) {
  free(poly->coefficients);
  free(poly);
}

poly_t*	mul(poly_t* a, poly_t* b) {
  int i, j, c, e;
  
  poly_t *r = empty_poly_of_degree(a->degree + b->degree);

  for (i = 0; i < a->degree; i++) {
    for (j = 0; j < b->degree; j++) {
      c = a->coefficients[i] * b->coefficients[j];
      e = i + j; // degree given by index.

      // check if coefficient is already present.
      if (r->coefficients[e] != 0)
        c += r->coefficients[e];

      // update coefficient at index of exponent.
      r->coefficients[e] = c; 
    }
  }
  return r;
}

void print_poly(poly_t* p) {
  int i, coef;
  char sign;
  for (i = p->degree; i > 0; i--) {
    coef = p->coefficients[i];
    sign = (coef > 0) ? '+' : '-';
    if (coef != 0) {
      if (i == p->degree && sign == '+')
        printf("%d", coef);
      else 
        printf("%c %d", sign, coef);
    }
    if (i > 1)
      printf("x^%d ", i);
    else if (i == 1)
      printf("x ");
  }
  printf("\n");  
}
