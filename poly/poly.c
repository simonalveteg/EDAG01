#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>

#include "error.h"
#include "poly.h"

typedef struct term_t {
    int coeff;
    int exp;
    struct term_t* next;
} term_t;

struct poly_t {
    term_t* first;
    term_t* last;
};

void free_poly(poly_t* poly) {
  if (poly == NULL) return;
  term_t *term = poly->first;
  term_t *next;

  while (term != NULL) {
    next = term->next;
    term->next = NULL;
    free(term);
    term = next;
  }
  poly->first = poly->last = NULL;
  free(poly);
}

void print_poly(poly_t* poly) {

  term_t *term = poly->first;
  int first = 1; // Flag to handle the first term.

  while (term != NULL) {
    char sign = (term->coeff > 0) ? '+' : '-';
    int coeff = abs(term->coeff);
    int exp = term->exp;
  
    // Print sign if not first element, since first is always positive.
    if (!first) {
      printf(" %c ", sign);
    } else {
      first = 0;
    }
    if (coeff > 1 || exp == 0)
      printf("%d", coeff);
    if (exp > 0)
      printf("x");
    if (exp > 1)
      printf("^%d", exp);

    term = term->next; // Move to the next term
  }
  printf("\n");
}

void add_term(poly_t* poly, int coeff, int exp) {
  // Create a new term and add it to the polynomial.
  term_t* new_term = malloc(sizeof(term_t));
  new_term->coeff = coeff;
  new_term->exp = exp;
  new_term->next = NULL;

  // Add first if no terms are present.
  if (poly->first == NULL && poly->last == NULL) {
    poly->first = new_term;
    poly->last = new_term;
    return;
  } 

  // update term with same exponent if possible.
  term_t *current = poly->first;
  while (current != NULL) {
    if (current->exp == exp) {
      current->coeff += coeff;
      free(new_term);
      return;
    }
    current = current->next;
  }
  
  // Add term last.
  poly->last->next = new_term;
  poly->last = new_term;
}

poly_t* new_poly_from_string(const char* str) {

  poly_t *poly = malloc(sizeof(poly_t));
  poly->first = poly->last = NULL;
  
  const char* s = str;
  char sign = '+';
  int num = 0;
  int c = 0;
  int e = 0;

  while (*s) {

    if (*s == '-' || *s == '+') {
      sign = *s;
      s++;
      continue;
    }

    if (isdigit(*s)) {
      num = 10*num + (*s - '0');
      s++;
      continue;
    }

    if (*s == 'x') {
      if (num == 0) {
        num = 1; // no coefficient, but x is present so it is 1.
      }
      c = (sign == '-') ? -num : num;
      e = 1; // at least one x present.
      num = 0;
      s++;
      continue;
    }

    if (*s == ' ') {
      // c is 0 => don't neew to add new term.
      if (c == 0) {
        s++;
        continue;
      }
      // update e if a number has been found.
      if (num != 0)
        e = num;
      
      // add term to poly.
      add_term(poly, c, e);
      // reset state.
      sign = '+';
      c = 0;
      e = 0;
      num = 0;
    }

    s++;
  }

  // end of string reached, add term if possible.
  if (num != 0) {
    add_term(poly, num, e);
  }
  
  return poly;
}

poly_t* mul(poly_t* a, poly_t* b) {
  term_t *at = a->first;
  term_t *bt = b->first;
  int c, e;

  poly_t *r = malloc(sizeof(poly_t));
  r->first = r->last = NULL;

  while (at != NULL) {
    while (bt != NULL) {
      c = at->coeff * bt->coeff;
      e = at->exp + bt->exp;
      add_term(r, c, e);
      bt = bt->next;
    }
    at = at->next;
    bt = b->first;
  }
  return r;
}

// TEST MAIN
/*static void poly_test(const char* a, const char* b)*/
/*{*/
/*	poly_t*		p;*/
/*	poly_t*		q;*/
/*	poly_t*		r;*/
/**/
/*	printf("Begin polynomial test of (%s) * (%s)\n", a, b);*/
/**/
/*	p = new_poly_from_string(a);*/
/*	q = new_poly_from_string(b);*/
/**/
/*	print_poly(p);*/
/*	print_poly(q);*/
/**/
/*	r = mul(p, q);*/
/**/
/*	print_poly(r);*/
/**/
/*	free_poly(p);*/
/*	free_poly(q);*/
/*	free_poly(r);*/
/**/
/*	printf("End polynomial test of (%s) * (%s)\n", a, b);*/
/*}*/
/**/
/*int main(void)*/
/*{*/
/*	poly_test("x^2 - 7x + 1", "3x + 2");*/
/*	putchar('\n');*/
/*	poly_test("x^10000000 + 2", "2x^2 + 3x + 4");*/
/**/
/*	return 0;*/
/*}*/
