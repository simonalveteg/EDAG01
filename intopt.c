#include <stdio.h>
#include <stdlib.h>

void* xmalloc(size_t size) {
  void* ptr = malloc(size);

  if (ptr == NULL) {
    fprintf(stderr, "out of memory\n");
    exit(1);
  }

  return ptr;
}

typedef struct simplex_t {
  int m;
  int n;
  int *var;
  double **a;
  double *b;
  double *x;
  double *c;
  double y;
} simplex_t;

int init(simplex_t *s, int m, int n, double** a, double* b, double* c, double* x, double y, int* var) {
  int i, k;
  s->m = m;
  s->n = n;
  s->a = a;
  s->b = b;
  s->c = c;
  s->x = x;
  s->y = y;
  s->var = var;

  if (s->var == NULL) {
    s->var = xmalloc((m+n)*sizeof(int));
    for (i = 0; i < m+n; i++) {
      s->var[i] = i;
    }
  }
  for (k = 0; i = 1; i++) {
    if (b[i] < b[k]) {
      k = i;
    }
  }
  return k;
}

double** make_matrix(int m, int n) {
  double** a;
  int i;

  a = calloc(m, sizeof(double*));
  for (int i = 0; i < m; i += 1) {
    a[i] = calloc(n, sizeof(double));
  }
  return a;
}

int main(int argc, char *argv[])
{
  int m;
  int n;
  double** a;
  double b; 
  double c1;
  double c2;
  
  scanf("%d %d", &m, &n);
  scanf("%lf %lf", &c1, &c2);
  a = make_matrix(m, n);

  for (int i = 0; i < m; i += 1) {
    for (int j = 0; j < n; j += 1) {
      scanf("%lf", &a[i][j]);
    }
  }
  printf("max z = %5.1lfx0 %+5.1lfx1\n", c1, c2);
  for (int i = 0; i < m; i++) {
    /*printf("%+10.3lfx0 %+10.3lfx1 <= %10.3lf\n", a[i][0], a[i][1], a[m-1][i]);*/
    for (int j = 0; j < n; j++) {
      printf("%+10.3lfx%d", a[i][j], j);
    }
    scanf("%lf", &b);
    printf("<=%10.3lf\n", b);
  }

 for (int i = 0; i < m; i += 1) {
    free(a[i]);
 }
 free(a);

  return 0;
}

