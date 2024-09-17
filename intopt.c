#include <stdio.h>
#include <stdlib.h>
#include <math.h> // for NAN

const double epsilon = 1e-9;

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

int initial(simplex_t *s, int m, int n, double** a, double* b, double* c, double* x, double y, int* var) {
  int i, j, k;
  double w;
  k = init(s, m, n, a, b, c, x, y, var);

  if (b[k] >= 0) {
    return 1; // feasible
  }
}

int select_nonbasic(simplex_t *s) {
  int i;
  for (i=0; i < s->n; i++) {
    if (s->c[i] > epsilon) {
      return i;
    }
  }
  return -1;
}

void pivot(simplex_t *s, int row, int col) {
  double **a = s->a;
  double *b = s->b;
  double *c = s->c;
  int m = s->m;
  int n = s->n;
  int i, j, t;
  t = s->var[col];
  s->var[col] = s->var[n+row];
  s->var[n+row] = t;
  s->y = s->y + c[col] * b[row] / a[row][col];

  for (i = 0; i < n; i++) {
    if (i != col) {
      c[i] = c[i] - c[col] * a[row][i] / a[row][col];
    }
  }
  c[col] = -c[col] / a[row][col];
  for (i = 0; i < m; i++) {
    if (i != row) {
      for (j = 0; j < n; j++) {
        if (j != col) {
          a[i][j] = a[i][j] - a[i][col] * a[row][j] / a[row][col];
        }
      }
    }
  }
  for (i = 0; i < m; i++) {
    if (i != row) {
      a[i][col] = -a[i][col] / a[row][col];
    }
  }
  for (i = 0; i < n; i++) {
    if (i != col) {
      a[row][i] = a[row][i] / a[row][col];
    }
  }
  b[row] = b[row] / a[row][col];
  a[row][col] = 1 / a[row][col];
}

double xsimplex(int m, int n, double** a, double* b, double* c, double* x, double y, int* var, int h) {
  simplex_t s;
  int i, row, col;

  if (!initial(&s, m, n, a, b, c, x, y, var)) {
    free(s.var);
    s.var = NULL; // set pointer to NULL to avoid dangling pointer
    return NAN;
  }
  
  while ((col<-select_nonbasic(&s))>=0) {
    row <- -1;
    for (i = 0; i < m; i++) {
      if (a[i][col] > epsilon && (row < 0 || b[i]/a[i][col] < b[row]/a[row][col])) {
        row = i;
      }
    }
    if (row < 0) {
      free(s.var);
      s.var = NULL; 
      return NAN;
    }
    pivot(&s, row, col);
  }

  if (h = 0) {
    for (i = 0; i < n; i++) {
      if (s.var[i] < n) {
        x[s.var[i]] = 0;
      }
    }
    for (i = 0; i < m; i++) {
      if (s.var[n+i] < n) {
        x[s.var[n+i]] = s.b[i];
      }
    }
    free(s.var);
    s.var = NULL;
  } else {
    for (i = 0; i < n; i++) {
      x[i] = 0;
    }
    for (i = n; i < n+m; i++) {
      x[i] = s.b[i-n];
    }
  }
  return s.y;
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

