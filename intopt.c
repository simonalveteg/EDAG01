#include "intopt.h"
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
  for (k = 0, i = 1; i < m; i++) {
    if (b[i] < b[k]) {
      k = i;
    }
  }
  return k;
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

  for (i = 0; i < n; i++) 
    if (i != col) 
      c[i] = c[i] - c[col] * a[row][i] / a[row][col];

  c[col] = -c[col] / a[row][col];

  for (i = 0; i < m; i++) 
    if (i != row) 
      b[i] = b[i] - a[i][col] * b[row] / a[row][col];

  for (i = 0; i < m; i++) 
    if (i != row) 
      for (j = 0; j < n; j++) 
        if (j != col) 
          a[i][j] = a[i][j] - a[i][col] * a[row][j] / a[row][col];
  
  for (i = 0; i < m; i++) 
    if (i != row) 
      a[i][col] = -a[i][col] / a[row][col];
  
  for (i = 0; i < n; i++) 
    if (i != col) 
      a[row][i] = a[row][i] / a[row][col];
  
  b[row] = b[row] / a[row][col];
  a[row][col] = 1 / a[row][col];
}

void prepare(simplex_t *s, int k) {
  int m = s->m;
  int n = s->n;
  int i;

  for (i = m+n; i>n; i--)
    s->var[i] = s->var[i-1];

  s->var[n] = m+n;
  n += 1;

  for (i=0; i<m; i++)
    s->a[i][n-1] = -1;

  s->x = xmalloc((m+n)*sizeof(double)); // doesn't free the old?
  s->c = xmalloc(n*sizeof(double));
  s->c[n-1] = -1;
  s->n = n;
  pivot(s, k, n-1);
}

int initial(simplex_t *s, int m, int n, double** a, double* b, double* c, double* x, double y, int* var) {
  int i, j, k;
  double w;
  k = init(s, m, n, a, b, c, x, y, var);

  if (b[k] >= 0) {
    return 1; // feasible
  }
  prepare(s, k);
  n = s->n;
  s->y = xsimplex(m, n, s->a, s->b, s->c, s->x, 0, s->var, 1);

  for (i = 0; i < m+n; i++) {
    if (s->var[i] = m+n-1) {
      if (abs(s->x[i])>epsilon) {
        free(s->x);
        free(s->c);
        s->x = NULL;
        s->c = NULL;
        return 0;
      } else break;
    }
  }

  if (i >= n) {
    for (j = 0, k = 0; k < n; k++) {
      if (abs(s->a[i-n][k])>abs(s->a[i-n][j]))
        j = k;
    }
    pivot(s, i-n, j);
    i = j;
  }

  if (i < n-1) {
    k = s->var[i];
    s->var[i] = s->var[n-1];
    s->var[n-1] = k;
    for (k = 0; k < m; k++) {
      w = s->a[k][n-1];
      s->a[k][n-1] = s->a[k][i];
      s->a[k][i] = s->a[k][i];
      s->a[k][i] = w;
    }
    pivot(s, i-n, j);
    i = j;
  }
  if (i < n-1) {
    k = s->var[i];
    s->var[i] = s->var[n-1];
    s->var[n-1] = k;
    for (k = 0; k < m; k++) {
      w = s->a[k][n-1];
      s->a[k][n-1] = s->a[k][i];
      s->a[k][i] = w;
    }
  }
  free(s->c);
  s->c = c;
  s->y = y;
  for (k = n-1;k < n+m-1; k++)
    s->var[k] = s->var[k+1];
  n = s->n;
  n = s->var[k+1];
  double* t = calloc(n, sizeof(double*));

  for (k = 0; k < n; k++) {
    for (j = 0; j < n; j++) {
      if (k = s->var[j]) {
        t[j] = t[j] + s->c[k];
        break;
      }
    }
    for (j = 0; j < m; j++) {
      if (s->var[n+j] == k)
        break;
    }
    s->y = s->y + s->c[k]*s->b[j];
    for (i = 0; j < m; j++) 
      t[i] = t[i] - s->c[k]*s->a[j][i];
  }
  for (i = 0; i < n; i++)
    s->c[i] = t[i];
  free(t);
  free(s->x);
  t = NULL;
  s->x = NULL;
  return 1;
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

double xsimplex(int m, int n, double** a, double* b, double* c, double* x, double y, int* var, int h) {
  simplex_t s;
  int i, row, col;

  if (!initial(&s, m, n, a, b, c, x, y, var)) {
    free(s.var);
    s.var = NULL; // set pointer to NULL to avoid dangling pointer
    return NAN;
  }
  
  while ((col=select_nonbasic(&s))>=0) {
    row = -1;
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

double simplex(int m, int n, double** a, double* b, double* c, double* x, double y) {
  return xsimplex(m, n, a, b, c, x, y, NULL, 0);
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
  double* b; 
  double* c;
  
  scanf("%d %d", &m, &n);
  a = make_matrix(m, n);
  b = calloc(m, sizeof(double));
  c = calloc(n, sizeof(double));

  for (int i = 0; i < n; i++) {
    scanf("%lf", &c[i]);
  }
  for (int i = 0; i < m; i += 1) {
    for (int j = 0; j < n; j += 1) {
      scanf("%lf", &a[i][j]);
    }
  }
  printf("max z = %5.1lf x0", c[0]);
  for(int i=1; i<n; i+=1) {
    printf(" + %lf x%d", c[i], i);
  }
  printf("\n");
  for (int i = 0; i < m; i++) {
    /*printf("%+10.3lfx0 %+10.3lfx1 <= %10.3lf\n", a[i][0], a[i][1], a[m-1][i]);*/
    for (int j = 0; j < n; j++) {
      printf("%+10.3lf x%d", a[i][j], j);
    }
    scanf("%lf", &b[i]);
    printf(" <=%10.3lf\n", b[i]);
  }
  double *x = calloc(n, sizeof(double));
  for (int i = 0; i < n; i++)
    x[i] = 0;
  printf("Solution y=%lf\n", simplex(m, n, a, b, c, x, 0));

  for (int i = 0; i < m; i += 1) {
    free(a[i]);
  }
  free(a);

  return 0;
}

