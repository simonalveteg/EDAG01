#include <stdio.h>
#include <stdlib.h>
#include <math.h> // for NAN

// header

typedef struct node_t {
  int m;
  int n;
  int k;
  int h;
  double xh;
  double ak;
  double bk;
  double *min;
  double *max;
  double **a;
  double *b;
  double *x;
  double *c;
  double z;
} node_t;

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

typedef struct list_t {
	node_t *node;
	struct list_t *next;
} list_t;

int init(simplex_t *s, int m, int n, double** a, double* b, double* c, double* x, double y, int* var);
void pivot(simplex_t *s, int row, int col);
void prepare(simplex_t *s, int k);
int initial(simplex_t *s, int m, int n, double** a, double* b, double* c, double* x, double y, int* var);
int select_nonbasic(simplex_t *s);
double xsimplex(int m, int n, double** a, double* b, double* c, double* x, double y, int* var, int h);
double simplex(int m, int n, double** a, double* b, double* c, double* x, double y);
double** make_matrix(int m, int n);
int main(int argc, char *argv[]);

node_t* initial_node(int m, int n, double **a, double *b, double *c);
node_t* extend(node_t *p, int m, int n, double **a, double *b, double *c, int k, double ak, double bk);
int is_integer(double *xp);
int integer(node_t *p);
void bound(node_t *p, list_t *h, double *zp, double *x);
/*int isfinite(double *x);*/
int branch(node_t *q, double z);
void succ(node_t *p, list_t *h, int m, int n, double **a, double *b, double *c, int k, double ak, double bk, double *zp, double *x);
double intopt(int m, int n, double **a, double *b, double *c, double *x);

// end of header

const double epsilon = 1e-9;

void* xmalloc(size_t size) {
  void* ptr = malloc(size);

  if (ptr == NULL) {
    fprintf(stderr, "out of memory\n");
    exit(1);
  }

  return ptr;
}

void bp(void) {}

node_t* initial_node(int m, int n, double **a, double *b, double *c) {
  int i, j;
  node_t *p = calloc(1, sizeof(node_t));
  p->a = make_matrix(m+1, n+1);
  p->b = calloc(m+1, sizeof(double));
  p->c = calloc(n+1, sizeof(double));
  p->x = calloc(m+n+1, sizeof(double));
  p->min = calloc(n, sizeof(double));
  p->max = calloc(n, sizeof(double));
  p->m = m;
  p->n = n;

  // Copy a, b, and c to p
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++)
      p->a[i][j] = a[i][j];

    p->b[i] = b[i];
  }

  for (i = 0; i < n; i++)
    p->c[i] = c[i];

  for (i = 0; i < n; i++) {
    p->min[i] = -INFINITY;
    p->max[i] = INFINITY;
  }
  return p;
} 

node_t* extend(node_t *p, int m, int n, double **a, double *b, double *c, int k, double ak, double bk) {
  node_t *q = calloc(1, sizeof(node_t));
  int i, j;
  q->k = k;
  q->ak = ak;
  q->bk = bk;

  if (ak > 0 && isfinite(p->max[k]))
    q->m = p->m;
  else if (ak < 0 && p->min[k] > 0)
    q->m = p->m;
  else 
    q->m = p->m + 1;

  q->n = p->n;
  q->h = -1;
  q->a = make_matrix(q->m+1, q->n+1);
  q->b = calloc(q->m+1, sizeof(double));
  q->c = calloc(q->n+1, sizeof(double));
  q->x = calloc(q->n+1, sizeof(double));
  q->min = calloc(q->n, sizeof(double));
  q->max = calloc(q->n, sizeof(double));
  
  // copy values to q
  for (i = 0; i < n; i++) {
    q->min[i] = p->min[i];
    q->max[i] = p->max[i];
    q->c[i] = c[i];
  }
  for (i = 0; i < m; i++) {
    for (j=0; j<n; j++)
			q->a[i][j] = a[i][j];

    q->b[i] = b[i];
  }
  
  if (ak > 0) {
    if (q->max[k] == INFINITY || bk < q->max[k])
      q->max[k] = bk;
  } else if (q->min[k] == -INFINITY || -bk > q->min[k])
    q->min[k] = -bk;

  for (i = m, j = 0; j < n; j++) {
    if (q->min[j] > -INFINITY) {
      q->a[i][j] = -1;
      q->b[i] = -q->min[j];
      i += 1;
    }
    if (q->max[j] < INFINITY) {
      q->a[i][j] = 1;
      q->b[i] = q->max[j];
      i += 1;
    }
  }
  return q;
}

int is_integer(double *xp) {
  double x = *xp;
  double r = round(x);
  if (fabs(r-x) < epsilon) {
    *xp = r;
    return 1;
  }
  return 0;
}

int integer(node_t *p) {
  int i;
  for (i = 0; i < p->n; i++)
    if (!is_integer(&p->x[i]))
      return 0;
  return 1;
}

int list_size_of(list_t *h) {
  int count = 0;
  list_t *current = h;

  while(current != NULL) {
    count ++;
    current = current->next;
  }
  return count;
}

list_t* list_at_index(list_t *h, int i) {
  list_t *current = h;
  for (int j = 0; j < i; j++)
    current = current->next;
  return current;
}

void list_add(list_t **h, node_t *q) {
  // create new element with node
  list_t *new_node = (list_t *)malloc(sizeof(list_t));
  new_node->node = q;
  new_node->next = NULL;

  // If the list is empty, the new node becomes the head
  if (*h == NULL) {
      *h = new_node;
      return;
  }

  // loop through list to get the last element
  list_t *current = *h;
  while(current->next != NULL) {
    current = current->next;
  }
  // set our new element as next
  current->next = new_node;
}

// Function to delete nodes based on the given condition
void delete_nodes(list_t **h, double pz) {
    // Check if the list is empty
    if (*h == NULL) {
        return; // Nothing to delete
    }

    // Check first element
    while (*h != NULL && (*h)->node->z < pz) {
        // Store old head
        list_t *old = *h;
        // Move head to the next node
        *h = (*h)->next; 
        // Free memory of the old head
        free(old->node); // Free the node data if necessary
        free(old);       // Free the list node
    }

    // Now, *h points to the first valid node or NULL
    list_t *current = *h;

    // Iterate through the list and delete nodes
    while (current != NULL && current->next != NULL) {
        if (current->next->node->z < pz) {
            // Store the next node
            list_t *next = current->next;
            // Remove the next node from the list
            current->next = next->next; 
            // Free memory
            free(next->node); // Free the node data if necessary
            free(next);       // Free the list node
        } else {
            // Move to the next node only if no deletion was made
            current = current->next; 
        }
    }
}

node_t* pop(list_t *h) {
  node_t *node = h->node;
  h = h->next;
  return node;
}

void bound(node_t *p, list_t *h, double *zp, double *x) {
  int i;
  if (p->z > *zp) {
    *zp = p->z;
    for (i = 0; i < p->n; i++) 
      x[i] = p->x[i];
    delete_nodes(&h, p->z);
  }
}

int branch(node_t *q, double z) {
  double min, max;
  int h;

  if (q->z < z) 
    return 0;

  for (h = 0; h < q->n; h++) {
    if (!is_integer(&q->x[h])) {
      if (q->min[h] = -INFINITY)
        min = 0;
      else 
        min = q->min[h];
      max = q->max[h];
      if (floor(q->x[h])<min || ceil(q->x[h])>max)
        continue;
      q->h = h;
      q->xh = q->x[h];
      //delete a,b,c,x of q (or recycle in another way??)
      free(q->a);
      free(q->b);
      free(q->c);
      free(q->x);
      return 1;
    }
  }
}


void succ(node_t *p, list_t *h, int m, int n, double **a, double *b, double *c, int k, double ak, double bk, double *zp, double *x){
  node_t *q = extend(p, m, n, a, b, c, k, ak, bk);
  if (q == NULL) {
    if (integer(q)) {
      bound(q, h, zp, x);
    }
    else if (branch(q, *zp)){
      list_add(&h, q);
      return;
    }
  }
  free(q); // onödigt då jag allokerar minnet i list_add
}

double intopt(int m, int n, double **a, double *b, double *c, double *x) {
  int i;
  node_t *p = initial_node(m, n, a, b, c);
  list_t *h = malloc(sizeof(list_t));
  h->node = p;
  h->next = NULL;
  double z = -INFINITY;
  p->z = simplex(p->m, p->n, p->a, p->b, p->c, p->x, 0);

  if (integer(p) || !isfinite(p->z)) {
    z = p->z;
    if (integer(p)) {
      // copy px to x 
      for (i = 0; i<n; i++)
				x[i] = p->x[i];
      // delete p 
      free(p);
      return z;
    }
  }
  branch(p, z);
  while (list_size_of(h) != 0) {
    // take p from h????
    p = pop(h);
		succ(p, h, m, n, a, b, c, p->h, 1, floor(p->xh), &z, x);
		succ(p, h, m, n, a, b, c, p->h, -1, -ceil(p->xh), &z, x);
    free(p);
  }
  if (z = -INFINITY)
    return NAN;
  else 
    return z;
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
    s->var = xmalloc((m+n+1)*sizeof(int));
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

int  local_array[10];
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

  for (i = 0; i < 11; i++)
    local_array[i] = i;

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

  s->x = calloc(m+n, sizeof(double)); // doesn't free the old?
  s->c = calloc(n, sizeof(double));
  s->c[n-1] = -1;
  s->n = n;
  pivot(s, k, n-1);
}

int initial(simplex_t *s, int m, int n, double** a, double* b, double* c, double* x, double y, int* var) {
  int i, j, k = 0;
  double w;
  k = init(s, m, n, a, b, c, x, y, var);

  if (b[k] >= 0) {
    return 1; // feasible
  }
  prepare(s, k);
  n = s->n;
  s->y = xsimplex(m, n, s->a, s->b, s->c, s->x, 0, s->var, 1);

  for (i = 0; i < m+n; i++) {
    if (s->var[i] == m+n-1) {
      if (fabs(s->x[i])>epsilon) {
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
      if (fabs(s->a[i-n][k])>fabs(s->a[i-n][j]))
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
  }
  free(s->c);
  s->c = c;
  s->y = y;
  for (k = n-1;k < n+m-1; k++)
    s->var[k] = s->var[k+1];
  n = s->n = s->n-1;
  double* t = calloc(n, sizeof(double*));

  for (k = 0; k < n; k++) {
    for (j = 0; j < n; j++) {
      if (k == s->var[j]) {
        t[j] = t[j] + s->c[k];
        goto next_k;
      }
    }
    for (j = 0; j < m; j++) {
      if (s->var[n+j] == k)
        break;
    }
    s->y = s->y + s->c[k]*s->b[j];
    for (i = 0; i < n; i++) 
      t[i] = t[i] - s->c[k]*s->a[j][i];

  next_k:;
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

  if (h == 0) {
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
  /*printf("Solution y=%lf\n", simplex(m, n, a, b, c, x, 0));*/
  printf("Solution y=%lf\n", intopt(m, n, a, b, c, x));

  for (int i = 0; i < m; i += 1)
    free(a[i]);
  
  free(a);
  free(b);
  free(c);
  free(x);

  return 0;
}

