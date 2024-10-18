#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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

typedef struct list_t list_t;
typedef struct link_t link_t;
list_t* create_list(node_t* first); 
static link_t* create_link(void* data); 
void insert(list_t* head, void* data);
void free_list(list_t* head); 
int remove_link_at_index(list_t* head, int index); 
int list_size(list_t* head); 
void* get_data_at_index(list_t* head, int index);
void delete_nodes(list_t* h, double zp); 
node_t* pop(list_t* head);

int init(simplex_t *s, int m, int n, double** a, double* b, double* c, double* x, double y, int* var);
void pivot(simplex_t *s, int row, int col);
void prepare(simplex_t *s, int k);
int initial(simplex_t *s, int m, int n, double** a, double* b, double* c, double* x, double y, int* var);
int select_nonbasic(simplex_t *s);
double xsimplex(int m, int n, double** a, double* b, double* c, double* x, double y, int* var, int h);
double simplex(int m, int n, double** a, double* b, double* c, double* x, double y);
double** make_matrix(int m, int n);
int main(void);

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

void free_node(node_t **p)
{
  if (*p==NULL) return;
	if ((*p)->a!=NULL){
		for (int i = 0; i < (*p)->m+1; i++) {
      free((*p)->a[i]);
      (*p)->a[i] = NULL;
		}
		free((*p)->a);
    (*p)->a = NULL;
	}
	if ((*p)->b!=NULL) {
		free((*p)->b);
    (*p)->b = NULL;
  }
	if ((*p)->c!=NULL) {
		free((*p)->c);
    (*p)->c = NULL;
  }
	if ((*p)->x!=NULL) {
		free((*p)->x);
    (*p)->x = NULL;
  }
	free((*p)->min);
  (*p)->min = NULL;
	free((*p)->max);
  (*p)->max = NULL;
	free(*p);
  *p = NULL;
}

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
  node_t *q = xmalloc(sizeof(*q));
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
  q->min = calloc(n, sizeof(double));
  q->max = calloc(n, sizeof(double));
  
  // copy values to q
  for (i = 0; i < n; i++) {
    q->min[i] = p->min[i];
    q->max[i] = p->max[i];
    q->c[i] = c[i];
  }
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++)
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
      i++;
    }
    if (q->max[j] < INFINITY) {
      q->a[i][j] = 1;
      q->b[i] = q->max[j];
      i++;
    }
  }
  return q;
}

int is_integer(double *xp) {
  double x = *xp;
  double r = lround(x);
  if (fabs(r-x) < epsilon) {
    *xp = r;
    return 1;
  } else return 0;
}

int integer(node_t *p) {
  int i;
  for (i = 0; i < p->n; i++) {
    if (!is_integer(&(p->x[i])))
      return 0;
  } return 1;
}


void bound(node_t *p, list_t *h, double *zp, double *x) {
  int i;
  if (p->z > *zp) {
    *zp = p->z;
    for (i = 0; i < p->n; i++) 
      x[i] = p->x[i];

    // remove and delete all nodes q in h with q->z < p->z
    delete_nodes(h, p->z);
  }
}

int branch(node_t *q, double z) {
  double min, max;
  int h, i;

  if (q->z < z) 
    return 0;

  for (h = 0; h < q->n; h++) {
    if (!is_integer(&(q->x[h]))) {
      if (q->min[h] == -INFINITY)
        min = 0;
      else 
        min = q->min[h];
      max = q->max[h];
      if (floor(q->x[h]) < min || ceil(q->x[h]) > max)
        continue;
      q->h = h;
      q->xh = q->x[h];
      //delete a,b,c,x of q (or recycle in another way??)
      if (q->a != NULL) {
        for (i = 0; i < q->m + 1; i++) {
          free(q->a[i]);
          q->a[i] = NULL;
        }
        free(q->a);
        q->a = NULL;
      }
      if (q->b != NULL) {
        free(q->b);
        q->b = NULL;
      }
      if (q->c != NULL) {
        free(q->c);
        q->c = NULL;
      }
      if (q->x != NULL) {
        free(q->x);
        q->x = NULL;
      }
      return 1;
    }
  }
}


void succ(node_t *p, list_t *h, int m, int n, double **a, double *b, double *c, int k, double ak, double bk, double *zp, double *x){
  node_t *q = extend(p, m, n, a, b, c, k, ak, bk);

  if (q == NULL) 
    return;

  q->z = simplex(q->m, q->n, q->a, q->b, q->c, q->x, 0);

  if (isfinite(q->z)) {
    if (integer(q)) {
      bound(q, h, zp, x);
    } else if (branch(q, *zp)){
      insert(h, q);
      return;
    }
  }
  free_node(&q);
}

double intopt(int m, int n, double **a, double *b, double *c, double *x) {
  int i;
  node_t *p = initial_node(m, n, a, b, c);
  list_t *h = create_list(p);
  double z = -INFINITY;

  p->z = simplex(p->m, p->n, p->a, p->b, p->c, p->x, 0);

  if (integer(p) || !isfinite(p->z)) {
    z = p->z;
    if (integer(p)) {
      // copy px to x 
      for (i = 0; i<n; i++)
				x[i] = p->x[i];
    }
    // delete p BENJAMIN GÖR DETTA I IF SATSEN OVAN??
    free_node(&p);
    /*free(h);*/
    /*h = NULL;*/
    free_list(h);
    return z;
  }
  
  branch(p, z);

  while (list_size(h) != 0) {
    int size = list_size(h);
    node_t *q = pop(h);
		succ(q, h, m, n, a, b, c, q->h, 1, floor(q->xh), &z, x);
		succ(q, h, m, n, a, b, c, q->h, -1, -ceil(q->xh), &z, x);
    free_node(&q);
  }
  free(h);
  h = NULL;
  if (z == -INFINITY)
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
  for (i = 0; i < m; i += 1) {
    a[i] = calloc(n, sizeof(double));
  }
  return a;
}

int main(void)
{
  int i, j;
  int m = 0;
  int n = 0;
  double** a;
  double* b; 
  double* c;
  
  scanf("%d %d", &m, &n);
  a = make_matrix(m, n);
  b = calloc(m, sizeof(double));
  c = calloc(n, sizeof(double));

  for (i = 0; i < n; i++) {
    scanf("%lf", &c[i]);
  }
  for (i = 0; i < m; i += 1) {
    for (j = 0; j < n; j += 1) {
      scanf("%lf", &a[i][j]);
    }
  }
  printf("max z = %5.1lf x0", c[0]);
  for(i=1; i<n; i+=1) {
    printf(" + %lf x%d", c[i], i);
  }
  printf("\n");
  for (i = 0; i < m; i++) {
    /*printf("%+10.3lfx0 %+10.3lfx1 <= %10.3lf\n", a[i][0], a[i][1], a[m-1][i]);*/
    for (j = 0; j < n; j++) {
      printf("%+10.3lf x%d", a[i][j], j);
    }
    scanf("%lf", &b[i]);
    printf(" <=%10.3lf\n", b[i]);
  }
  double *x = calloc(n+m+1, sizeof(double));
  for (i = 0; i < n; i++)
    x[i] = 0;
  /*printf("Solution y=%lf\n", simplex(m, n, a, b, c, x, 0));*/
  printf("Solution y=%lf\n", intopt(m, n, a, b, c, x));

  for (i = 0; i < m; i += 1) {
    free(a[i]);
    a[i] = NULL;
  }
  
  free(a);
  a = NULL;
  free(b);
  b = NULL;
  free(c);
  c = NULL;
  free(x);
  x = NULL;

  return 0;
}


struct list_t {
    link_t* first;
    link_t* last;
};

struct link_t {
    link_t *next;
    void *data;
};

list_t* create_list(node_t *first) {
    list_t *head;
    link_t *link = create_link(first);

    head = xmalloc(sizeof(list_t));
    head->first = head->last  = link;
    return head;
}

static link_t* create_link(void* data) {
    link_t *link;

    link = xmalloc(sizeof(link_t));
    link->next = NULL;
    link->data = data;
    return link;
}

void insert(list_t* head, void* data) {
    link_t* link;
    link = create_link(data);

    if (head->first == NULL) {
        head->first = link;
    } else {
        head->last->next = link;
    }
    head->last = link;
}

void free_list(list_t* head) {
    link_t* p;
    link_t* q;

    p = head->first;
    if (head != NULL) {
        free(head);
        head = NULL;
    }

    while(p != NULL) {
        q = p->next;
        if (p != NULL) {
            free(p);
            p = NULL;
        }
        p = q;
    }
}

int list_size(list_t* head) {
  int size = 0;
  link_t* current = head->first;

  while (current != NULL) {
    size++;
    current = current->next;
  }

  return size;
}

void delete_nodes(list_t* h, double zp) {
    if (h == NULL) return;  // Nothing to delete

    link_t* current = h->first;
    link_t* prev = NULL;

    while (current != NULL) {
        node_t* node = current->data;

        if (node->z < zp) {
            // We are deleting this node

            // Check if current is the last node before setting next_node
            link_t* next_node = current->next;

            if (prev == NULL) {
                // If deleting the first node
                h->first = next_node;
            } else {
                // Skip over the deleted node
                prev->next = next_node;
            }

            // Check if the node being deleted is the last node
            if (current == h->last) {
                h->last = prev;  // Set h->last to prev, or NULL if no prev exists
            }

            // Free memory
            free(current->data);
            free(current);

            current = next_node;
        } else {
            // Move to the next node if no deletion was made
            prev = current;
            current = current->next;
        }
    }

    // If the list is now empty, make sure h->last is NULL
    if (h->first == NULL) {
        h->last = NULL;
    }
}

node_t* pop(list_t* head) {
  if (head == NULL || head->first == NULL) return NULL; // List is empty.

  link_t* prev = NULL;
  link_t* current = head->first;

  // Case 1: List has only one element
  if (head->first == head->last) {
    node_t* data = current->data;  // Get the data before freeing
    free(current);  // Free the only node
    head->first = head->last = NULL;  // Reset list pointers
    return data;
  }

  // Case 2: List has more than one element
  // Traverse to the last element
  while (current != head->last) {
    prev = current;
    current = current->next;
  }

  // At this point, prev is the second-to-last node, and current is the last node

  // Detach the last node safely
  prev->next = NULL;
  head->last = prev;  // Update head->last to the previous node

  node_t* data = current->data;  // Get the data before freeing
  free(current);  // Free the last node

  return data;
}

