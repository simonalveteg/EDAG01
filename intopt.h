
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

int init(simplex_t *s, int m, int n, double** a, double* b, double* c, double* x, double y, int* var);

void pivot(simplex_t *s, int row, int col);

void prepare(simplex_t *s, int k);

int initial(simplex_t *s, int m, int n, double** a, double* b, double* c, double* x, double y, int* var);

int select_nonbasic(simplex_t *s);

double xsimplex(int m, int n, double** a, double* b, double* c, double* x, double y, int* var, int h);

double simplex(int m, int n, double** a, double* b, double* c, double* x, double y);

double** make_matrix(int m, int n);

int main(int argc, char *argv[]);


