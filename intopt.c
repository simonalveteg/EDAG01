#include <stdio.h>
#include <stdlib.h>

double** make_matrix(int m, int n) {
  double** a;
  int i;

  a = calloc(m, sizeof(double*));
  for (int i = 0; i < m; i += 1) {
    a[i] = calloc(n-1, sizeof(double));
  }
  return a;
}

int main(int argc, char *argv[])
{
  int a;
  int b;
  double** c;
  
  scanf("%d %d", &a, &b);
  c = make_matrix(a, b);

  for (int i = 0; i < a; i += 1) {
    for (int j = 0; j < b; j += 1) {
      scanf("%lf", &c[i][j]);
    }
  }
  for (int i = 0; i < a - 1; i++) {
    if (i == 0) {
      printf("max z = ");
      for (int j = 0; j < b; j++) {
        printf("%+5.1lfx%d", c[i][j], j);
      }
      printf("\n");
      continue;
    }
    for (int j = 0; j < b; j++) {
      printf("%+10.3lfx%d", c[i][j], j);
    }
    printf("<=%10.3lf", c[a-1][i-1]);
    printf("\n");
  }
  free(c);

  return 0;
}

