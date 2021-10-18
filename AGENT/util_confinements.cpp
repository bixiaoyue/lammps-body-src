#include "util_confinements.h"

int test_arma()
{
  mat A(4, 5, fill::randu);
  mat B(4, 5, fill::randu);

  cout << A*B.t() << endl;

  return 0;
}