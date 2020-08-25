#define MY_UTIL

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include "math_extra.h"
#include "myutil.h"

void project_pt_plane(const double* q,
                                        const double* p, const double* n,
                                        double* q_proj, double &d)
{
  double dot, ans[3], n_p[3];
  n_p[0] = n[0]; n_p[1] = n[1]; n_p[2] = n[2];
  MathExtra::sub3(q, p, ans);
  dot = MathExtra::dot3(ans, n_p);
  MathExtra::scale3(dot, n_p);
  MathExtra::sub3(q, n_p, q_proj);
  MathExtra::sub3(q, q_proj, ans);
  d = MathExtra::len3(ans);
}

void distance_bt_edges(const double* x1,
                  const double* x2, const double* x3, const double* x4,
                  double* h1, double* h2, double& t1, double& t2, double& r)
{
  double u[3],v[3],n[3],dot;

  // set the default returned values

  t1 = -2;
  t2 = 2;
  r = 0;

  // find the edge unit directors and their dot product

  MathExtra::sub3(x2, x1, u);
  MathExtra::norm3(u);
  MathExtra::sub3(x4, x3, v);
  MathExtra::norm3(v);
  dot = MathExtra::dot3(u,v);
  dot = fabs(dot);

  // check if two edges are parallel
  // find the two ends of the overlapping segment, if any

  if (fabs(dot - 1.0) < 0.01) {

    double s1,s2,x13[3],x23[3],x13h[3];
    double t13,t23,t31,t41,x31[3],x41[3];
    t13=t23=t31=t41=0.0;
    
    MathExtra::sub3(x1,x3,x13); // x13 = x1 - x3
    MathExtra::sub3(x2,x3,x23); // x23 = x2 - x3

    s1 = MathExtra::dot3(x13,v);
    x13h[0] = x13[0] - s1*v[0];
    x13h[1] = x13[1] - s1*v[1];
    x13h[2] = x13[2] - s1*v[2];
    r = MathExtra::len3(x13h);
    
    // x13 is the projection of x1 on x3-x4

    x13[0] = x3[0] + s1*v[0];
    x13[1] = x3[1] + s1*v[1];
    x13[2] = x3[2] + s1*v[2];

    // x23 is the projection of x2 on x3-x4

    s2 = MathExtra::dot3(x23,v);
    x23[0] = x3[0] + s2*v[0];
    x23[1] = x3[1] + s2*v[1];
    x23[2] = x3[2] + s2*v[2];
    
    // find the fraction of the projection points on the edges

    if (fabs(x4[0] - x3[0]) > 0)
      t13 = (x13[0] - x3[0])/(x4[0] - x3[0]);
    else if (fabs(x4[1] - x3[1]) > 0)
      t13 = (x13[1] - x3[1])/(x4[1] - x3[1]);
    else if (fabs(x4[2] - x3[2]) > 0)
      t13 = (x13[2] - x3[2])/(x4[2] - x3[2]);

    if (fabs(x4[0] - x3[0]) > 0)
      t23 = (x23[0] - x3[0])/(x4[0] - x3[0]);
    else if (fabs(x4[1] - x3[1]) > 0)
      t23 = (x23[1] - x3[1])/(x4[1] - x3[1]);
    else if (fabs(x4[2] - x3[2]) > 0)
      t23 = (x23[2] - x3[2])/(x4[2] - x3[2]);

    if (fabs(x23[0] - x13[0]) > 0)
      t31 = (x3[0] - x13[0])/(x23[0] - x13[0]);
    else if (fabs(x23[1] - x13[1]) > 0)
      t31 = (x3[1] - x13[1])/(x23[1] - x13[1]);
    else if (fabs(x23[2] - x13[2]) > 0)
      t31 = (x3[2] - x13[2])/(x23[2] - x13[2]);

    // x31 is the projection of x3 on x1-x2

    x31[0] = x1[0] + t31*(x2[0] - x1[0]);
    x31[1] = x1[1] + t31*(x2[1] - x1[1]);
    x31[2] = x1[2] + t31*(x2[2] - x1[2]);

    if (fabs(x23[0] - x13[0]) > 0)
      t41 = (x4[0] - x13[0])/(x23[0] - x13[0]);
    else if (fabs(x23[1] - x13[1]) > 0)
      t41 = (x4[1] - x13[1])/(x23[1] - x13[1]);
    else if (fabs(x23[2] - x13[2]) > 0)
      t41 = (x4[2] - x13[2])/(x23[2] - x13[2]);

    // x41 is the projection of x4 on x1-x2

    x41[0] = x1[0] + t41*(x2[0] - x1[0]);
    x41[1] = x1[1] + t41*(x2[1] - x1[1]);
    x41[2] = x1[2] + t41*(x2[2] - x1[2]);

    // determine two ends from the overlapping segments

    int n1 = 0;
    int n2 = 0;
    if (t13 >= 0 && t13 <= 1) {
      h1[0] = x1[0];
      h1[1] = x1[1];
      h1[2] = x1[2];
      h2[0] = x13[0];
      h2[1] = x13[1];
      h2[2] = x13[2];
      t1 = 0;
      t2 = t13;
      n1++;
      n2++;
    }
    if (t23 >= 0 && t23 <= 1) {
      if (n1 == 0) {
        h1[0] = x2[0];
        h1[1] = x2[1];
        h1[2] = x2[2];
        h2[0] = x23[0];
        h2[1] = x23[1];
        h2[2] = x23[2];
        t1 = 1;
        t2 = t23;
        n1++;
        n2++;
      } else {
        h1[0] = (x1[0]+x2[0])/2;
        h1[1] = (x1[1]+x2[1])/2;
        h1[2] = (x1[2]+x2[2])/2;
        h2[0] = (x13[0]+x23[0])/2; 
        h2[1] = (x13[1]+x23[1])/2; 
        h2[2] = (x13[2]+x23[2])/2;
        t1 = 0.5;
        t2 = (t13+t23)/2;
        n1++;
        n2++;
      }
    }

    if (n1 == 0 && n2 == 0) {
      if (t31 >= 0 && t31 <= 1) {
        h1[0] = x31[0];
        h1[1] = x31[1];
        h1[2] = x31[2];
        h2[0] = x3[0];
        h2[1] = x3[1];
        h2[2] = x3[2];
        t1 = t31;
        t2 = 0;
        n1++;
        n2++;
      }
      if (t41 >= 0 && t41 <= 1) {
        if (n1 == 0) {
          h1[0] = x41[0];
          h1[1] = x41[1];
          h1[2] = x41[2];
          h2[0] = x4[0];
          h2[1] = x4[1];
          h2[2] = x4[2];
          t1 = t41;
          t2 = 1;
          n1++;
          n2++;
        } else {
          h1[0] = (x31[0]+x41[0])/2;
          h1[1] = (x31[1]+x41[1])/2;
          h1[2] = (x31[2]+x41[2])/2;
          h2[0] = (x3[0]+x4[0])/2; 
          h2[1] = (x3[1]+x4[1])/2; 
          h2[2] = (x3[2]+x4[2])/2;
          t1 = (t31+t41)/2;
          t2 = 0.5;
          n1++;
          n2++;
        }
      }
    }   

    // if n1 == 0 and n2 == 0 at this point,
    // which means no overlapping segments bt two parallel edges,
    // return the default values of t1 and t2

    return;

  } 

  // find the vector n perpendicular to both edges
 
  MathExtra::cross3(u, v, n);
  MathExtra::norm3(n);

  // find the intersection of the line (x3,x4) and the plane (x1,x2,n)
  // s = director of the line (x3,x4)
  // n_p = plane normal vector of the plane (x1,x2,n)

  double s[3], n_p[3];
  MathExtra::sub3(x4, x3, s);
  MathExtra::sub3(x2, x1, u);
  MathExtra::cross3(u, n, n_p);
  MathExtra::norm3(n_p);

  // solve for the intersection between the line and the plane

  double m[3][3], invm[3][3], p[3], ans[3];
  m[0][0] = -s[0];
  m[0][1] = u[0];
  m[0][2] = n[0];

  m[1][0] = -s[1];
  m[1][1] = u[1];
  m[1][2] = n[1];

  m[2][0] = -s[2];
  m[2][1] = u[2];
  m[2][2] = n[2];

  MathExtra::sub3(x3, x1, p);
  MathExtra::invert3(m, invm);
  MathExtra::matvec(invm, p, ans);

  t2 = ans[0];
  h2[0] = x3[0] + s[0] * t2;
  h2[1] = x3[1] + s[1] * t2;
  h2[2] = x3[2] + s[2] * t2;

  project_pt_plane(h2, x1, n, h1, r);

  if (fabs(x2[0] - x1[0]) > 0)
    t1 = (h1[0] - x1[0])/(x2[0] - x1[0]);
  else if (fabs(x2[1] - x1[1]) > 0)
    t1 = (h1[1] - x1[1])/(x2[1] - x1[1]);
  else if (fabs(x2[2] - x1[2]) > 0)
    t1 = (h1[2] - x1[2])/(x2[2] - x1[2]);
}


double brute_force_distance(Segment S1, Segment S2) {
  double r = MAXFLOAT;
  double nseg = 2000;
  double dt = 1/nseg;
  Vector v1 = S1.P1 - S1.P0;
  Vector v2 = S2.P1 - S2.P0;
  for (int i = 0; i <= nseg; i++) {
    Point h1 = i * dt * v1 + S1.P0;
    for (int j = 0; j <= nseg; j++) {
      Point h2 = j * dt * v2 + S2.P0;
      r = fmin(r, d(h1, h2));
    }
  }
  return r;
}


int main() {

    srand(time(NULL));

    for (int i = 1; i < 2; i++) {
      Point P1, P2;
      Point P3, P4;
      Segment S1(P1, P2), S2(P3, P4);
      S1.print();
      S2.print();

      double r2 = dist3D_Segment_to_Segment(S1, S2);
      printf("Dan Sunday's solution is: %f\n", r2);

      double h1[3], h2[3];
      double t1, t2;
      double r = -1;
      double x1[3] = {P1.x, P1.y, P1.z}, x2[3] = {P2.x, P2.y, P2.z};
      double x3[3] = {P3.x, P3.y, P3.z}, x4[3] = {P4.x, P4.y, P4.z};
      distance_bt_edges(x1, x2, x3, x4, h1, h2, t1, t2, r);
      printf("LAMMPS version's solution is : %f\n", r);

      double brute_r = brute_force_distance(S1, S2);
      printf("brute force solution is: %f\n", brute_r);
    }
    return 0;
}