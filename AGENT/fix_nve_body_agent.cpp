/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Modified by Changhao Li (czl478@psu.edu, changhaoli1997@gmail.com) for
   modeling of bacteria film. Last modified date: 06/30/2020
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include "body_rounded_polyhedron.h"
#include "fix_nve_body_agent.h"
#include "math_extra.h"
#include "atom.h"
#include "atom_vec_body.h"
#include "force.h"
#include "update.h"
#include "memory.h"
#include "error.h"


using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixNVEBodyAgent::FixNVEBodyAgent(LAMMPS *lmp, int narg, char **arg) : FixNVE(lmp, narg, arg)
{
  if (narg != 3 && narg != 6 && narg != 8 && narg != 10 && narg != 12 && narg!= 14)
    error->all(FLERR, "Invalid fix nve/body/agent command");

  // looking for args for body growth and proliferation
  // default values are 0, if 0 is set, proliferation and growth will be skiped
  growth_rate = 0;
  L_critical = 0;
  del_height = 1e6;
  for (int i = 0; i < narg - 1; i++)
  {
    if (strcmp(arg[i], "grow") == 0)
    {
      growth_rate = force->numeric(FLERR, arg[i + 1]);
      standard_dev = force->numeric(FLERR, arg[i + 2]);
    }
    if (strcmp(arg[i], "maxlen") == 0)
    {
      L_critical = force->numeric(FLERR, arg[i + 1]);
    }
    if (strcmp(arg[i], "damp") == 0)
    {
      nu_0 = force->numeric(FLERR, arg[i + 1]);
    }
    if (strcmp(arg[i], "noise") == 0)
    {
      noise_level = force->numeric(FLERR, arg[i + 1]);
    }
    if (strcmp(arg[i], "zdel") == 0)
    {
      del_height = force->numeric(FLERR, arg[i + 1]);
    }
  }

  // set up random seed and normal distribution
  srand(time(NULL));
  distribution = std::normal_distribution<double>(growth_rate, standard_dev);
}

/* ---------------------------------------------------------------------- */

void FixNVEBodyAgent::init()
{
  avec = (AtomVecBody *)atom->style_match("body");
  if (!avec)
    error->all(FLERR, "Fix nve/body requires atom style body");

  // check that all particles are bodies
  // no point particles allowed

  int *body = atom->body;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  list_growth_rate = std::vector<double>(nlocal * 2, 0);

  for (int i = 0; i < nlocal; i++)
  {
    list_growth_rate[i] = distribution(generator); // initiate all rates, including gel particles
    if (mask[i] & groupbit)
    {
      if (body[i] < 0)
      {
        error->one(FLERR, "Fix nve/body requires bodies");
      }
      list_growth_rate[i] = distribution(generator);
    }
  }

  FixNVE::init();
}

/* ---------------------------------------------------------------------- */

int FixNVEBodyAgent::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= INITIAL_INTEGRATE_RESPA;
  mask |= FINAL_INTEGRATE_RESPA;

  mask |= PRE_EXCHANGE;

  return mask;
}

/* ---------------------------------------------------------------------- */

void FixNVEBodyAgent::initial_integrate(int vflag)
{
  double dtfm;
  double omega[3];
  double *quat, *inertia;

  AtomVecBody::Bonus *bonus = avec->bonus;
  int *body = atom->body;
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **angmom = atom->angmom;
  double **torque = atom->torque;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;

  // bodies will grow during every timestep
  // here the growth law is exponential. dL = growth_rate * L in every timestep

  grow_all_body();

  // set timestep here since dt may have changed or come via rRESPA

  dtq = 0.5 * dtv;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
    {
      dtfm = dtf / rmass[i];

      inertia = bonus[body[i]].inertia;
      quat = bonus[body[i]].quat;
      // compute omega at the start of step from angmom at start of step and current q
      MathExtra::mq_to_omega(angmom[i], quat, inertia, omega);
      apply_damping_force(i, omega, f, torque);

      // adding random noise to force vector and moment vector

      // add_noise(v[i], omega, noise_level);
      add_noise(f[i], torque[i], noise_level);

      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];
      x[i][0] += dtv * v[i][0];
      x[i][1] += dtv * v[i][1];
      x[i][2] += dtv * v[i][2];

      // update angular momentum by 1/2 step

      angmom[i][0] += dtf * torque[i][0];
      angmom[i][1] += dtf * torque[i][1];
      angmom[i][2] += dtf * torque[i][2];

      // compute omega at 1/2 step from angmom at 1/2 step and current q
      // update quaternion a full step via Richardson iteration
      // returns new normalized quaternion

      MathExtra::mq_to_omega(angmom[i], quat, inertia, omega);

      MathExtra::richardson(quat, angmom[i], omega, inertia, dtq);
    }
}

/* ---------------------------------------------------------------------- */

void FixNVEBodyAgent::pre_exchange()
{
  double** x = atom->x;
  int *body = atom->body;
  int *type = atom->type;
  AtomVecBody::Bonus *bonus = avec->bonus;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;

  // change cell color(type), delete cell if cell is detached from substrate
  int i = 0;
  while (i < atom->nlocal)
  {
    if (atom->mask[i] & groupbit)
    {
    double vertices[6];
    body2space(bonus[body[i]].dvalue, bonus[body[i]].quat, vertices);
    double L = length(bonus[body[i]].dvalue);
    double nx = (vertices[3] - vertices[0]) / L;
    double ny = (vertices[4] - vertices[1]) / L;
    double nz = (vertices[5] - vertices[2]) / L;

    Cell icell = Cell(x[i][0], x[i][1], x[i][2], nx, ny, nz, L);

    if (cell_vertical(&icell))
    {
      type[i] = 2;
    } else {
      type[i] = 1;
    }

    if (out_of_z_plane(&icell, del_height))
    {
      printf("Cell %d is out of substrate, it will be removed! \n", i);
      if (true)
      {
        avec->deep_copy_body((atom->nlocal)-1, i, true);
      }
    }
    }
    i++;
  }
}

/* ---------------------------------------------------------------------- */

void FixNVEBodyAgent::final_integrate()
{
  double dtfm;

  AtomVecBody::Bonus *bonus = avec->bonus;
  int *body = atom->body;
  double **v = atom->v;
  double **f = atom->f;
  double **angmom = atom->angmom;
  double **torque = atom->torque;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
    {
      dtfm = dtf / rmass[i];

      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];

      angmom[i][0] += dtf * torque[i][0];
      angmom[i][1] += dtf * torque[i][1];
      angmom[i][2] += dtf * torque[i][2];
    }

  // cells will proliferate if length > critical length
  proliferate_all_body();
}

/* ----------------------------------------------------------------------
  compute damping force of ith body
---------------------------------------------------------------------- */

void FixNVEBodyAgent::apply_damping_force(int ibody, double *omega, double **f, double **torque)
{
  AtomVecBody::Bonus *bonus = avec->bonus;
  double *v = (atom->v)[ibody];
  int *body = atom->body;
  double L = length(bonus[body[ibody]].dvalue);
  double R = radius(bonus[body[ibody]].dvalue, 2);

  // adding damping force, applying on mass center
  f[ibody][0] += -nu_0 * (L+4/3*R) * v[0];
  f[ibody][1] += -nu_0 * (L+4/3*R) * v[1];
  f[ibody][2] += -nu_0 * (L+4/3*R) * v[2];

  // adding damping moment
  torque[ibody][0] += -1.0 / 6.0 * nu_0 * omega[0] * pow(L+4/3*R, 3);
  torque[ibody][1] += -1.0 / 6.0 * nu_0 * omega[1] * pow(L+4/3*R, 3);
  torque[ibody][2] += -1.0 / 6.0 * nu_0 * omega[2] * pow(L+4/3*R, 3);

  // printf("fx: %f, fy: %f, fz: %f \n", f[ibody][0], f[ibody][1], f[ibody][2]);
}

/* ----------------------------------------------------------------------
  transfer the body-frame relative displacement to space-frame
  specifially for rods
  don't forget to add the coords of the rod center!
---------------------------------------------------------------------- */

void FixNVEBodyAgent::body2space(double *coords, double *quat, double *space_coords)
{
  double p[3][3];
  MathExtra::quat_to_mat(quat, p);
  for (int m = 0; m < 2; m++)
  {
    MathExtra::matvec(p, &coords[3 * m], &space_coords[3 * m]);
  }
}

/* ----------------------------------------------------------------------
  grow single body, index is ibody
  growth_ratio is alpha in 2018 Nat. Phys. Paper
---------------------------------------------------------------------- */

void FixNVEBodyAgent::grow_single_body(int ibody, double growth_rate)
{
  AtomVecBody::Bonus *bonus = avec->bonus;
  double *rmass = atom->rmass;
  int *body = atom->body;
  double r = radius(bonus[body[ibody]].dvalue, 2);
  double L = length(bonus[body[ibody]].dvalue);

  if (bonus[body[ibody]].ndouble != 12) {printf("warning: ndouble != 12, herendouble = %d\n", bonus[body[ibody]].ndouble);}

  // add Gaussian random noise to growth rate: current sigma = 0.2<alpha>

  // dL/L = alpha * (4/3*R + L)/L * dt
  // dV/V = alpha * dt
  double length_ratio = 1 + 2 * dtf * growth_rate * (4.0 / 3.0 * r + L) / L;
  double growth_ratio = 1 + 2 * dtf * growth_rate;

  // update the coords of vertices, mass, rotation inertia
  for (int j = 0; j < 6; j++)
    bonus[body[ibody]].dvalue[j] *= length_ratio; // coords of vertices (in body frame)
  bonus[body[ibody]].dvalue[6+2] *= length_ratio; // enclosing radius (not rounded radius)
  rmass[ibody] *= growth_ratio;                   // mass ~ V
  for (int j = 0; j < 3; j++)
    bonus[body[ibody]].inertia[j] *= pow(growth_ratio, 3); // rotational inertia ~ m*L^2
}

/* ----------------------------------------------------------------------
  translate single body
---------------------------------------------------------------------- */

void FixNVEBodyAgent::translate_single_body(int ibody, double *vec)
{
  double **x = atom->x;
  for (int j = 0; j < 3; j++)
    x[ibody][j] += vec[j];
}

/* ----------------------------------------------------------------------
  return rounded radius for rod
---------------------------------------------------------------------- */

double FixNVEBodyAgent::radius(double *data, int nvert)
{
  return data[nvert * 3 + 2 + 1];
}

/* ----------------------------------------------------------------------
  return length for rod
---------------------------------------------------------------------- */

double FixNVEBodyAgent::length(double *data)
{
  return sqrt(pow(data[3] - data[0], 2) + pow(data[4] - data[1], 2) + pow(data[5] - data[2], 2));
}

/* ----------------------------------------------------------------------
  return the volume of the body
---------------------------------------------------------------------- */

double FixNVEBodyAgent::volume(double r, double L)
{
  double pi = 3.1415926;
  return 4 / 3 * pi * pow(r, 3) + pi * pow(r, 2) * L;
}

/* ----------------------------------------------------------------------
  if the cell is not in contact with z=height plane, return true
---------------------------------------------------------------------- */

bool FixNVEBodyAgent::out_of_z_plane(Cell* icell, double height)
{
  double relative_z = icell->get_z() - height;
  icell->set_z(relative_z);

  return !cell_contact(icell);
}

/* ----------------------------------------------------------------------
  return the rotational inertia of the body (unfinished)
---------------------------------------------------------------------- */

double *rot_inertia(double r, double L, double *mom)
{
  return NULL;
}

/* ----------------------------------------------------------------------
  grow pre-selected bodies by ratio (1+growth_rate) in single timestep
---------------------------------------------------------------------- */

void FixNVEBodyAgent::grow_all_body(double given_growth_ratio)
{
  AtomVecBody::Bonus *bonus = avec->bonus;
  int *body = atom->body;
  double **x = atom->x;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;

  if (growth_rate != 0)
  {
    for (int i = 0; i < nlocal; i++)
    {
      int nvertices = bonus[body[i]].ivalue[0];
      if ((mask[i] & groupbit) && (nvertices == 2))
      {
        grow_single_body(i, list_growth_rate[i]);
      }
    }
  }
}

/* ----------------------------------------------------------------------
  proliferate pre-selected bodies if their length exceeds L_critial
  L_critial = 0 means no proliferation
---------------------------------------------------------------------- */

void FixNVEBodyAgent::proliferate_all_body()
{
  AtomVecBody::Bonus *bonus = avec->bonus;
  int *body = atom->body;
  double **x = atom->x;
  double **angmom = atom->angmom;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;

  // loop to see if cells will proliferate

  if (L_critical != 0)
  {
    for (int i = 0; i < nlocal; i++)
    {
      int nvertices = bonus[body[i]].ivalue[0];
      if ((mask[i] & groupbit) && (nvertices == 2))
      {
        double *p = bonus[body[i]].dvalue;
        double L = length(p);    // length of the rod
        double r = radius(p, 2); // rounded radius of the rod
        if (L >= L_critical)
        {
          // find new center locations of two daughter bodies
          // c1 is the relative displacements of the first daughter body, c2 is the second

          double temp[6];
          for (int j = 0; j < 6; j++)
            temp[j] = bonus[body[i]].dvalue[j] * (L / 2 + r) / L;
          double cc[6];
          body2space(temp, bonus[body[i]].quat, cc);
          double *c1 = cc;
          double *c2 = cc + 3;

          double x_old[3];
          for (int j = 0; j < 3; j++)
            x_old[j] = x[i][j];

          // mother body -> the first daughter body

          double prolif_ratio = (L / 2 - r) / L;
          double extra_ratio = r / L;
          translate_single_body(i, c1);
          for (int j = 0; j < 6; j++)
            bonus[body[i]].dvalue[j] *= prolif_ratio;
          bonus[body[i]].dvalue[6+2] *= prolif_ratio;
          rmass[i] *= 0.50;
          for (int j = 0; j < 3; j++)
          {
            bonus[body[i]].inertia[j] *= pow(0.50, 3);
            angmom[i][j] *= pow(0.50, 3);
          }

          // insert the second daughter body

          avec->add_body(i);
          atom->avec->grow_reset();
          int new_body_index = atom->nlocal - 1; // append the new cell to the last

          // printf("body array: ");
          // for (int i = 0; i < atom->nlocal; i++) {
          //   printf("%d ", body[i]);
          // }
          // printf("\n");

          // update the second center
          for (int j = 0; j < 3; j++)
          {
            x[new_body_index][j] = x_old[j] + c2[j];
          }

          set_force(i, 0, 0, 0, 0, 0, 0);
          set_force(new_body_index, 0, 0, 0, 0, 0, 0);

          int p = new_body_index;
          // printf("mother cell %d: length %e, center %e %e %e, mass %e, inertia %e %e %e \n", i, length(bonus[body[i]].dvalue), x[i][0], x[i][1], x[i][2], rmass[i], bonus[body[i]].inertia[0], bonus[body[i]].inertia[1], bonus[body[i]].inertia[2]);
          // printf("mother force: %e %e %e %e %e %e \n", atom->f[i][0], atom->f[i][1], atom->f[i][2], atom->torque[i][0], atom->torque[i][1], atom->torque[i][2]);
          // printf("daughter cell %d: dindex: %d, iindex: %d, length %e, center %e %e %e, mass %e, inertia %e %e %e \n", p, bonus[body[p]].dindex, bonus[body[p]].iindex, length(bonus[body[p]].dvalue), x[p][0], x[p][1], x[p][2], rmass[i], bonus[body[i]].inertia[0], bonus[body[i]].inertia[1], bonus[body[i]].inertia[2]);
          // printf("daughter force: %e %e %e %e %e %e \n", atom->f[p][0], atom->f[p][1], atom->f[p][2], atom->torque[p][0], atom->torque[p][1], atom->torque[p][2]);

          // generate new random growth rate
          double new_rate = distribution(generator);
          if (new_rate > growth_rate + 3*standard_dev) new_rate = growth_rate + 3*standard_dev;
          if (new_rate < growth_rate - 3*standard_dev) new_rate = growth_rate - 3*standard_dev;
          list_growth_rate[new_body_index] = new_rate;
        }
      }
    }

    // update vector for random growth rate

    if (atom->nlocal >= list_growth_rate.size())
    {
      list_growth_rate.resize(2 * atom->nlocal, 0);
    }
  }
}

/* ----------------------------------------------------------------------
  Adding noise to force vector and moment vector
---------------------------------------------------------------------- */

void FixNVEBodyAgent::add_noise(double *f, double *mom, double given_noise_level)
{
  f[0] += 0 * noise_level * (rand() / (static_cast<double>(RAND_MAX)) - 0.5);
  f[1] += 0 * noise_level * (rand() / (static_cast<double>(RAND_MAX)) - 0.5);
  f[2] += 0 * noise_level * (rand() / (static_cast<double>(RAND_MAX)) - 0.5);
  mom[0] += given_noise_level * (rand() / (static_cast<double>(RAND_MAX)) - 0.5);
  mom[1] += given_noise_level * (rand() / (static_cast<double>(RAND_MAX)) - 0.5);
  mom[2] += given_noise_level * (rand() / (static_cast<double>(RAND_MAX)) - 0.5);
}

/* ----------------------------------------------------------------------
  Manually set force and moment values for given atom
---------------------------------------------------------------------- */

void FixNVEBodyAgent::set_force(int ibody, double fx, double fy, double fz, double tx, double ty, double tz)
{
  double* f = atom->f[ibody];
  double *torque = atom->torque[ibody];

  f[0] = fx;
  f[1] = fy;
  f[2] = fz;
  torque[0] = tx;
  torque[1] = ty;
  torque[2] = tz;
}