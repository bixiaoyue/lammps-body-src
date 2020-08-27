/* ----------------------------------------------------------------------
   Author: Farzan Beroz, Princeton University
   Github: https://github.com/farzanb/Agent-based-biofilm-model
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Modified by Changhao Li (czl478@psu.edu) for modeling of bacteria film
   Last modified: 08/13/2020
------------------------------------------------------------------------- */

#include "util_fibrin.h"

// File name for output
string my_name;

Cell::Cell(double mx, double my, double mz, double mnx, double mny, double mnz, double ml)
{
	x = mx;
	y = my;
	z = mz;
	nx = mnx;
	ny = mny;
	nz = mnz;

	dx = noise_level * (rand() / (static_cast<double>(RAND_MAX)) - 0.5);
	dy = noise_level * (rand() / (static_cast<double>(RAND_MAX)) - 0.5);
	dz = noise_level * (rand() / (static_cast<double>(RAND_MAX)) - 0.5);
	dnx = noise_level * (rand() / (static_cast<double>(RAND_MAX)) - 0.5);
	dny = noise_level * (rand() / (static_cast<double>(RAND_MAX)) - 0.5);
	dnz = noise_level * (rand() / (static_cast<double>(RAND_MAX)) - 0.5);

	current_nz = nz;

	l = ml;

	Lfi = Lf;

	confine = 0;
	fpop = 0;

	lineage = "";

	ai = A;
}

vector<Cell *> cells;

my6Vec cell_surface_gforce(Cell *cell_1, double kn)
{
	// Return the surface torque acting on cell 1

	double t;

	double E_1 = kn;
	double E_2 = E_1;

	double A_1 = E_1 / 300;
	double A_2 = A_1;

	double rix = cell_1->get_x();
	double riy = cell_1->get_y();
	double riz = cell_1->get_z();
	double rinx = cell_1->get_nx();
	double riny = cell_1->get_ny();
	double rinz = cell_1->get_nz();

	double nzi = cell_1->get_nz();

	double check = 1;

	double l = cell_1->get_l();
	double z = cell_1->get_z();

	if (rinz < 0)
	{
		rinz = -rinz;
		check = -1;
	}

	if (nzi < -1)
	{
		nzi = -1;
	}
	if (nzi > 1)
	{
		nzi = 1;
	}

	if (nzi > 0)
	{
		t = asin(nzi);
	}
	else
	{
		t = asin(-nzi);
	}

	if (rinz < -1)
	{
		rinz = -1;
	}
	if (rinz > 1)
	{
		rinz = 1;
	}

	double h1 = z + l * nzi / 2.0;
	double h2 = z - l * nzi / 2.0;

	if (h1 > h2)
	{
		h1 = h2;
	}

	// printf("%f, %f\n", h1, h2);

	my6Vec F;
	F.x = 0;
	F.y = 0;
	F.z = 0;
	F.nx = 0;
	F.ny = 0;
	F.nz = 0;

	if (h1 >= R || abs(h1 - R) < 1e-6)
	{
		// No contact
	}
	else if (h1 > (R - l * sin(t)))
	{
		// Partial contact
		F.z -= (A_2 * M_PI * rinz) / R - (4 * E_2 * sqrt(R) * rinz * pow(R + (l * rinz) / 2. - riz, 1.5)) / 3. -
			   (E_1 * (1 - pow(rinz, 2)) * pow(-R - (l * rinz) / 2. + riz, 2)) / rinz +
			   (A_1 * (1 - pow(rinz, 2)) * sqrt(1 - (-(l * rinz) / 2. + riz) / R)) / (R * rinz);
		F.nz -= -(A_2 * l * M_PI * rinz) / (2. * R) - (A_2 * M_PI * (R + (l * rinz) / 2. - riz)) / R +
				(2 * E_2 * l * sqrt(R) * rinz * pow(R + (l * rinz) / 2. - riz, 1.5)) / 3. +
				(8 * E_2 * sqrt(R) * pow(R + (l * rinz) / 2. - riz, 2.5)) / 15. +
				(E_1 * l * (1 - pow(rinz, 2)) * pow(-R - (l * rinz) / 2. + riz, 2)) / (2. * rinz) +
				(2 * E_1 * pow(-R - (l * rinz) / 2. + riz, 3)) / 3. +
				(E_1 * (1 - pow(rinz, 2)) * pow(-R - (l * rinz) / 2. + riz, 3)) / (3. * pow(rinz, 2)) -
				(A_1 * l * (1 - pow(rinz, 2)) * sqrt(1 - (-(l * rinz) / 2. + riz) / R)) / (2. * R * rinz) +
				(4 * A_1 * pow(1 - (-(l * rinz) / 2. + riz) / R, 1.5)) / 3. +
				(2 * A_1 * (1 - pow(rinz, 2)) * pow(1 - (-(l * rinz) / 2. + riz) / R, 1.5)) / (3. * pow(rinz, 2));
	}
	else if (h1 > 0)
	{
		// Full contact
		if (abs(t) < 1e-6)
		{
			// linear version. Contact forces are linearized for small elevation angle
			// to avoid numerical instabilities
			F.z -= -2 * E_1 * l * R + (A_1 * l) / (2. * pow(R, 1.5) * sqrt(R - riz)) + 2 * E_1 * l * riz +
				   pow(rinz, 2) * (2 * E_1 * l * R + (A_1 * pow(l, 3)) / (64. * pow(R, 1.5) * pow(R - riz, 2.5)) +
								   (2 * E_2 * l * pow(R, 2.5)) / (3. * pow(R - riz, 1.5)) - (A_1 * l) / (2. * pow(R, 1.5) * sqrt(R - riz)) -
								   (8 * E_2 * l * pow(R, 1.5)) / (3. * sqrt(R - riz)) - 2 * E_1 * l * riz -
								   (4 * E_2 * l * pow(R, 1.5) * riz) / (3. * pow(R - riz, 1.5)) + (8 * E_2 * l * sqrt(R) * riz) / (3. * sqrt(R - riz)) +
								   (2 * E_2 * l * sqrt(R) * pow(riz, 2)) / (3. * pow(R - riz, 1.5)));
			F.nz -= 2 * rinz * ((E_1 * pow(l, 3)) / 12. - (A_2 * l * M_PI) / R - E_1 * l * pow(R, 2) + (A_1 * pow(l, 3)) / (96. * pow(R, 1.5) * pow(R - riz, 1.5)) + (4 * E_2 * l * pow(R, 2.5)) / (3. * sqrt(R - riz)) + (A_1 * l * sqrt(R - riz)) / pow(R, 1.5) + 2 * E_1 * l * R * riz - (8 * E_2 * l * pow(R, 1.5) * riz) / (3. * sqrt(R - riz)) - E_1 * l * pow(riz, 2) + (4 * E_2 * l * sqrt(R) * pow(riz, 2)) / (3. * sqrt(R - riz)));
		}
		else
		{
			// nonlinear version
			F.z -= (E_1 * l * (1 - pow(rinz, 2)) * (3 * l * rinz + 6 * (-R - (l * rinz) / 2. + riz))) / 3. -
				   (2 * A_1 * (1 - pow(rinz, 2)) * (-(l * rinz) / (2. * sqrt(R - (l * rinz) / 2. - riz)) + sqrt(R - (l * rinz) / 2. - riz) - sqrt(R + (l * rinz) / 2. - riz) - (1 / (2. * sqrt(R - (l * rinz) / 2. - riz)) - 1 / (2. * sqrt(R + (l * rinz) / 2. - riz))) * (-R - (l * rinz) / 2. + riz))) / (3. * pow(R, 1.5) * rinz) -
				   (8 * E_2 * sqrt(R) * rinz * (2 * l * rinz * sqrt(R - (l * rinz) / 2. - riz) - 2 * (-sqrt(R - (l * rinz) / 2. - riz) + sqrt(R + (l * rinz) / 2. - riz)) * (-R - (l * rinz) / 2. + riz) - (1 / (2. * sqrt(R - (l * rinz) / 2. - riz)) - 1 / (2. * sqrt(R + (l * rinz) / 2. - riz))) * pow(-R - (l * rinz) / 2. + riz, 2) - (l * rinz * (-2 * R + l * rinz + 2 * (-(l * rinz) / 2. + riz))) / (2. * sqrt(R - (l * rinz) / 2. - riz)))) / 15.;
			F.nz -= (-2 * A_2 * l * M_PI * rinz) / R - (2 * A_1 * (1 - pow(rinz, 2)) * ((l * (-sqrt(R - (l * rinz) / 2. - riz) + sqrt(R + (l * rinz) / 2. - riz))) / 2. - (pow(l, 2) * rinz) / (4. * sqrt(R - (l * rinz) / 2. - riz)) + l * sqrt(R - (l * rinz) / 2. - riz) - (l / (4. * sqrt(R - (l * rinz) / 2. - riz)) + l / (4. * sqrt(R + (l * rinz) / 2. - riz))) * (-R - (l * rinz) / 2. + riz))) / (3. * pow(R, 1.5) * rinz) + (4 * A_1 * (l * rinz * sqrt(R - (l * rinz) / 2. - riz) - (-sqrt(R - (l * rinz) / 2. - riz) + sqrt(R + (l * rinz) / 2. - riz)) * (-R - (l * rinz) / 2. + riz))) / (3. * pow(R, 1.5)) + (2 * A_1 * (1 - pow(rinz, 2)) * (l * rinz * sqrt(R - (l * rinz) / 2. - riz) - (-sqrt(R - (l * rinz) / 2. - riz) + sqrt(R + (l * rinz) / 2. - riz)) * (-R - (l * rinz) / 2. + riz))) / (3. * pow(R, 1.5) * pow(rinz, 2)) -
					(8 * E_2 * sqrt(R) * rinz * (l * (-sqrt(R - (l * rinz) / 2. - riz) + sqrt(R + (l * rinz) / 2. - riz)) * (-R - (l * rinz) / 2. + riz) - (l / (4. * sqrt(R - (l * rinz) / 2. - riz)) + l / (4. * sqrt(R + (l * rinz) / 2. - riz))) * pow(-R - (l * rinz) / 2. + riz, 2) - (pow(l, 2) * rinz * (-2 * R + l * rinz + 2 * (-(l * rinz) / 2. + riz))) / (4. * sqrt(R - (l * rinz) / 2. - riz)) + l * sqrt(R - (l * rinz) / 2. - riz) * (-2 * R + l * rinz + 2 * (-(l * rinz) / 2. + riz)))) / 15. -
					(8 * E_2 * sqrt(R) * (-((-sqrt(R - (l * rinz) / 2. - riz) + sqrt(R + (l * rinz) / 2. - riz)) * pow(-R - (l * rinz) / 2. + riz, 2)) + l * rinz * sqrt(R - (l * rinz) / 2. - riz) * (-2 * R + l * rinz + 2 * (-(l * rinz) / 2. + riz)))) / 15. +
					(E_1 * l * (1 - pow(rinz, 2)) * (-(pow(l, 2) * rinz) / 2. - 3 * l * (-R - (l * rinz) / 2. + riz) + l * (-3 * R + l * rinz + 3 * (-(l * rinz) / 2. + riz)))) / 3. -
					(2 * E_1 * l * rinz * (3 * pow(-R - (l * rinz) / 2. + riz, 2) + l * rinz * (-3 * R + l * rinz + 3 * (-(l * rinz) / 2. + riz)))) / 3.;
		}
	}
	else
	{
		cout << "ERROR: Torque is below the ground!" << endl;
		return F;
	}

	F.nz = check * F.nz;

	return F;
}

/* ----------------------------------------------------------------------
   compute surface damping force and torque for cell-substrate interaction
------------------------------------------------------------------------- */

my6Vec compute_surface_damping_force(Cell *cell_1, double *v, double *omega, double nu_1, double A)
{
	double x = cell_1->get_x();
	double y = cell_1->get_y();
	double z = cell_1->get_z();
	double nx = cell_1->get_nx();
	double ny = cell_1->get_ny();
	double nz = cell_1->get_nz();
	double L = cell_1->get_l();

	double d1 = R - (z - L / 2 * nz);
	double d2 = R - (z + L / 2 * nz);

	my6Vec damping;

	double quad_points[5] = {-0.906, -0.538, 0, 0.538, 0.906};
	double weights[5] = {0.237, 0.479, 0.569, 0.479, 0.237};

	// make sure d1 < d2
	if (d1 > d2)
	{
		double temp = 0;
		temp = d1;
		d1 = d2;
		d2 = temp;
		// nx = -nx; ny = -ny; nz = -nz;
	}

	if (d1 < 0 || abs(d1) < 1e-6)
	{
		// no contact
	}
	else
	{
		double fx[5] = {0}, fy[5] = {0}, fz[5] = {0}, mx[5] = {0}, my[5] = {0}, mz[5] = {0};

		// printf("omega: %e %e %e\n", omega[0], omega[1], omega[2]);
		// printf("nx: %f, ny: %f, nz: %f \n", nx, ny, nz);
		for (int i = 0; i < 5; i++)
		{
			double rr = quad_points[i] * L / 2;
			fx[i] = -L * nu_1 / R * contact_area_density(cell_1, rr) * (v[0] + cross(omega[0], omega[1], omega[2], rr * nx, rr * ny, rr * nz, 0));
			fy[i] = -L * nu_1 / R * contact_area_density(cell_1, rr) * (v[1] + cross(omega[0], omega[1], omega[2], rr * nx, rr * ny, rr * nz, 1));
			fz[i] = 0;
			// printf("fx: %f, fy: %f \n", fx[i], fy[i]);
			mx[i] = cross(rr * nx, rr * ny, rr * nz, fx[i], fy[i], fz[i], 0);
			my[i] = cross(rr * nx, rr * ny, rr * nz, fx[i], fy[i], fz[i], 1);
			mz[i] = cross(rr * nx, rr * ny, rr * nz, fx[i], fy[i], fz[i], 2);
			// printf("mx: %f, my: %f, mz: %f \n", mx[i], my[i], mz[i]);
		}

		damping.x = gaussian_quadrature(fx, weights, 5);
		damping.y = gaussian_quadrature(fy, weights, 5);
		damping.z = 0;
		damping.nx = gaussian_quadrature(mx, weights, 5);
		damping.ny = gaussian_quadrature(my, weights, 5);
		damping.nz = gaussian_quadrature(mz, weights, 5);

		// printf("Velocities: %f %f %f %f %f %f \n", v[0], v[1], v[2], omega[0], omega[1], omega[2]);
		// printf("Forces: %f %f %f %f %f %f\n", damping.x, damping.y, damping.z, damping.nx, damping.ny, damping.nz);
	}
	return damping;
}

/* ----------------------------------------------------------------------
   returns equivalent contact area density at relative location rr
------------------------------------------------------------------------- */

double contact_area_density(Cell *cell_1, double rr)
{
	double x = cell_1->get_x();
	double y = cell_1->get_y();
	double z = cell_1->get_z();
	double nx = cell_1->get_nx();
	double ny = cell_1->get_ny();
	double nz = cell_1->get_nz();
	double L = cell_1->get_l();

	double sin2_theta = pow(nz, 2);
	double cos2_theta = 1 - sin2_theta;

	double delta = R - (z + rr * nz);

	if (delta > 0)
		return sqrt(R) * cos2_theta * sqrt(delta) + M_PI * R * sin2_theta;
	else
		return 0;
}

/* ----------------------------------------------------------------------
   Add the cell-substrate force to the overall force vector
   cell_surface_gforce() computes generalized force on the normal vector
   needed to be translated to moment in space fixed frame
------------------------------------------------------------------------- */

void contact_forces_new(int ibody, Cell *cell, double *v, double *omega, double **f, double **torque, double kn, double cn, double ct, int shift_flag)
{
	my6Vec force_and_torque = cell_surface_gforce(cell, kn);
	double cell_force[3] = {force_and_torque.x, force_and_torque.y, force_and_torque.z};

	shift_vector(cell_force, 3, -shift_flag);

	f[ibody][0] += cell_force[0];
	f[ibody][1] += cell_force[1];
	f[ibody][2] += cell_force[2];

	// translate generalized force on z component of the normal vector
	double L = cell->get_l();
	double nx = cell->get_nx(), ny = cell->get_ny(), nz = cell->get_nz();
	double cos_phi = sqrt(nx * nx + ny * ny);
	double proj_xy = L * cos_phi;
	double cell_torque[3] = {cross(nx, ny, nz, 0, 0, force_and_torque.nz, 0) * cos_phi, cross(nx, ny, nz, 0, 0, force_and_torque.nz, 1) * cos_phi, 0};

	shift_vector(cell_torque, 3, -shift_flag);

	torque[ibody][0] += cell_torque[0];
	torque[ibody][1] += cell_torque[1];
	torque[ibody][2] += cell_torque[2];

	// compute surface damping force and momentum
	double nu_1 = cn;
	double A = ct;
	my6Vec damping = compute_surface_damping_force(cell, v, omega, nu_1, A);

	f[ibody][0] += damping.x;
	f[ibody][1] += damping.y;
	f[ibody][2] += damping.z;
	torque[ibody][0] += damping.nx;
	torque[ibody][1] += damping.ny;
	torque[ibody][2] += damping.nz;

	// printf("%e %e %e %e %e %e \n", force_and_torque.x, force_and_torque.y, force_and_torque.z, force_and_torque.nx, force_and_torque.ny, force_and_torque.nz);
	// printf("damping: %e %e %e %e %e %e \n", damping.x, damping.y, damping.z, damping.nx, damping.ny, damping.nz);
}

/* ----------------------------------------------------------------------
   shifting given vector by 1 position
   flag > 0: to right
   flag < 0: to left
   flag = 0: do nothing
------------------------------------------------------------------------- */

void shift_vector(double *v, int n, int flag)
{
	if (n > 0)
	{
		double temp[n];
		for (int i = 0; i < n; i++)
		{
			temp[i] = v[i];
		}
		if (flag > 0)
		{
			for (int i = 1; i < n; i++)
			{
				v[i] = temp[i - 1];
			}
			v[0] = temp[n - 1];
		}
		else if (flag < 0)
		{
			for (int i = 1; i < n; i++)
			{
				v[i - 1] = temp[i];
			}
			v[n - 1] = temp[0];
		}
	}
}

/* ----------------------------------------------------------------------
   Gaussian quadrature of n points
------------------------------------------------------------------------- */

double gaussian_quadrature(double *f, double *weights, int n)
{
	double result = 0;
	for (int i = 0; i < n; i++)
	{
		result += f[i] * weights[i];
	}
	return result;
}