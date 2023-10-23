// #pragma once
#include "xBurgers.hpp"
using namespace std;

void CADSolver::Boundary_uu(void)
{
	// #pragma omp barrier

	switch (BoundaryScheme)
	{
	case 0:
		Boundary_uu_periodic();
		break;
	case 1:
		Boundary_uu_Expand();
		break;
	case 2:
		Boundary_uu_reflection_x();
		break;
	case 24:
		Boundary_uu_Viscous_shock();
		break;

	case 22:
		Boundary_uu_DoubleMach();
		break;

	case 23:
		Boundary_uu_High_jet();
		break;
	case 26:
		Boundary_uu_Couette_flow();
		break;

	default:
		cout << "Boundary_uu is not avaliable for this Solver." << endl;
		break;
	}

	// #pragma omp barrier
}

void CADSolver::Boundary_uu_periodic(void)
{
	mydouble uuu;
	// #pragma omp parallel for num_threads(NumberThreads)
	for (long jCell = 0; jCell < nCellYo_2L; jCell++) // X方向
	{
		for (long iCell = 0; iCell < nb; iCell++)
			for (long nn = 0; nn < 4; nn++)
			{

				uuu = D_un[nCell + iCell][jCell][nn];
				D_un[iCell][jCell][nn] = uuu;
			}

		for (long iCell = nCell_L; iCell < nCell_2L; iCell++)
			for (long nn = 0; nn < 4; nn++)
			{

				uuu = D_un[iCell - nCell][jCell][nn];
				D_un[iCell][jCell][nn] = uuu;
			}
	}

	// // #pragma omp parallel for num_threads(NumberThreads)
	for (long iCell = 0; iCell < nCell_2L; iCell++) // y方向
	{
		for (long jCell = 0; jCell < nb; jCell++)
			for (long nn = 0; nn < 4; nn++)
			{

				uuu = D_un[iCell][jCell + nCellYo][nn];
				D_un[iCell][jCell][nn] = uuu;
			}

		for (long jCell = nCellYo_L; jCell < nCellYo_2L; jCell++)
			for (long nn = 0; nn < 4; nn++)
			{

				uuu = D_un[iCell][jCell - nCellYo][nn];
				D_un[iCell][jCell][nn] = uuu;
			}
	}
}

void CADSolver::Boundary_uu_Expand(void)
{
	mydouble uuu;
	// // #pragma omp parallel for num_threads(NumberThreads)
	for (long jCell = 0; jCell < nCellYo_2L; jCell++)
	{
		for (long iCell = 0; iCell < nb; iCell++)
			for (long nn = 0; nn < 4; nn++)
			{

				uuu = D_un[nb][jCell][nn];
				D_un[iCell][jCell][nn] = uuu;
			}

		for (long iCell = nCell_L; iCell < nCell_2L; iCell++)
			for (long nn = 0; nn < 4; nn++)
			{

				uuu = D_un[nCell_L - 1][jCell][nn];
				D_un[iCell][jCell][nn] = uuu;
			}
	}

	// // #pragma omp parallel for num_threads(NumberThreads)
	for (long iCell = 0; iCell < nCell_2L; iCell++)
	{
		for (long jCell = 0; jCell < nb; jCell++)
			for (long nn = 0; nn < 4; nn++)
			{

				uuu = D_un[iCell][nb][nn];
				D_un[iCell][jCell][nn] = uuu;
			}

		for (long jCell = nCellYo_L; jCell < nCellYo_2L; jCell++)
			for (long nn = 0; nn < 4; nn++)
			{

				uuu = D_un[iCell][nCellYo_L - 1][nn];
				D_un[iCell][jCell][nn] = uuu;
			}
	}
}

void CADSolver::Boundary_uu_reflection_x(void)
{
	mydouble uuu;
	// // #pragma omp parallel for num_threads(NumberThreads)
	for (long jCell = 0; jCell < nCellYo_2L; jCell++)
	{
		for (long iCell = 0; iCell < nb; iCell++)
		{
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_un[nb * 2 - iCell - 1][jCell][nn];
				D_un[iCell][jCell][nn] = uuu;
			}
			uuu = D_un[nb * 2 - iCell - 1][jCell][1];
			D_un[iCell][jCell][1] = -uuu;
		}

		for (long iCell = nCell_L; iCell < nCell_2L; iCell++)
		{
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_un[nCell_L * 2 - 1 - iCell][jCell][nn];
				D_un[iCell][jCell][nn] = uuu;
			}
			uuu = D_un[nCell_L * 2 - 1 - iCell][jCell][1];
			D_un[iCell][jCell][1] = -uuu;
		}
	}

	// // #pragma omp parallel for num_threads(NumberThreads)
	for (long iCell = 0; iCell < nCell_2L; iCell++)
	{
		for (long jCell = 0; jCell < nb; jCell++)
			for (long nn = 0; nn < 4; nn++)
			{

				uuu = D_un[iCell][nb][nn];
				D_un[iCell][jCell][nn] = uuu;
			}

		for (long jCell = nCellYo_L; jCell < nCellYo_2L; jCell++)
			for (long nn = 0; nn < 4; nn++)
			{

				uuu = D_un[iCell][nCellYo_L - 1][nn];
				D_un[iCell][jCell][nn] = uuu;
			}
	}
}

void CADSolver::Boundary_uu_Viscous_shock(void)
{

#pragma omp parallel for num_threads(NumberThreads)
	for (long jCell = 0; jCell < nCellYo_2L; jCell++)
	{
		mydouble uuu, density, velocity_x, velocity_y, pressure;
		long nb1;
		for (long iCell = 0; iCell < nb; iCell++) //_________lift__________
		{
			nb1 = nb * 2 - iCell - 1;
			density = D_un[nb1][jCell][0], velocity_x = D_un[nb1][jCell][1] / density;
			velocity_y = D_un[nb1][jCell][2] / density;
			pressure = (Gamma - 1.0) * (D_un[nb1][jCell][3] - 0.5 * density * (velocity_x * velocity_x + velocity_y * velocity_y));

			D_un[iCell][jCell][0] = density;
			D_un[iCell][jCell][1] = 0.0;
			D_un[iCell][jCell][2] = 0.0;
			D_un[iCell][jCell][3] = pressure / Gamma_1;
		}

		for (long iCell = nCell_L; iCell < nCell_2L; iCell++) //_________right__________
		{

			nb1 = nCell_L * 2 - 1 - iCell;
			density = D_un[nb1][jCell][0], velocity_x = D_un[nb1][jCell][1] / density;
			velocity_y = D_un[nb1][jCell][2] / density;
			pressure = (Gamma - 1.0) * (D_un[nb1][jCell][3] - 0.5 * density * (velocity_x * velocity_x + velocity_y * velocity_y));
			D_un[iCell][jCell][0] = density;
			D_un[iCell][jCell][1] = 0.0;
			D_un[iCell][jCell][2] = 0.0;
			D_un[iCell][jCell][3] = pressure / Gamma_1;
		}
	}

#pragma omp parallel for num_threads(NumberThreads)
	for (long iCell = 0; iCell < nCell_2L; iCell++)
	{
		mydouble uuu, density, velocity_x, velocity_y, pressure;
		long nb1;
		for (long jCell = 0; jCell < nb; jCell++) //_________bottom__________
		{
			nb1 = nb * 2 - jCell - 1;
			density = D_un[iCell][nb1][0], velocity_x = D_un[iCell][nb1][1] / density;
			velocity_y = D_un[iCell][nb1][2] / density;
			pressure = (Gamma - 1.0) * (D_un[iCell][nb1][3] - 0.5 * density * (velocity_x * velocity_x + velocity_y * velocity_y));

			// D_un[iCell][jCell][0] = density;
			// D_un[iCell][jCell][1] = 0.0;
			// D_un[iCell][jCell][2] = 0.0;
			// D_un[iCell][jCell][3] = pressure / Gamma_1;
			D_un[iCell][jCell][0] = density;
			D_un[iCell][jCell][1] = -density * velocity_x;
			D_un[iCell][jCell][2] = 0.0;
			D_un[iCell][jCell][3] = D_un[iCell][nb1][3] - 0.5 * density * velocity_y * velocity_y;
		}

		for (long jCell = nCellYo_L; jCell < nCellYo_2L; jCell++) //_________top__________
		{

			nb1 = nCellYo_L * 2 - 1 - jCell;
			density = D_un[iCell][nb1][0], velocity_x = D_un[iCell][nb1][1] / density;
			velocity_y = D_un[iCell][nb1][2] / density;
			pressure = (Gamma - 1.0) * (D_un[iCell][nb1][3] - 0.5 * density * (velocity_x * velocity_x + velocity_y * velocity_y));

			D_un[iCell][jCell][0] = density;
			D_un[iCell][jCell][1] = D_un[iCell][nb1][1];
			D_un[iCell][jCell][2] = 0.0;
			D_un[iCell][jCell][3] = D_un[iCell][nb1][3] - 0.5 * density * (velocity_y * velocity_y);
		}
	}
}

void CADSolver::Boundary_uu_Couette_flow(void)
{
	// #pragma omp parallel for num_threads(NumberThreads)
	// 	for (long jCell = 0; jCell < nCellYo_2L; jCell++) // X方向
	// 	{
	// 		mydouble uuu;
	// 		for (long iCell = 0; iCell < nb; iCell++)
	// 			for (long nn = 0; nn < 4; nn++)
	// 			{

	// 				uuu = D_un[nCell + iCell][jCell][nn];
	// 				D_un[iCell][jCell][nn] = uuu;
	// 			}

	// 		for (long iCell = nCell_L; iCell < nCell_2L; iCell++)
	// 			for (long nn = 0; nn < 4; nn++)
	// 			{

	// 				uuu = D_un[iCell - nCell][jCell][nn];
	// 				D_un[iCell][jCell][nn] = uuu;
	// 			}
	// 	}

#pragma omp parallel for num_threads(NumberThreads)
	for (long jCell = 0; jCell < nCellYo_2L; jCell++) // X方向
	{
		mydouble pp, dd, uu, vv, yp;
		mydouble H, u1, mu, ka, T0, T1, Ta;
		for (long iCell = 0; iCell < nb; iCell++)
		{
			yp = D_Coord[iCell][jCell][1];
			H = 2., u1 = 0.1;
			uu = yp * u1 / H, vv = 0.;
			T0 = 0.8, T1 = 0.85;
			Ta = T0 + yp * (T1 - T0) / H + Pr * Gamma_1 * u1 * u1 * yp / (2. * H) * (1 - yp / H);
			pp = 1. / Gamma, dd = Gamma * pp / Ta;
			D_un[iCell][jCell][0] = dd;
			D_un[iCell][jCell][1] = dd * uu;
			D_un[iCell][jCell][2] = dd * vv;
			D_un[iCell][jCell][3] = pp / Gamma_1 + 0.5 * dd * (uu * uu + vv * vv);
		}
		for (long iCell = nCell_L; iCell < nCell_2L; iCell++)
		{
			yp = D_Coord[iCell][jCell][1];
			H = 2., u1 = 0.1;
			uu = yp * u1 / H, vv = 0.;
			T0 = 0.8, T1 = 0.85;
			Ta = T0 + yp * (T1 - T0) / H + Pr * Gamma_1 * u1 * u1 * yp / (2. * H) * (1 - yp / H);
			pp = 1. / Gamma, dd = Gamma * pp / Ta;
			D_un[iCell][jCell][0] = dd;
			D_un[iCell][jCell][1] = dd * uu;
			D_un[iCell][jCell][2] = dd * vv;
			D_un[iCell][jCell][3] = pp / Gamma_1 + 0.5 * dd * (uu * uu + vv * vv);
		}
	}

#pragma omp parallel for num_threads(NumberThreads)
	for (long iCell = 0; iCell < nCell_2L; iCell++)
	{
		mydouble pp, dd, uu, vv, yp;
		mydouble H, u1, mu, ka, T0, T1, Ta;
		for (long jCell = 0; jCell < nb + 1; jCell++) //_________bottom__________
		{

			H = 2., u1 = 0.1;
			uu = 0., vv = 0.;
			T0 = 0.8, T1 = 0.85;
			// Ta = 0.8;
			// T0 + yp *(T1 - T0) / H + Pr *Gamma_1 *u1 *u1 *yp / (2. * H) * (1 - yp / H);

			pp = 1. / Gamma, dd = Gamma * pp / T0;

			D_un[iCell][jCell][0] = dd;
			D_un[iCell][jCell][1] = dd * uu;
			D_un[iCell][jCell][2] = dd * vv;
			D_un[iCell][jCell][3] = pp / Gamma_1 + 0.5 * dd * (uu * uu + vv * vv);
		}

		for (long jCell = nCellYo_L - 1; jCell < nCellYo_2L; jCell++) //_________top__________
		{
			yp = D_Coord[iCell][nCellYo_L - 1][1];
			H = 2., u1 = 0.1;
			uu = yp * u1 / H, vv = 0.;
			T0 = 0.8, T1 = 0.85;
			Ta = T0 + yp * (T1 - T0) / H + Pr * Gamma_1 * u1 * u1 * yp / (2. * H) * (1 - yp / H);
			pp = 1. / Gamma, dd = Gamma * pp / Ta;
			D_un[iCell][jCell][0] = dd;
			D_un[iCell][jCell][1] = dd * uu;
			D_un[iCell][jCell][2] = dd * vv;
			D_un[iCell][jCell][3] = pp / Gamma_1 + 0.5 * dd * (uu * uu + vv * vv);
		}
	}
}

void CADSolver::Boundary_uu_DoubleMach(void)
{

#pragma omp parallel for num_threads(NumberThreads)
	for (long jCell = 0; jCell < nCellYo_2L; jCell++)
	{
		mydouble alfa, dd, uu, vv, pp;
		alfa = 60. / 180. * PI_Number;
		dd = 8.0, uu = 8.25 * sin(alfa), vv = -8.25 * cos(alfa), pp = 116.5;

		for (long iCell = 0; iCell < nb; iCell++) //_________lift__inlet________
		{

			D_un[iCell][jCell][0] = dd;
			D_un[iCell][jCell][1] = dd * uu;
			D_un[iCell][jCell][2] = dd * vv;
			D_un[iCell][jCell][3] = pp / Gamma_1 + 0.5 * dd * (uu * uu + vv * vv);
		}

		for (long iCell = nCell_L; iCell < nCell_2L; iCell++) //_________right_out_________
		{

			for (long nn = 0; nn < 4; nn++)
			{
				D_un[iCell][jCell][nn] = D_un[nCell_L - 1][jCell][nn];
			}
		}
	}

#pragma omp parallel for num_threads(NumberThreads)
	for (long iCell = 0; iCell < nCell_2L; iCell++)
	{
		mydouble alfa, dd, uu, vv, pp, x0, x1, xp, yp, density, velocity_x, velocity_y, pressure;
		long nb1;

		alfa = 60. / 180. * PI_Number;
		dd = 8.0, uu = 8.25 * sin(alfa), vv = -8.25 * cos(alfa), pp = 116.5;
		x0 = 1. / 6.;

		for (long jCell = 0; jCell < nb; jCell++) //_________bottom__________
		{
			xp = D_Coord[iCell][jCell][0];
			if (xp < x0)
			{
				D_un[iCell][jCell][0] = dd;
				D_un[iCell][jCell][1] = dd * uu;
				D_un[iCell][jCell][2] = dd * vv;
				D_un[iCell][jCell][3] = pp / Gamma_1 + 0.5 * dd * (uu * uu + vv * vv);
			}
			else
			{

				nb1 = nb * 2 - jCell - 1;
				density = D_un[iCell][nb1][0], velocity_x = D_un[iCell][nb1][1] / density;
				velocity_y = D_un[iCell][nb1][2] / density;

				D_un[iCell][jCell][0] = density;
				D_un[iCell][jCell][1] = density * velocity_x;
				D_un[iCell][jCell][2] = 0.0;
				D_un[iCell][jCell][3] = D_un[iCell][nb1][3] - 0.5 * density * velocity_y * velocity_y;
			}
		}

		for (long jCell = nCellYo_L; jCell < nCellYo_2L; jCell++) //_________top__________
		{
			yp = D_Coord[iCell][jCell][1];
			xp = D_Coord[iCell][jCell][0];
			x1 = x0 + yp / tan(alfa) + 10. * TimeNow / sin(alfa);

			if (xp < x1)
			{
				dd = 8.0, uu = 8.25 * sin(alfa), vv = -8.25 * cos(alfa), pp = 116.5;
			}
			else
			{
				dd = 1.4, uu = 0.0, vv = 0.0, pp = 1.0;
			}

			D_un[iCell][jCell][0] = dd;
			D_un[iCell][jCell][1] = dd * uu;
			D_un[iCell][jCell][2] = dd * vv;
			D_un[iCell][jCell][3] = pp / Gamma_1 + 0.5 * dd * (uu * uu + vv * vv);
		}
	}
}

void CADSolver::Boundary_uu_High_jet(void)
{

#pragma omp parallel for num_threads(NumberThreads)
	for (long jCell = 0; jCell < nCellYo_2L; jCell++)
	{
		mydouble uuu, yp, uu, vv, pp, dd;
		for (long iCell = 0; iCell < nb; iCell++)
		{
			yp = D_Coord[iCell][jCell][1];
			if (abs(yp) <= 0.050001)
			{
				uu = 30.0, vv = 0., pp = 0.4127, dd = 5.;
			}
			else
			{
				uu = 0.0, vv = 0., pp = 0.4127, dd = 5.;
			}

			D_un[iCell][jCell][0] = dd;
			D_un[iCell][jCell][1] = dd * uu;
			D_un[iCell][jCell][2] = dd * vv;
			D_un[iCell][jCell][3] = pp / Gamma_1 + 0.5 * dd * (uu * uu + vv * vv);
		}

		for (long iCell = nCell_L; iCell < nCell_2L; iCell++)
			for (long nn = 0; nn < 4; nn++)
			{

				uuu = D_un[nCell_L - 1][jCell][nn];
				D_un[iCell][jCell][nn] = uuu;
			}
	}

#pragma omp parallel for num_threads(NumberThreads)
	for (long iCell = 0; iCell < nCell_2L; iCell++)
	{
		mydouble uuu;
		for (long jCell = 0; jCell < nb; jCell++)
			for (long nn = 0; nn < 4; nn++)
			{

				uuu = D_un[iCell][nb][nn];
				D_un[iCell][jCell][nn] = uuu;
			}

		for (long jCell = nCellYo_L; jCell < nCellYo_2L; jCell++)
			for (long nn = 0; nn < 4; nn++)
			{

				uuu = D_un[iCell][nCellYo_L - 1][nn];
				D_un[iCell][jCell][nn] = uuu;
			}
	}
}

void CADSolver::Boundary_dutxy_all(void)
{
	switch (BoundaryScheme)
	{
	case 0:
		Boundary_dutxy_periodic();
		break;
	case 1:
		Boundary_dutxy_expand();
		break;

	case 2:
		Boundary_dutxy_reflection_x();
		break;

	default:
		cout << "BoundaryScheme is not avaliable for this Solver." << endl;
		break;
	}
}

void CADSolver::Boundary_dutxy_periodic(void)
{
	mydouble uuu;
	// #pragma omp parallel for num_threads(NumberThreads) // x方向
	for (long jCell = 0; jCell < nCellYo_2L; jCell++)
	{
		for (long iCell = 0; iCell < nb; iCell++)
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_utx0[nCell + iCell][jCell][nn];
				D_utx0[iCell][jCell][nn] = uuu;
			}

		for (long iCell = nCell_L; iCell < nCell_2L; iCell++)
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_utx0[iCell - nCell][jCell][nn];
				D_utx0[iCell][jCell][nn] = uuu;
			}
	}

	// #pragma omp parallel for num_threads(NumberThreads) // y方向
	for (long iCell = 0; iCell < nCell_2L; iCell++)
	{
		for (long jCell = 0; jCell < nb; jCell++)
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_uty0[iCell][jCell + nCellYo][nn];
				D_uty0[iCell][jCell][nn] = uuu;
			}

		for (long jCell = nCellYo_L; jCell < nCellYo_2L; jCell++)
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_uty0[iCell][jCell - nCellYo][nn];
				D_uty0[iCell][jCell][nn] = uuu;
			}
	}
}

void CADSolver::Boundary_dutxy_expand(void)
{
	mydouble uuu;
	// #pragma omp parallel for num_threads(NumberThreads) // x方向
	for (long jCell = 0; jCell < nCellYo_2L; jCell++)
	{
		for (long iCell = 0; iCell < nb; iCell++)
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_utx0[nb][jCell][nn];
				D_utx0[iCell][jCell][nn] = uuu;
			}

		for (long iCell = nCell_L; iCell < nCell_2L; iCell++)
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_utx0[nCell_L - 1][jCell][nn];
				D_utx0[iCell][jCell][nn] = uuu;
			}
	}

	// #pragma omp parallel for num_threads(NumberThreads) // y方向
	for (long iCell = 0; iCell < nCell_2L; iCell++)
	{
		for (long jCell = 0; jCell < nb; jCell++)
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_uty0[iCell][nb][nn];
				D_uty0[iCell][jCell][nn] = uuu;
			}

		for (long jCell = nCellYo_L; jCell < nCellYo_2L; jCell++)
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_uty0[iCell][nCellYo_L - 1][nn];
				D_uty0[iCell][jCell][nn] = uuu;
			}
	}
}

void CADSolver::Boundary_dutxy_reflection_x(void)
{
	mydouble uuu;
	// #pragma omp parallel for num_threads(NumberThreads) // x方向
	for (long jCell = 0; jCell < nCellYo_2L; jCell++)
	{
		for (long iCell = 0; iCell < nb; iCell++)
		{
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_utx0[nb * 2 - iCell - 1][jCell][nn];
				D_utx0[iCell][jCell][nn] = uuu;
			}
			uuu = D_utx0[nb * 2 - iCell - 1][jCell][1];
			D_utx0[iCell][jCell][1] = uuu;
		}
		for (long iCell = nCell_L; iCell < nCell_2L; iCell++)
		{
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_utx0[nCell_L * 2 - 1 - iCell][jCell][nn];
				D_utx0[iCell][jCell][nn] = uuu;
			}
			uuu = D_utx0[nCell_L * 2 - 1 - iCell][jCell][1];
			D_utx0[iCell][jCell][1] = uuu;
		}
	}

	// #pragma omp parallel for num_threads(NumberThreads) // y方向
	for (long iCell = 0; iCell < nCell_2L; iCell++)
	{
		for (long jCell = 0; jCell < nb; jCell++)
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_uty0[iCell][nb][nn];
				D_uty0[iCell][jCell][nn] = uuu;
			}

		for (long jCell = nCellYo_L; jCell < nCellYo_2L; jCell++)
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_uty0[iCell][nCellYo_L - 1][nn];
				D_uty0[iCell][jCell][nn] = uuu;
			}
	}
}

void CADSolver::Boundary_dutxy_LWA_periodic(void)
{
	mydouble uuu;

	// #pragma omp parallel for num_threads(NumberThreads) // x方向
	for (long jCell = 0; jCell < nCellYo_2L; jCell++)
	{
		for (long iCell = 0; iCell < nb; iCell++)
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_fxat[nCell + iCell][jCell][nn][0];
				D_fxat[iCell][jCell][nn][0] = uuu;
			}

		for (long iCell = nCell_L; iCell < nCell_2L; iCell++)
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_fxat[iCell - nCell][jCell][nn][0];
				D_fxat[iCell][jCell][nn][0] = uuu;
			}
	}

	// #pragma omp parallel for num_threads(NumberThreads) // y方向
	for (long iCell = 0; iCell < nCell_2L; iCell++)
	{
		for (long jCell = 0; jCell < nb; jCell++)
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_fyat[iCell][jCell + nCellYo][nn][0];
				D_fyat[iCell][jCell][nn][0] = uuu;
			}

		for (long jCell = nCellYo_L; jCell < nCellYo_2L; jCell++)
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_fyat[iCell][jCell - nCellYo][nn][0];
				D_fyat[iCell][jCell][nn][0] = uuu;
			}
	}
}
