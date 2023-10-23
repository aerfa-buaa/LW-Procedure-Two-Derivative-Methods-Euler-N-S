// #pragma once
#include "xBurgers.hpp"
using namespace std;

void CADSolver::TimeIntegration_RK11(unsigned short iRK_Step)
{

	mydouble local;
	mydouble du, un0, un1, du0, un2;

	// #pragma omp parallel for num_threads(NumberThreads)
	for (long jCell = nb; jCell < nCellYo_L; jCell++)
	{
		mydouble local;
		for (long iCell = nb; iCell < nCell_L; iCell++)
		{
			for (short iVar = 0; iVar < nVar; iVar++)
			{

				du = D_du[iCell][jCell][iVar];
				un0 = D_un[iCell][jCell][iVar];

				local = un0 - DT * du;

				D_un[iCell][jCell][iVar] = local;
			}
		}
	}
}

void CADSolver::TimeIntegration(unsigned short iRK_Step)
{

	mydouble local;
	mydouble du, un0, un1, du0, un2;

	if (iRK_Step == 0)
	{
		// #pragma omp parallel for num_threads(NumberThreads)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{
			mydouble local;
			for (long iCell = nb; iCell < nCell_L; iCell++)
			{
				for (short iVar = 0; iVar < nVar; iVar++)
				{

					du = D_du[iCell][jCell][iVar];
					un0 = D_un[iCell][jCell][iVar];

					local = un0 - DT * du;

					D_un[iCell][jCell][iVar] = local;

					D_un0[iCell][jCell][iVar] = un0;
					D_du0[iCell][jCell][iVar] = du;
				}
			}
		}
	}
	else if (iRK_Step == 1)
	{
		// #pragma omp parallel for num_threads(NumberThreads)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{
			mydouble local;
			for (long iCell = nb; iCell < nCell_L; iCell++)
			{

				for (short iVar = 0; iVar < nVar; iVar++)
				{

					du = D_du[iCell][jCell][iVar];
					du0 = D_du0[iCell][jCell][iVar];

					un0 = D_un0[iCell][jCell][iVar];
					un1 = D_un[iCell][jCell][iVar];

					local = 0.75 * un0 + 0.25 * (un1 - DT * du);

					D_un1[iCell][jCell][iVar] = un1;
					D_un[iCell][jCell][iVar] = local;
				}
			}
		}
	}
	else
	{
		// #pragma omp parallel for num_threads(NumberThreads)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{
			mydouble local;
			for (long iCell = nb; iCell < nCell_L; iCell++)
			{
				for (short iVar = 0; iVar < nVar; iVar++)
				{
					du = D_du[iCell][jCell][iVar];
					du0 = D_du0[iCell][jCell][iVar];

					un0 = D_un0[iCell][jCell][iVar];
					un1 = D_un1[iCell][jCell][iVar];
					un2 = D_un[iCell][jCell][iVar];
					local = 1.0 / 3.0 * un0 + 2.0 / 3.0 * (un2 - DT * du);
					D_un[iCell][jCell][iVar] = local;
				}
			}
		}
	}
}

void CADSolver::TimeIntegration_TD23_ALW(unsigned short iRK_Step)
{

	// // #pragma omp barrier

	// #pragma omp parallel for private(a21, v2, v1, vv1, vv2, aa21)
	mydouble du, ddu, du0, ddu0, uu, uu0, local;
	mydouble a21, v2, v1, vv1, vv2, aa21;
	// a21 = 0.594223212099088, v2 = 0.306027487008159;
	a21 = 2. / 3., v2 = 3. / 8., aa21 = a21 * a21 * 0.5;
	v1 = 1. - v2, vv1 = 1. / 6. * (3. - 1. / a21 - 3. * a21 * v2), vv2 = 1. / (6. * a21) - (a21 * v2) / 2.;
	// #pragma omp parallel for num_threads(NumberThreads)
	for (long jCell = nb; jCell < nCellYo_L; jCell++)
	{
		for (long iCell = nb; iCell < nCell_L; iCell++)
		{
			for (short iVar = 0; iVar < nVar; iVar++)
			{
				if (iRK_Step == 0)
				{
					du = D_du[iCell][jCell][iVar];
					ddu = D_ddu[iCell][jCell][iVar];
					uu0 = D_un[iCell][jCell][iVar];

					uu = uu0 - a21 * DT * du - aa21 * DT * DT * ddu;

					D_un[iCell][jCell][iVar] = uu;
					D_un0[iCell][jCell][iVar] = uu0;
					D_du0[iCell][jCell][iVar] = du;
					D_ddu0[iCell][jCell][iVar] = ddu;
				}
				else
				{
					du = D_du[iCell][jCell][iVar];
					ddu = D_ddu[iCell][jCell][iVar];
					du0 = D_du0[iCell][jCell][iVar];
					ddu0 = D_ddu0[iCell][jCell][iVar];
					uu0 = D_un0[iCell][jCell][iVar];

					local = uu0 - DT * (v1 * du0 + v2 * du) -
							DT * DT * (vv1 * ddu0 + vv2 * ddu);
					D_un[iCell][jCell][iVar] = local;
				}
			}
		}
	}

	// // #pragma omp barrier
}

void CADSolver::TimeIntegration_TD23_ALW_Hy(unsigned short iRK_Step)
{

	mydouble local;
	mydouble du, ddu, du0, ddu0, uu, uu0, Tu1 = 0., Tu2 = 0., TL_n1 = 0., TL_n2 = 0., TLt_n1 = 0., TLt_n2 = 0.;
	mydouble a21, v2, ww1, ww2, w11, w22, v1, vv1, vv2;
	mydouble h_b1, h_b2, h_v1, h_v2, h_vv1, h_vv2;
	a21 = 0.594223212099088, v2 = 0.306027487008159;
	v1 = 1. - v2, vv1 = 1. / 6. * (3. - 1. / a21 - 3. * a21 * v2), vv2 = 1. / (6. * a21) - (a21 * v2) / 2.;

	h_b1 = 0.2013 * pow(DT_w1, -1.592); // h_b1 = 0.;
	h_b2 = 1.0 - h_b1;
	h_vv1 = 0.0;
	h_v1 = 1. / 3. * (h_b1 + 1. / pow(DT_w1, 3));
	h_v2 = 1. - 1. / (3.0 * DT_w1 * DT_w1) + (2.0 * h_b1 * DT_w1) / 3.0;
	h_vv2 = (2. + 3.0 * DT_w1 - h_b1 * pow(DT_w1, 3)) / (6.0 * DT_w1);

	if (iRK_Step == 0)
	{
		// #pragma omp parallel for private(du, ddu, du0, ddu0, uu, uu0, a21, v2, v1, vv1, vv2, local)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{
			for (long iCell = nb; iCell < nCell_L; iCell++)
			{
				for (short iVar = 0; iVar < nVar; iVar++)
				{
					du = D_du[iCell][jCell][iVar];
					ddu = D_ddu[iCell][jCell][iVar];
					uu0 = D_un[iCell][jCell][iVar];

					uu = uu0 - a21 * DT * du - a21 * a21 * 0.5 * DT * DT * ddu;

					D_un[iCell][jCell][iVar] = uu;
					D_un0[iCell][jCell][iVar] = uu0;
					D_du0[iCell][jCell][iVar] = du;
					D_ddu0[iCell][jCell][iVar] = ddu;
				}
			}
		}
	}
	else if (iRK_Step == 1)
	{
		// #pragma omp parallel for private(du, ddu, du0, ddu0, uu, uu0, a21, v2, v1, vv1, vv2, local)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{

			for (long iCell = nb; iCell < nCell_L; iCell++)
			{

				for (short iVar = 0; iVar < nVar; iVar++)
				{
					uu0 = D_un0[iCell][jCell][iVar];
					du0 = D_du0[iCell][jCell][iVar];
					ddu0 = D_ddu0[iCell][jCell][iVar];
					if (flagh_hy[iCell][jCell][4] > 0.5)
					{
						du = D_du[iCell][jCell][iVar];
						ddu = D_ddu[iCell][jCell][iVar];

						local = uu0 - DT * (v1 * du0 + v2 * du) -
								DT * DT * (vv1 * ddu0 + vv2 * ddu);
					}
					else
					{

						Tu1 = D_un1[iCell][jCell][iVar];
						TL_n1 = D_dun1[iCell][jCell][iVar];
						TLt_n1 = D_ddun1[iCell][jCell][iVar];

						local = h_b1 * Tu1 + h_b2 * uu0 -
								DT * (h_v1 * DT_w1 * TL_n1 + h_v2 * du0) -
								DT * DT * (h_vv1 * DT_w1 * DT_w1 * TLt_n1 + h_vv2 * ddu0);
					}
					D_un[iCell][jCell][iVar] = local;
					D_un1[iCell][jCell][iVar] = uu0;
					D_dun1[iCell][jCell][iVar] = du0;
					D_ddun1[iCell][jCell][iVar] = ddu0;
				}
			}
		}
	}
}

void CADSolver::TimeIntegration_TDMSRK223_ALW(unsigned short iRK_Step)
{

	mydouble local, b1, b2, aa21;
	mydouble du, ddu, du0, ddu0, uu, uu0, Tu1 = 0., Tu2 = 0., TL_n1 = 0., TL_n2 = 0., TLt_n1 = 0., TLt_n2 = 0.;
	mydouble a21, aa22, a22, ww1, ww2, w11, w22, v1, v2, v3, vv1, vv2, vv3, d31, d32, w1;
	w1 = DT_w1;
	a22 = 0.5416 * pow(w1, -0.04681);
	b1 = 0.1397 * exp(0.8707 * w1) - 0.2405 * pow(w1, 1.319);
	v1 = 0.6519 * exp(-2.359 * w1) + 0.01983 * pow(w1, 1.807);
	v3 = 0.3203 * exp(0.4813 * w1) * pow(w1, -0.2711);

	d31 = 0., vv1 = 0., aa21 = 0., a21 = 0.;
	b2 = 1. - b1, d32 = 1. - d31;
	v2 = 1. - v3 + b1 * w1 - v1 * w1;
	aa22 = a22 * a22 * 0.5;
	vv3 = (1. - 3. * a22 * a22 * v3 + (b1 - 3. * v1) * w1 * w1 * w1) / (6. * a22);
	vv2 = -((1. + 3. * a22 * (-1. + a22 * v3) + 3. * a22 * (b1 - 2. * v1) * w1 * w1 + (b1 - 3. * v1) * w1 * w1 * w1) / (6. * a22));

	if (iRK_Step == 0)
	{
		// #pragma omp parallel for private(du, ddu, du0, ddu0, uu, uu0, a21, v2, v1, vv1, vv2, local)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{
			for (long iCell = nb; iCell < nCell_L; iCell++)
			{
				for (short iVar = 0; iVar < nVar; iVar++)
				{
					du = D_du[iCell][jCell][iVar];
					ddu = D_ddu[iCell][jCell][iVar];
					uu0 = D_un[iCell][jCell][iVar];

					Tu1 = D_un1[iCell][jCell][iVar];
					TL_n1 = D_dun1[iCell][jCell][iVar];
					TLt_n1 = D_ddun1[iCell][jCell][iVar];

					local = d31 * Tu1 + d32 * uu0 - DT * (a21 * w1 * TL_n1 + a22 * du) -
							DT * DT * (aa21 * w1 * w1 * TLt_n1 + aa22 * ddu);

					D_un[iCell][jCell][iVar] = local;
					D_un0[iCell][jCell][iVar] = uu0;
					D_du0[iCell][jCell][iVar] = du;
					D_ddu0[iCell][jCell][iVar] = ddu;
				}
			}
		}
	}
	else if (iRK_Step == 1)
	{
		// #pragma omp parallel for private(du, ddu, du0, ddu0, uu, uu0, a21, v2, v1, vv1, vv2, local)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{

			for (long iCell = nb; iCell < nCell_L; iCell++)
			{

				for (short iVar = 0; iVar < nVar; iVar++)
				{
					du = D_du[iCell][jCell][iVar];
					ddu = D_ddu[iCell][jCell][iVar];

					uu0 = D_un0[iCell][jCell][iVar];
					du0 = D_du0[iCell][jCell][iVar];
					ddu0 = D_ddu0[iCell][jCell][iVar];

					Tu1 = D_un1[iCell][jCell][iVar];
					TL_n1 = D_dun1[iCell][jCell][iVar];
					TLt_n1 = D_ddun1[iCell][jCell][iVar];

					local = b1 * Tu1 + b2 * uu0 - DT * (v1 * w1 * TL_n1 + v2 * du0 + v3 * du) -
							DT * DT * (vv1 * w1 * w1 * TLt_n1 + vv2 * ddu0 + vv3 * ddu);

					D_un[iCell][jCell][iVar] = local;

					D_un1[iCell][jCell][iVar] = uu0;
					D_dun1[iCell][jCell][iVar] = du0;
					D_ddun1[iCell][jCell][iVar] = ddu0;
				}
			}
		}
	}
}

void CADSolver::TimeIntegration_TDMS23_ALW(unsigned short iRK_Step)
{

	mydouble du, ddu, uu0, Tu1, TL_n1, TLt_n1, local;
	mydouble h_b1, h_b2, h_v1, h_v2, h_vv1, h_vv2;
	// mydouble a21, v2, v1, vv1, vv2, b0, b1, b2, v0, vv0;
	// v0 = 0.8, b1 = -0.2;
	// v1 = 1 - b1 - v0;
	// vv0 = (7 - b1 - 3 * v0) / 6.;
	// vv1 = (2 - 2 * b1 - 3 * v0) / 6.;
	// b0 = -1 - b1;
	// vv2 = 0.0, v2 = 0., b2 = 0.;

	h_b1 = 0.4091 * pow(DT_w1, -1.861); // h_b1 = 0.;  //Tylor
	h_b2 = 1.0 - h_b1;
	h_vv1 = 0.0;
	h_v1 = 1. / 3. * (h_b1 + 1. / pow(DT_w1, 3));
	h_v2 = 1. - 1. / (3.0 * DT_w1 * DT_w1) + (2.0 * h_b1 * DT_w1) / 3.0;
	h_vv2 = (2. + 3.0 * DT_w1 - h_b1 * pow(DT_w1, 3)) / (6.0 * DT_w1);

	// #pragma omp parallel for private(du, ddu, uu0, Tu1, TL_n1, TLt_n1, local, h_b1, h_b2, h_v1, h_v2, h_vv1, h_vv2)
	for (long jCell = nb; jCell < nCellYo_L; jCell++)
	{
		for (long iCell = nb; iCell < nCell_L; iCell++)
		{

			for (short iVar = 0; iVar < nVar; iVar++)
			{

				uu0 = D_un[iCell][jCell][iVar];
				du = D_du[iCell][jCell][iVar];
				ddu = D_ddu[iCell][jCell][iVar];

				Tu1 = D_un0[iCell][jCell][iVar];
				TL_n1 = D_du0[iCell][jCell][iVar];
				TLt_n1 = D_ddu0[iCell][jCell][iVar];

				local = h_b1 * Tu1 + h_b2 * uu0 -
						DT * (h_v1 * DT_w1 * TL_n1 + h_v2 * du) -
						DT * DT * (h_vv1 * DT_w1 * DT_w1 * TLt_n1 + h_vv2 * ddu);

				D_un[iCell][jCell][iVar] = local;
				D_un0[iCell][jCell][iVar] = uu0;
				D_du0[iCell][jCell][iVar] = du;
				D_ddu0[iCell][jCell][iVar] = ddu;
			}
		}
	}
}

void CADSolver::TimeIntegration_TDMS34_ALW(unsigned short iRK_Step) // Tylor SSP
{

	mydouble local;

	mydouble a21, v2, v1, vv1, vv2, b0, b1, b2, b3, v0, vv0, uu, vv, w1, w2, v3, vv3, w11, w22, ww1, ww2;

	w2 = DT_w1, w1 = DT_w2;
	if (w1 > 1.3)
		w1 = 1.3;
	if (w2 > 1.3)
		w2 = 1.3;
	b1 = -0.283 + 0.088 * w1 * w1 + 6.246 * exp(-1.227 * w1 - 1.354 * w2);
	b2 = 0.205 + 0.3156 * pow(w1, -1.508) * pow(w2, 1.065) - 3.1 * exp(-2.123 * w1);

	v1 = 22.23 + 8.837 * w1 * w1 - 14.05 * exp(0.7255 * w1 - 0.1072 * w2) - 9.275 * sin(0.5058 * w2);

	b3 = 1. - b1 - b2;
	vv2 = 0.;
	uu = w1 + w2, vv = (b1 + b2) * w2;
	vv1 = (3. + 12. * pow(uu, 2) * v1 * pow(w1, 2) - b1 * pow(uu, 3) * (3. * w1 - w2) + 4. * w2 + b2 * pow(w2, 4)) / (12. * uu * pow(w1, 2) * (3. * w1 + w2));

	v2 = (1. + w2 + w1 * (1. + pow(uu, 3) * v1 - 6. * pow(uu, 2) * vv1 * w1 + b2 * pow(w2, 3))) / ((3. * w1 - w2) * pow(w2, 3));

	vv3 = (3. + (b1 - 12. * vv1) * pow(w1, 4) + 2. * (b1 - 6. * vv1) * pow(w1, 3) * w2 + 2. * w1 * (2. + 3. * w2 - vv * pow(w2, 2)) + w2 * (8. + 6. * w2 - vv * pow(w2, 2))) / (12. * uu * w2);

	v3 = (-1. - 2. * uu - (b1 - 2. * v1) * pow(w1, 4) + 2. * (3. + 2. * vv) * w1 * pow(w2, 2) + (2. + vv) * pow(w2, 3) + pow(w1, 3) * (-4. * b1 * w2 + 6. * v1 * w2)) / (2. * pow(w2, 2) * (3. * w1 + w2));
	w2 = DT_w1, w1 = DT_w2;
	// #pragma omp parallel for private(vv1, vv2, vv3, b1, b2, b3, v2, v1, v3, w1, w2, local)
	for (long jCell = nb; jCell < nCellYo_L; jCell++)
	{
		mydouble du, ddu, du0, ddu0, uu0, Tu1 = 0., Tu2 = 0., TL_n1 = 0., TL_n2 = 0., TLt_n1 = 0., TLt_n2 = 0.;
		for (long iCell = nb; iCell < nCell_L; iCell++)
		{

			for (short iVar = 0; iVar < nVar; iVar++)
			{
				uu0 = D_un[iCell][jCell][iVar];
				du = D_du[iCell][jCell][iVar];
				ddu = D_ddu[iCell][jCell][iVar];

				Tu1 = D_un1[iCell][jCell][iVar];
				TL_n1 = D_dun1[iCell][jCell][iVar];
				TLt_n1 = D_ddun1[iCell][jCell][iVar];

				Tu2 = D_un2[iCell][jCell][iVar];
				TL_n2 = D_dun2[iCell][jCell][iVar];
				TLt_n2 = D_ddun2[iCell][jCell][iVar];

				local = b1 * Tu2 + b2 * Tu1 + b3 * uu0 -
						DT * (v1 * w1 * TL_n2 + v2 * w2 * TL_n1 + v3 * du) -
						DT * DT * (vv1 * w1 * w1 * TLt_n2 + vv2 * w2 * w2 * TLt_n1 + vv3 * ddu);

				D_un[iCell][jCell][iVar] = local;

				D_un2[iCell][jCell][iVar] = D_un1[iCell][jCell][iVar];
				D_dun2[iCell][jCell][iVar] = D_dun1[iCell][jCell][iVar];
				D_ddun2[iCell][jCell][iVar] = D_ddun1[iCell][jCell][iVar];

				D_un1[iCell][jCell][iVar] = uu0;
				D_dun1[iCell][jCell][iVar] = du;
				D_ddun1[iCell][jCell][iVar] = ddu;
			}
		}
	}
}

void CADSolver::TimeIntegration_TDMS34_ALW_1(unsigned short iRK_Step) // SD SSP
{

	mydouble local;
	mydouble du, ddu, du0, ddu0, uu, uu0, Tu1 = 0., Tu2 = 0., TL_n1 = 0., TL_n2 = 0., TLt_n1 = 0., TLt_n2 = 0.;
	mydouble a21, v2, v1, vv1, vv2, b0, b1, b2, v0, vv0;
	mydouble h_b1, h_b2, h_v1, h_v2, h_vv1, h_vv2, h_b3, h_v3, h_vv3, h_w1, h_w2;

	h_w2 = DT_w1, h_w1 = DT_w2;
	// if (h_w1 > 1.3)
	// 	h_w1 = 1.3;
	// if (h_w2 > 1.3)
	// 	h_w2 = 1.3;
	h_b1 = 22.27 - 22.26 * pow(h_w1, 0.00209) * exp(-0.0006788 * h_w1) * pow(h_w2, -0.000214);
	h_b2 = 0.06399 + 0.2848 * pow(h_w1, -0.2672) * pow(h_w2, -1.599);

	h_b3 = 1. - h_b1 - h_b2, h_v1 = 0., h_vv2 = 0.;

	h_vv1 = (3. + 4. * h_w2 + h_b2 * pow(h_w2, 4) - h_b1 * (3. * h_w1 - h_w2) * pow(h_w1 + h_w2, 3)) /
			(12. * h_w1 * h_w1 * (h_w1 + h_w2) * (3. * h_w1 + h_w2));
	h_v3 = (-2 - 3 * h_w2 + pow(h_w2, 3) + 6 * h_vv1 * h_w1 * h_w1 * (h_w1 + h_w2) * (4 * h_w1 + h_w2) + h_b1 * h_w1 * h_w1 * (h_w1 + h_w2) * (2 * h_w1 + 3 * h_w2)) / pow(h_w2, 3);
	h_v2 = (-1 - h_w2 + h_b1 * h_w1 * pow(h_w1 + h_w2, 3) + 6 * h_vv1 * h_w1 * h_w1 * (h_w1 + h_w2) * (2 * h_w1 + h_w2)) / pow(h_w2, 4);
	h_vv3 = (-h_vv1) * h_w1 * h_w1 + h_v2 * h_w2 * h_w2 + 0.5 * (1. - h_b1 * (h_w1 + h_w2) * (h_w1 + h_w2) - h_b2 * h_w2 * h_w2);

	h_w2 = DT_w1, h_w1 = DT_w2;

	// #pragma omp parallel for private(du, ddu, du0, ddu0, uu, uu0, v2, v1, vv1, vv2, b0, b1, b2, v0, vv0, local)
	for (long jCell = nb; jCell < nCellYo_L; jCell++)
	{
		for (long iCell = nb; iCell < nCell_L; iCell++)
		{

			for (short iVar = 0; iVar < nVar; iVar++)
			{
				uu0 = D_un[iCell][jCell][iVar];
				du = D_du[iCell][jCell][iVar];
				ddu = D_ddu[iCell][jCell][iVar];

				Tu1 = D_un1[iCell][jCell][iVar];
				TL_n1 = D_dun1[iCell][jCell][iVar];
				TLt_n1 = D_ddun1[iCell][jCell][iVar];

				Tu2 = D_un2[iCell][jCell][iVar];
				TL_n2 = D_dun2[iCell][jCell][iVar];
				TLt_n2 = D_ddun2[iCell][jCell][iVar];

				local = h_b1 * Tu2 + h_b2 * Tu1 + h_b3 * uu0 -
						DT * (h_v1 * h_w1 * TL_n2 + h_v2 * h_w2 * TL_n1 + h_v3 * du) -
						DT * DT * (h_vv1 * h_w1 * h_w1 * TLt_n2 + h_vv2 * h_w2 * h_w2 * TLt_n1 + h_vv3 * ddu);

				D_un[iCell][jCell][iVar] = local;

				D_un2[iCell][jCell][iVar] = D_un1[iCell][jCell][iVar];
				D_dun2[iCell][jCell][iVar] = D_dun1[iCell][jCell][iVar];
				D_ddun2[iCell][jCell][iVar] = D_ddun1[iCell][jCell][iVar];

				D_un1[iCell][jCell][iVar] = uu0;
				D_dun1[iCell][jCell][iVar] = du;
				D_ddun1[iCell][jCell][iVar] = ddu;
			}
		}
	}
}

void CADSolver::TimeIntegration_TD24_ALW(unsigned short iRK_Step)
{

	mydouble local;
	mydouble du, ddu, du0, ddu0, uu, uu0;
	mydouble a21, v2, ww1, ww2, w11, w22, v1, vv1, vv2;

	a21 = 0.5, v1 = 1., v2 = 0., vv1 = 1. / 6., vv2 = 1. / 3.; // C0 = 1. / 3., C1 = 2. / 3.;
	if (iRK_Step == 0)
	{
		// #pragma omp parallel for private(du, ddu, du0, ddu0, uu, uu0, a21, v2, v1, vv1, vv2, local)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{
			for (long iCell = nb; iCell < nCell_L; iCell++)
			{
				for (short iVar = 0; iVar < nVar; iVar++)
				{
					du = D_du[iCell][jCell][iVar];
					ddu = D_ddu[iCell][jCell][iVar];
					uu0 = D_un[iCell][jCell][iVar];
					// ddu = 0.;
					uu = uu0 - a21 * DT * du - a21 * a21 * 0.5 * DT * DT * ddu;

					// uu = uu0 - DT * du - 0.5 * DT * DT * ddu;
					// uu = uu0- DT * du;

					D_un[iCell][jCell][iVar] = uu;
					D_un0[iCell][jCell][iVar] = uu0;
					D_du0[iCell][jCell][iVar] = du;
					D_ddu0[iCell][jCell][iVar] = ddu;
				}
			}
		}
	}
	else if (iRK_Step == 1)
	{
		// #pragma omp parallel for private(du, ddu, du0, ddu0, uu, uu0, a21, v2, v1, vv1, vv2, local)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{

			for (long iCell = nb; iCell < nCell_L; iCell++)
			{

				for (short iVar = 0; iVar < nVar; iVar++)
				{
					du = D_du[iCell][jCell][iVar];
					ddu = D_ddu[iCell][jCell][iVar];

					du0 = D_du0[iCell][jCell][iVar];
					ddu0 = D_ddu0[iCell][jCell][iVar];

					uu0 = D_un0[iCell][jCell][iVar];

					local = uu0 - DT * (v1 * du0 + v2 * du) -
							DT * DT * (vv1 * ddu0 + vv2 * ddu);
					// local = 0.75 * Element[iCell][jCell]->GetSolution_Old(0, iVar) +
					// 		1.0 / 4.0 * (Element[iCell][jCell]->GetSolution(0, 0, iVar) - DT * Element[iCell][jCell]->GetSolution_du(0, iVar));
					D_un[iCell][jCell][iVar] = local;
				}
			}
		}
	}
}

void CADSolver::TimeIntegration_TD34_ALW(unsigned short iRK_Step)
{

	mydouble local;
	mydouble du, ddu, du0, ddu0, uu, uu0, du1, ddu1, du2, ddu2;
	mydouble a21, v2, ww1, ww2, w11, w22, v1, vv1, vv2, a31, a32, aa31, aa32, v3, vv3, aa21;

	a21 = 0.532230958301739;
	a31 = 0.457385979790557;
	a32 = 0.207902718086617;
	aa31 = 0.0553261314403881;
	v3 = 0.374215233358609;

	aa21 = pow(a21, 2) / 2.;
	aa32 = (pow(a31, 2) - 2. * a21 * a32 + 2. * a31 * a32 + pow(a32, 2) - 2. * aa31) / 2.;
	uu = a31 + a32;
	v2 = -(((3. * pow(a21, 2) * a32 + 6. * a21 * aa31 - 3. * a21 * pow(uu, 2) + pow(uu, 3)) * v3) / pow(a21, 3));
	v1 = 1. - v2 - v3;
	vv1 = -(-1. + 2. * uu - 2. * a21 * (-1. + 3. * uu + a21 * (a21 - 3. * uu) * v2) + 2. * (3. * a21 - uu) * pow(uu, 2) * v3) / (12. * a21 * uu);
	vv2 = -(-1. + 2. * uu + 4. * pow(a21, 3) * v2 - 6. * pow(a21, 2) * uu * v2 - 2. * pow(uu, 3) * v3) / (12. * a21 * (a21 - uu));
	vv3 = -(-1. + 2. * a21 - 2. * pow(a21, 3) * v2 - 6. * a21 * pow(uu, 2) * v3 + 4. * pow(uu, 3) * v3) / (12. * uu * (-a21 + uu));

	if (iRK_Step == 0)
	{
		// #pragma omp parallel for private(du, ddu, du0, ddu0, uu, uu0, a21, v2, v1, vv1, vv2, local)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{
			for (long iCell = nb; iCell < nCell_L; iCell++)
			{
				for (short iVar = 0; iVar < nVar; iVar++)
				{
					du = D_du[iCell][jCell][iVar];
					ddu = D_ddu[iCell][jCell][iVar];
					uu0 = D_un[iCell][jCell][iVar];
					// ddu = 0.;
					uu = uu0 - a21 * DT * du - DT * DT * aa21 * ddu;

					// uu = uu0 - DT * du - 0.5 * DT * DT * ddu;
					// uu = uu0- DT * du;

					D_un[iCell][jCell][iVar] = uu;
					D_un0[iCell][jCell][iVar] = uu0;
					D_du0[iCell][jCell][iVar] = du;
					D_ddu0[iCell][jCell][iVar] = ddu;
				}
			}
		}
	}
	else if (iRK_Step == 1)
	{
		// #pragma omp parallel for private(du, ddu, du0, ddu0, uu, uu0, a21, v2, v1, vv1, vv2, local)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{

			for (long iCell = nb; iCell < nCell_L; iCell++)
			{

				for (short iVar = 0; iVar < nVar; iVar++)
				{
					du = D_du[iCell][jCell][iVar];
					ddu = D_ddu[iCell][jCell][iVar];

					du0 = D_du0[iCell][jCell][iVar];
					ddu0 = D_ddu0[iCell][jCell][iVar];

					uu0 = D_un0[iCell][jCell][iVar];

					// local = uu0 - DT * (v1 * du0 + v2 * du) -
					// 		DT * DT * (vv1 * ddu0 + vv2 * ddu);
					local = uu0 - DT * (a31 * du0 + a32 * du) - DT * DT * (aa31 * ddu0 + aa32 * ddu);
					// local = 0.75 * Element[iCell][jCell]->GetSolution_Old(0, iVar) +
					// 		1.0 / 4.0 * (Element[iCell][jCell]->GetSolution(0, 0, iVar) - DT * Element[iCell][jCell]->GetSolution_du(0, iVar));
					D_un[iCell][jCell][iVar] = local;
					D_du1[iCell][jCell][iVar] = du;
					D_ddu1[iCell][jCell][iVar] = ddu;
				}
			}
		}
	}
	else if (iRK_Step == 2)
	{
		// #pragma omp parallel for private(du, ddu, du0, ddu0, uu, uu0, a21, v2, v1, vv1, vv2, local)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{

			for (long iCell = nb; iCell < nCell_L; iCell++)
			{

				for (short iVar = 0; iVar < nVar; iVar++)
				{
					du = D_du[iCell][jCell][iVar];
					ddu = D_ddu[iCell][jCell][iVar];

					du0 = D_du0[iCell][jCell][iVar];
					ddu0 = D_ddu0[iCell][jCell][iVar];

					du1 = D_du1[iCell][jCell][iVar];
					ddu1 = D_ddu1[iCell][jCell][iVar];

					uu0 = D_un0[iCell][jCell][iVar];

					local = uu0 - DT * (v1 * du0 + v2 * du1 + v3 * du) -
							DT * DT * (vv1 * ddu0 + vv2 * ddu1 + vv3 * ddu);
					// local = 0.75 * Element[iCell][jCell]->GetSolution_Old(0, iVar) +
					// 		1.0 / 4.0 * (Element[iCell][jCell]->GetSolution(0, 0, iVar) - DT * Element[iCell][jCell]->GetSolution_du(0, iVar));
					D_un[iCell][jCell][iVar] = local;
				}
			}
		}
	}
}

void CADSolver::TimeIntegration_TDMSRK224_ALW(unsigned short iRK_Step)
{

	mydouble local, b1, b2, aa21;
	mydouble du, ddu, du0, ddu0, uu, uu0, Tu1 = 0., Tu2 = 0., TL_n1 = 0., TL_n2 = 0., TLt_n1 = 0., TLt_n2 = 0.;
	mydouble a21, aa22, a22, ww1, ww2, w11, w22, v1, v2, v3, vv1, vv2, vv3, d31, d32, w1, a21_d31, w1_2, w1_3, w1_4;
	w1 = DT_w1;
	a21 = 0.07982 * pow(w1, -2.688) + 0.01096;
	b1 = 0.06514 * pow(w1, -1.418) + 0.0125;
	d31 = 0.06017 * pow(w1, -1.835) + 0.02575;
	v3 = 0.3985 * exp(-0.2794 * w1) * pow(w1, 0.2483);
	vv3 = 0.1651 * exp(0.1884 * w1) * pow(w1, -0.2552);

	b2 = 1. - b1, d32 = 1. - d31, aa21 = 0.0;
	a21_d31 = pow(3. * a21 - d31, 1. / 3.) * w1;
	w1_2 = w1 * w1, w1_3 = w1_2 * w1, w1_4 = w1_3 * w1;

	a22 = w1 * (d31 - a21) + a21_d31;
	aa22 = 0.5 * ((2. * a21 - d31) * w1 * w1 + a21_d31 * a21_d31);
	v1 = (1. + 2. * w1 + b1 * w1_4 - 2. * a21_d31 * (6. * vv3 * w1 + a21_d31 * (6. * vv3 + 3. * v3 * w1 + 2. * v3 * a21_d31))) / (2. * w1_4);
	v2 = 1. - v3 + (b1 - v1) * w1;
	vv1 = v1 * 0.5 - (1. + b1 * w1_3 - 6. * vv3 * a21_d31 - 3. * v3 * a21_d31 * a21_d31) / (6. * w1_3);
	vv2 = vv1 * w1_2 + (2. + (3. - 6. * vv3) * w1 - b1 * w1_3 - 6. * a21_d31 * (2. * vv3 + v3 * w1 + v3 * a21_d31)) / (6. * w1);

	if (iRK_Step == 0)
	{
		// #pragma omp parallel for private(du, ddu, du0, ddu0, uu, uu0, a21, v2, v1, vv1, vv2, local)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{
			for (long iCell = nb; iCell < nCell_L; iCell++)
			{
				for (short iVar = 0; iVar < nVar; iVar++)
				{
					du = D_du[iCell][jCell][iVar];
					ddu = D_ddu[iCell][jCell][iVar];
					uu0 = D_un[iCell][jCell][iVar];

					Tu1 = D_un1[iCell][jCell][iVar];
					TL_n1 = D_dun1[iCell][jCell][iVar];
					TLt_n1 = D_ddun1[iCell][jCell][iVar];

					local = d31 * Tu1 + d32 * uu0 - DT * (a21 * w1 * TL_n1 + a22 * du) -
							DT * DT * (aa21 * w1 * w1 * TLt_n1 + aa22 * ddu);

					D_un[iCell][jCell][iVar] = local;
					D_un0[iCell][jCell][iVar] = uu0;
					D_du0[iCell][jCell][iVar] = du;
					D_ddu0[iCell][jCell][iVar] = ddu;
				}
			}
		}
	}
	else if (iRK_Step == 1)
	{
		// #pragma omp parallel for private(du, ddu, du0, ddu0, uu, uu0, a21, v2, v1, vv1, vv2, local)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{

			for (long iCell = nb; iCell < nCell_L; iCell++)
			{

				for (short iVar = 0; iVar < nVar; iVar++)
				{
					du = D_du[iCell][jCell][iVar];
					ddu = D_ddu[iCell][jCell][iVar];

					uu0 = D_un0[iCell][jCell][iVar];
					du0 = D_du0[iCell][jCell][iVar];
					ddu0 = D_ddu0[iCell][jCell][iVar];

					Tu1 = D_un1[iCell][jCell][iVar];
					TL_n1 = D_dun1[iCell][jCell][iVar];
					TLt_n1 = D_ddun1[iCell][jCell][iVar];

					local = b1 * Tu1 + b2 * uu0 - DT * (v1 * w1 * TL_n1 + v2 * du0 + v3 * du) -
							DT * DT * (vv1 * w1 * w1 * TLt_n1 + vv2 * ddu0 + vv3 * ddu);

					D_un[iCell][jCell][iVar] = local;

					D_un1[iCell][jCell][iVar] = uu0;
					D_dun1[iCell][jCell][iVar] = du0;
					D_ddun1[iCell][jCell][iVar] = ddu0;
				}
			}
		}
	}
}

void CADSolver::TimeIntegration_TDMSRK224_ALW_Hy(unsigned short iRK_Step)
{

	mydouble local, b1, b2, aa21;
	mydouble du, ddu, du0, ddu0, uu, uu0, Tu1 = 0., Tu2 = 0., TL_n1 = 0., TL_n2 = 0., TLt_n1 = 0., TLt_n2 = 0.;
	mydouble a21, aa22, a22, ww1, ww2, w11, w22, v1, v2, v3, vv1, vv2, vv3, d31, d32, w1, a21_d31, w1_2, w1_3, w1_4;
	mydouble h_b1, h_b2, h_v1, h_v2, h_vv1, h_vv2, h_b3, h_v3, h_vv3, h_w1, h_w2;

	h_w2 = DT_w1, h_w1 = DT_w2;

	h_b1 = 22.27 - 22.26 * pow(h_w1, 0.00209) * exp(-0.0006788 * h_w1) * pow(h_w2, -0.000214);
	h_b2 = 0.06399 + 0.2848 * pow(h_w1, -0.2672) * pow(h_w2, -1.599);

	h_b3 = 1. - h_b1 - h_b2, h_v1 = 0., h_vv2 = 0.;

	h_vv1 = (3. + 4. * h_w2 + h_b2 * pow(h_w2, 4) - h_b1 * (3. * h_w1 - h_w2) * pow(h_w1 + h_w2, 3)) /
			(12. * h_w1 * h_w1 * (h_w1 + h_w2) * (3. * h_w1 + h_w2));
	h_v3 = (-2 - 3 * h_w2 + pow(h_w2, 3) + 6 * h_vv1 * h_w1 * h_w1 * (h_w1 + h_w2) * (4 * h_w1 + h_w2) + h_b1 * h_w1 * h_w1 * (h_w1 + h_w2) * (2 * h_w1 + 3 * h_w2)) / pow(h_w2, 3);
	h_v2 = (-1 - h_w2 + h_b1 * h_w1 * pow(h_w1 + h_w2, 3) + 6 * h_vv1 * h_w1 * h_w1 * (h_w1 + h_w2) * (2 * h_w1 + h_w2)) / pow(h_w2, 4);
	h_vv3 = (-h_vv1) * h_w1 * h_w1 + h_v2 * h_w2 * h_w2 + 0.5 * (1. - h_b1 * (h_w1 + h_w2) * (h_w1 + h_w2) - h_b2 * h_w2 * h_w2);
	//__________________________________________________

	w1 = DT_w1;
	a21 = 0.07982 * pow(w1, -2.688) + 0.01096;
	b1 = 0.06514 * pow(w1, -1.418) + 0.0125;
	d31 = 0.06017 * pow(w1, -1.835) + 0.02575;
	v3 = 0.3985 * exp(-0.2794 * w1) * pow(w1, 0.2483);
	vv3 = 0.1651 * exp(0.1884 * w1) * pow(w1, -0.2552);

	b2 = 1. - b1, d32 = 1. - d31, aa21 = 0.0;
	a21_d31 = pow(3. * a21 - d31, 1. / 3.) * w1;
	w1_2 = w1 * w1, w1_3 = w1_2 * w1, w1_4 = w1_3 * w1;

	a22 = w1 * (d31 - a21) + a21_d31;
	aa22 = 0.5 * ((2. * a21 - d31) * w1 * w1 + a21_d31 * a21_d31);
	v1 = (1. + 2. * w1 + b1 * w1_4 - 2. * a21_d31 * (6. * vv3 * w1 + a21_d31 * (6. * vv3 + 3. * v3 * w1 + 2. * v3 * a21_d31))) / (2. * w1_4);
	v2 = 1. - v3 + (b1 - v1) * w1;
	vv1 = v1 * 0.5 - (1. + b1 * w1_3 - 6. * vv3 * a21_d31 - 3. * v3 * a21_d31 * a21_d31) / (6. * w1_3);
	vv2 = vv1 * w1_2 + (2. + (3. - 6. * vv3) * w1 - b1 * w1_3 - 6. * a21_d31 * (2. * vv3 + v3 * w1 + v3 * a21_d31)) / (6. * w1);

	if (iRK_Step == 0)
	{
		// #pragma omp parallel for private(du, ddu, du0, ddu0, uu, uu0, a21, v2, v1, vv1, vv2, local)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{
			for (long iCell = nb; iCell < nCell_L; iCell++)
			{
				for (short iVar = 0; iVar < nVar; iVar++)
				{
					du = D_du[iCell][jCell][iVar];
					ddu = D_ddu[iCell][jCell][iVar];
					uu0 = D_un[iCell][jCell][iVar];

					local = uu0 - a21 * DT * du - a21 * a21 * 0.5 * DT * DT * ddu;

					D_un[iCell][jCell][iVar] = local;
					D_un0[iCell][jCell][iVar] = uu0;
					D_du0[iCell][jCell][iVar] = du;
					D_ddu0[iCell][jCell][iVar] = ddu;
				}
			}
		}
	}
	else if (iRK_Step == 1)
	{
		// #pragma omp parallel for private(du, ddu, du0, ddu0, uu, uu0, a21, v2, v1, vv1, vv2, local)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{

			for (long iCell = nb; iCell < nCell_L; iCell++)
			{

				for (short iVar = 0; iVar < nVar; iVar++)
				{
					uu0 = D_un0[iCell][jCell][iVar];
					du0 = D_du0[iCell][jCell][iVar];
					ddu0 = D_ddu0[iCell][jCell][iVar];
					if (flagh_hy[iCell][jCell][4] > 0.5)
					{

						du = D_du[iCell][jCell][iVar];
						ddu = D_ddu[iCell][jCell][iVar];

						uu0 = D_un0[iCell][jCell][iVar];
						du0 = D_du0[iCell][jCell][iVar];
						ddu0 = D_ddu0[iCell][jCell][iVar];

						Tu1 = D_un1[iCell][jCell][iVar];
						TL_n1 = D_dun1[iCell][jCell][iVar];
						TLt_n1 = D_ddun1[iCell][jCell][iVar];

						local = b1 * Tu1 + b2 * uu0 - DT * (v1 * w1 * TL_n1 + v2 * du0 + v3 * du) -
								DT * DT * (vv1 * w1 * w1 * TLt_n1 + vv2 * ddu0 + vv3 * ddu);
					}
					else
					{

						Tu1 = D_un1[iCell][jCell][iVar];
						TL_n1 = D_dun1[iCell][jCell][iVar];
						TLt_n1 = D_ddun1[iCell][jCell][iVar];

						Tu2 = D_un2[iCell][jCell][iVar];
						TL_n2 = D_dun2[iCell][jCell][iVar];
						TLt_n2 = D_ddun2[iCell][jCell][iVar];

						local = h_b1 * Tu2 + h_b2 * Tu1 + h_b3 * uu0 -
								DT * (h_v1 * h_w1 * TL_n2 + h_v2 * h_w2 * TL_n1 + h_v3 * du0) -
								DT * DT * (h_vv1 * h_w1 * h_w1 * TLt_n2 + h_vv2 * h_w2 * h_w2 * TLt_n1 + h_vv3 * ddu0);
					}
					D_un[iCell][jCell][iVar] = local;

					D_un2[iCell][jCell][iVar] = D_un1[iCell][jCell][iVar];
					D_dun2[iCell][jCell][iVar] = D_dun1[iCell][jCell][iVar];
					D_ddun2[iCell][jCell][iVar] = D_ddun1[iCell][jCell][iVar];

					D_un1[iCell][jCell][iVar] = D_un0[iCell][jCell][iVar];
					D_dun1[iCell][jCell][iVar] = D_du0[iCell][jCell][iVar];
					D_ddun1[iCell][jCell][iVar] = D_ddu0[iCell][jCell][iVar];
				}
			}
		}
	}
}

void CADSolver::TimeIntegration_RK22(unsigned short iRK_Step)
{

	mydouble local;
	mydouble du, ddu, du0, ddu0, uu, uu0, un0, un1;
	mydouble a21, v2, ww1, ww2, w11, w22, v1, vv1, vv2;
	a21 = 0.594223212099088, v2 = 0.306027487008159;
	ww1 = 0.0, ww2 = 0.0, w11 = 0.0, w22 = 0.0;
	v1 = 1. - v2, vv1 = 1. / 6. * (3. - 1. / a21 - 3. * a21 * v2), vv2 = 1. / (6. * a21) - (a21 * v2) / 2.;

	if (iRK_Step == 0)
	{
		// // #pragma omp parallel for num_threads(NumberThreads)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{
			for (long iCell = nb; iCell < nCell_L; iCell++)
			{
				for (short iVar = 0; iVar < nVar; iVar++)
				{

					du = D_du[iCell][jCell][iVar];
					un0 = D_un[iCell][jCell][iVar];

					local = un0 - DT * du;

					D_un[iCell][jCell][iVar] = local;

					D_un0[iCell][jCell][iVar] = un0;
					D_du0[iCell][jCell][iVar] = du;
				}
			}
		}
	}
	else if (iRK_Step == 1)
	{
		// // #pragma omp parallel for num_threads(NumberThreads)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{

			for (long iCell = nb; iCell < nCell_L; iCell++)
			{

				for (short iVar = 0; iVar < nVar; iVar++)
				{

					du = D_du[iCell][jCell][iVar];
					du0 = D_du0[iCell][jCell][iVar];

					un0 = D_un0[iCell][jCell][iVar];
					un1 = D_un[iCell][jCell][iVar];

					// local = Element[iCell][jCell]->GetSolution_Old(0, iVar) -
					// 		DT * (v1 * du0 + v2 * du) -
					// 		DT * DT * (vv1 * ddu0 + vv2 * ddu);
					// local = 0.75 * Element[iCell][jCell]->GetSolution_Old(0, iVar) +
					// 		1.0 / 4.0 * (Element[iCell][jCell]->GetSolution(0, 0, iVar) - DT * Element[iCell][jCell]->GetSolution_du(0, iVar));

					local = 0.5 * (un0 + un1) - 0.5 * DT * du;
					D_un[iCell][jCell][iVar] = local;
				}
			}
		}
	}
}

void CADSolver::TimeIntegration_RK54(unsigned short iRK_Step)
{

	mydouble local;
	mydouble du, du0, un0, un1, un2, un3, un4, du1;

	if (iRK_Step == 0)
	{
		// #pragma omp parallel for num_threads(NumberThreads)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{
			mydouble local;
			for (long iCell = nb; iCell < nCell_L; iCell++)
			{
				for (short iVar = 0; iVar < nVar; iVar++)
				{

					du = D_du[iCell][jCell][iVar];
					un0 = D_un[iCell][jCell][iVar];

					local = un0 - 0.391752226571890 * DT * du;

					D_un[iCell][jCell][iVar] = local;

					D_un0[iCell][jCell][iVar] = un0;
					D_du0[iCell][jCell][iVar] = du;
				}
			}
		}
	}
	else if (iRK_Step == 1)
	{
		// #pragma omp parallel for num_threads(NumberThreads)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{
			mydouble local;
			for (long iCell = nb; iCell < nCell_L; iCell++)
			{

				for (short iVar = 0; iVar < nVar; iVar++)
				{
					un0 = D_un0[iCell][jCell][iVar];
					du0 = D_du0[iCell][jCell][iVar];

					un1 = D_un[iCell][jCell][iVar];
					du = D_du[iCell][jCell][iVar];

					local = 0.444370493651235 * un0 + 0.555629506348765 * un1 - 0.368410593050371 * DT * du;

					D_un1[iCell][jCell][iVar] = un1;
					D_un[iCell][jCell][iVar] = local;
				}
			}
		}
	}
	else if (iRK_Step == 2)
	{
		// #pragma omp parallel for num_threads(NumberThreads)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{
			mydouble local;
			for (long iCell = nb; iCell < nCell_L; iCell++)
			{
				for (short iVar = 0; iVar < nVar; iVar++)
				{
					du = D_du[iCell][jCell][iVar];
					du0 = D_du0[iCell][jCell][iVar];

					un0 = D_un0[iCell][jCell][iVar];
					un1 = D_un1[iCell][jCell][iVar];
					un2 = D_un[iCell][jCell][iVar];
					local = 1.0 / 3.0 * un0 + 2.0 / 3.0 * (un2 - DT * du);

					local = 0.620101851488403 * un0 + 0.379898148511597 * un2 - 0.251891774271694 * DT * du;

					D_un2[iCell][jCell][iVar] = un2;
					D_un[iCell][jCell][iVar] = local;
				}
			}
		}
	}
	else if (iRK_Step == 3)
	{
		// #pragma omp parallel for num_threads(NumberThreads)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{
			mydouble local;
			for (long iCell = nb; iCell < nCell_L; iCell++)
			{
				for (short iVar = 0; iVar < nVar; iVar++)
				{
					du = D_du[iCell][jCell][iVar];
					un0 = D_un0[iCell][jCell][iVar];
					un3 = D_un[iCell][jCell][iVar];

					local = 0.178079954393132 * un0 + 0.821920045606868 * un3 - 0.544974750228521 * DT * du;
					D_un3[iCell][jCell][iVar] = un3;
					D_du1[iCell][jCell][iVar] = du;
					D_un[iCell][jCell][iVar] = local;
				}
			}
		}
	}
	else if (iRK_Step == 4)
	{
		// #pragma omp parallel for num_threads(NumberThreads)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{
			mydouble local;
			for (long iCell = nb; iCell < nCell_L; iCell++)
			{
				for (short iVar = 0; iVar < nVar; iVar++)
				{

					du = D_du[iCell][jCell][iVar];
					un0 = D_un0[iCell][jCell][iVar];

					un4 = D_un[iCell][jCell][iVar];
					un3 = D_un3[iCell][jCell][iVar];
					un2 = D_un2[iCell][jCell][iVar];

					du1 = D_du1[iCell][jCell][iVar];

					local = 0.517231671970585 * un2 + 0.096059710526147 * un3 - 0.063692468666290 * DT * du1 +
							0.386708617503269 * un4 - 0.226007483236906 * DT * du;
					D_un[iCell][jCell][iVar] = local;
				}
			}
		}
	}
}

void CADSolver::store_du_duu()
{

	mydouble local;
	mydouble du, ddu, du0, ddu0, uu, uu0;
	mydouble a21, v2, ww1, ww2, w11, w22, v1, vv1, vv2;

	// a21 = 0.5, v1 = 1., v2 = 0., vv1 = 1. / 6., vv2 = 1. / 3.; // C0 = 1. / 3., C1 = 2. / 3.;

	// #pragma omp parallel for private(du, ddu, du0, ddu0, uu, uu0, a21, v2, v1, vv1, vv2, local)
	for (long jCell = nb; jCell < nCellYo_L; jCell++)
	{
		for (long iCell = nb; iCell < nCell_L; iCell++)
		{
			for (short iVar = 0; iVar < nVar; iVar++)
			{

				D_un2[iCell][jCell][iVar] = D_un1[iCell][jCell][iVar];
				D_dun2[iCell][jCell][iVar] = D_dun1[iCell][jCell][iVar];
				D_ddun2[iCell][jCell][iVar] = D_ddun1[iCell][jCell][iVar];

				D_un1[iCell][jCell][iVar] = D_un0[iCell][jCell][iVar];
				D_dun1[iCell][jCell][iVar] = D_du0[iCell][jCell][iVar];
				D_ddun1[iCell][jCell][iVar] = D_ddu0[iCell][jCell][iVar];
			}
		}
	}
}
