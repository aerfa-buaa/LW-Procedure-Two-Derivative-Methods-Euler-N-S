// #pragma once
#include "xBurgers.hpp"

using namespace std;

void CSolver::Driver(void) {}

void CADSolver::Driver(void)
{

	InputControl_xy();

	MesherAndAllocation_xy();

	// Mesher_out();

	Initial_date();

	cout << " Solver Cell Number X Orientation = " << nCell << endl;
	cout << " Solver Cell Number Y Orientation = " << nCellYo << endl;
	omp_set_num_threads(NumberThreads);

	Boundary_uu();
	uu_to_cc();

	// out_date();
	if (Viscous_flag == 1)
		Viscous_mu();
	auto start = std::chrono::high_resolution_clock::now();
	auto start1 = std::chrono::high_resolution_clock::now();
	ExtIter = 1;
	TotalIter = 0;
	iSave = 0;
	TimeNow = 0.0;
	T_flag = 0.;
	while (ExtIter < StepMax)
	{

		ComputeDT();
		switch (TimeScheme)
		{
		case 11:
			RK11();
			break;
		case 33:
			RK33();
			break;

		case 54:
			RK54();
			break;

		case 22:
			RK22();
			break;
		case 230:
			TD23_old();
			break;
		case 231:
			TD23_New();
			break;

		case 2300:
			TD23();
			break;

		case 240:
			TD24();
			break;

 
		case 341:
			TD34_New();
			break;
		case 340:
			TD34_Old();
			break;

		case 923:
			TDMS23_New();
			break;

		case 9230:
			TDMS23_Old();
			break;

		case 934:
			TDMS34_New();
			break;

		case 9340:
			TDMS34_Old();
			break;

		default:
			cout << "TimeScheme is not avaliable for this Solver." << endl;
			break;
		}
		TimeNow += DT;

		// FinishTime = mydouble(clock()) / mydouble(CLOCKS_PER_SEC);
		// FinishTime = omp_get_wtime();
		// CurrentTotalTime = FinishTime - StartTime;
		if (ExtIter % show_n == 0)
		{
			cout << "Step = " << ExtIter << ";      Time = " << TimeNow << endl;
			auto stop1 = std::chrono::high_resolution_clock::now();
			auto duration1 = std::chrono::duration_cast<std::chrono::milliseconds>(stop1 - start1);
			std::cout << "Elapsed time: " << duration1.count() << " milliseconds" << std::endl;
		}
		ExtIter = ExtIter + 1;
		if (T_flag > 2.0)
			break;

		if (ska > 5)
		{
			cout << "Step = " << ExtIter << ";      Time = " << TimeNow << endl;
			cout << "NAN__NAN__NAN__NAN__NAN__NAN__= " << endl;
			break;
		}
	}
	cout << "Step = " << ExtIter << ";      Time = " << TimeNow << endl;
	auto stop_n = std::chrono::high_resolution_clock::now();
	// Calculate the elapsed time
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop_n - start);
	// Output the elapsed time in milliseconds
	std::cout << "Elapsed time: " << duration.count() << " milliseconds" << std::endl;
	cout << "TimeScheme = " << TimeScheme << ";  CFL=" << CFL << endl;
	// #pragma omp barrier

	comput_err1();
	out_date();
	Delete_date();
}

void CADSolver::RK33(void)
{

	for (iRK_Step = 0; iRK_Step < 3; iRK_Step++)
	{

		Non_viscous_Splitting();
		if (Viscous_flag == 1)
			Viscous_Splitting();
		TimeIntegration(iRK_Step);
		Boundary_uu();
		uu_to_cc();
	}
}

void CADSolver::RK54(void)
{

	for (iRK_Step = 0; iRK_Step < 5; iRK_Step++)
	{

		Non_viscous_Splitting();
		if (Viscous_flag == 1)
			Viscous_Splitting();

		TimeIntegration_RK54(iRK_Step);

		Boundary_uu();
		uu_to_cc();
	}
}

void CADSolver::RK11(void)
{

	for (iRK_Step = 0; iRK_Step < 1; iRK_Step++)
	{

		Non_viscous_Splitting();
		TimeIntegration_RK11(iRK_Step);

		Boundary_uu();
		uu_to_cc();
	}
}

void CADSolver::RK22(void)
{

	for (iRK_Step = 0; iRK_Step < 2; iRK_Step++)
	{

		Non_viscous_Splitting();
		TimeIntegration_RK22(iRK_Step);

		Boundary_uu();
		uu_to_cc();
	}
}

void CADSolver::TD23(void)
{

	for (iRK_Step = 0; iRK_Step < 2; iRK_Step++)
	{

		Non_viscous_Splitting();

		comput_du_ALW();
		scheme_du_dx_New();
		TimeIntegration_TD23_ALW(iRK_Step);

		Boundary_uu();
		uu_to_cc();
	}
}

 
void CADSolver::TD23_old(void)
{

	for (iRK_Step = 0; iRK_Step < 2; iRK_Step++)
	{

		Non_viscous_Splitting();
 
		comput_du_LW_Old();
		scheme_du_dx_New();
 
		TimeIntegration_TD23_ALW(iRK_Step);

		Boundary_uu();
		uu_to_cc();
	}
}

void CADSolver::TD23_New(void)
{

	for (iRK_Step = 0; iRK_Step < 2; iRK_Step++)
	{

		Non_viscous_Splitting();
		if (Viscous_flag == 1)
			Viscous_Splitting();

		comput_du_LW_New();
		scheme_du_dx_New();
		if (Viscous_flag == 1)
			Viscous_Splitting_LW();

		TimeIntegration_TD23_ALW(iRK_Step);

		Boundary_uu();
		uu_to_cc();
	}
}

void CADSolver::TD34_New(void)
{

	for (iRK_Step = 0; iRK_Step < 3; iRK_Step++)
	{
		Non_viscous_Splitting();
		if (Viscous_flag == 1)
			Viscous_Splitting();

		comput_du_LW_New();
		scheme_du_dx_New();
		if (Viscous_flag == 1)
			Viscous_Splitting_LW();

		TimeIntegration_TD34_ALW(iRK_Step);

		Boundary_uu();
		uu_to_cc();
	}
}

void CADSolver::TD34_Old(void)
{

	for (iRK_Step = 0; iRK_Step < 3; iRK_Step++)
	{

		Non_viscous_Splitting();

		comput_du_LW_Old();
		scheme_du_dx_New();
		TimeIntegration_TD34_ALW(iRK_Step);

		Boundary_uu();
		uu_to_cc();
	}
}

void CADSolver::TD24(void)
{

	for (iRK_Step = 0; iRK_Step < 2; iRK_Step++)
	{

		Non_viscous_Splitting();

		comput_du_ALW();
		scheme_du_dx_New();
		TimeIntegration_TD24_ALW(iRK_Step);

		Boundary_uu();
		uu_to_cc();
	}
}

void CADSolver::TD24_Old(void)
{

	for (iRK_Step = 0; iRK_Step < 2; iRK_Step++)
	{

		Non_viscous_Splitting();

		// comput_du_ALW();
		comput_du_LW_Old();
		scheme_du_dx_New();
		TimeIntegration_TD24_ALW(iRK_Step);

		Boundary_uu();
		uu_to_cc();
	}
}

void CADSolver::TDMS23_Old(void)
{

	if (ExtIter < 4)
	{
		TD23_old();
		store_du_duu();
	}
	else if (T_flag < 2.) //
	{

		Non_viscous_Splitting();

		// comput_du_LW_New();
		comput_du_LW_Old();
		scheme_du_dx_New();
		TimeIntegration_TDMS23_ALW(iRK_Step);

		Boundary_uu();
		uu_to_cc();
	}
	else if (T_flag > 3.)
	{
		TD23_old();
	}
}

void CADSolver::TDMS23_New(void)
{

	if (ExtIter < 4)
	{
		TD23_New();
		store_du_duu();
	}
	else if (T_flag < 2.) // else if (DT_w1 < 1.5 && DT_w1 > 0.75) //
	{

		Non_viscous_Splitting();
		if (Viscous_flag == 1)
			Viscous_Splitting();

		comput_du_LW_New();
		scheme_du_dx_New();
		if (Viscous_flag == 1)
			Viscous_Splitting_LW();

		TimeIntegration_TDMS23_ALW(iRK_Step);

		Boundary_uu();
		uu_to_cc();
	}
	if (T_flag > 3.) // else //
	{
		TD23_New();
		store_du_duu();
	}
}

void CADSolver::TDMS34_New(void)
{

	if (ExtIter < 5)
	{
		TD34_New();
		store_du_duu();
	}
	else if (T_flag < 2.) //
	{

		Non_viscous_Splitting();
		if (Viscous_flag == 1)
			Viscous_Splitting();

		comput_du_LW_New();
		scheme_du_dx_New();
		if (Viscous_flag == 1)
			Viscous_Splitting_LW();

		TimeIntegration_TDMS34_ALW(iRK_Step);

		Boundary_uu();
		uu_to_cc();
	}
	else if (T_flag > 3.)
	{
		TD34_New();
	}
}

void CADSolver::TDMS34_Old(void)
{

	if (ExtIter < 5)
	{
		TD34_Old();
		store_du_duu();
	}
	else if (T_flag < 2.) //
	{

		Non_viscous_Splitting();

		// comput_du_LW_New();
		comput_du_LW_Old();
		scheme_du_dx_New();
		TimeIntegration_TDMS34_ALW(iRK_Step);

		Boundary_uu();
		uu_to_cc();
	}
	else if (T_flag > 3.)
	{
		TD34_Old();
	}
}
