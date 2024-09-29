/* A Novel Lax-Wendroff Type Procedure for Two-Derivative Time-Stepping Schemes 
on Euler and Navier-Stokes Equations.
 * Copyright (C) 2023
 * Author QinXueyu
 */

#include "xBurgers.hpp"

using namespace std;

int main()
{

	cout << "----------------------------------------" << endl;
	cout << "------------Start  Order----------------" << endl;

	CSolver *realsolver = NULL;

	realsolver = new CADSolver();

	realsolver->Driver();

	if (realsolver != NULL)
	{
		delete realsolver;
	}

	cout << "----------------------------------------" << endl;
	cout << "----------------end-------------------" << endl;
	cin.get();
	return 0;
}