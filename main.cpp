#include <iostream>
#include <cmath>

#include "Matrix.h"

using namespace std;

#define F 3
#define E 0
#define C 5
#define D 3
#define DESIRED_NORM 10e-9

#include "gauss.h"
#include "jacobi.h"
#include "fstream"

#include "matrix_generators.h"

void iterativeMethods()
{
    int N = 9*C*D;
    Matrix A = generateAMatrix(N, 5+E, -1, -1);
    Matrix b = generateBMatrix(N, F);

    vector<double> solutionsGauss;
    Matrix xGauss = gauss(A, b, DESIRED_NORM, solutionsGauss);

    vector<double> solutionsJacobi;
    Matrix xJacobi = jacobi(A, b, DESIRED_NORM, solutionsJacobi);

    ofstream gaussStepsFile("graph_data/gauss_steps.txt");
    cout<<"Gauss-Seidel solution:"<<endl<<xGauss<<endl;
    cout<<"Number of steps: "<<solutionsGauss.size()<<endl;
    cout<<"Norm during Gauss-Seidel iterative steps:"<<endl;
    for(int i=0;i<solutionsGauss.size();i++)
    {
        cout<<solutionsGauss[i]<<endl;
        gaussStepsFile<<solutionsGauss[i]<<endl;
    }
    cout<<endl;

    ofstream jacobiStepsFile("graph_data/jacobi_steps.txt");
    cout<<"Jacobi solution:"<<endl<<xGauss<<endl<<endl;
    cout<<"Number of steps: "<<solutionsJacobi.size()<<endl;
    cout<<"Norm during Jacobi iterative steps:"<<endl;
    for(int i=0;i<solutionsJacobi.size();i++)
    {
        cout<<solutionsJacobi[i]<<endl;
        jacobiStepsFile<<solutionsJacobi[i]<<endl;
    }
    cout<<endl;

    Matrix dif = xGauss-xJacobi;
    cout<<"Gauss-Seidel vs Jacobi difference:"<<endl<<dif<<endl<<"Norm: "<<dif.norm()<<endl<<endl;
}

int main() {
    iterativeMethods();
    return 0;
}
