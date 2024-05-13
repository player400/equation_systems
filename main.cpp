#include <iostream>
#include <cmath>
#include <chrono>
#include "fstream"

using Time = std::chrono::steady_clock;

#include "Matrix.h"

using namespace std;

#define F 3
#define E 0
#define C 5
#define D 3
#define DESIRED_NORM 10e-9

const int speedTestMatrixSizes [] = {200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000};
const int skipping [] = {1, 5, 10, 30, 100, 300};

#include "gauss.h"
#include "jacobi.h"
#include "lu.h"

#include "matrix_generators.h"

void iterativeMethods(Matrix& A, Matrix& b, string folder_name)
{
    vector<double> solutionsGauss;
    Matrix xGauss = gauss(A, b, DESIRED_NORM, solutionsGauss);

    vector<double> solutionsJacobi;
    Matrix xJacobi = jacobi(A, b, DESIRED_NORM, solutionsJacobi);

    ofstream gaussStepsFile("graph_data/"+folder_name+"/gauss_steps.txt");
    cout<<"Gauss-Seidel solution:"<<endl<<xGauss<<endl;
    cout<<"Number of steps: "<<solutionsGauss.size()<<endl;
    cout<<"Residuum vector norm during Gauss-Seidel iterative steps:"<<endl;
    for(int i=0;i<solutionsGauss.size();i++)
    {
        cout<<solutionsGauss[i]<<endl;
        gaussStepsFile<<solutionsGauss[i]<<endl;
    }
    cout<<endl;

    ofstream jacobiStepsFile("graph_data/"+folder_name+"/jacobi_steps.txt");
    cout<<"Jacobi solution:"<<endl<<xGauss<<endl<<endl;
    cout<<"Number of steps: "<<solutionsJacobi.size()<<endl;
    cout<<"Residuum vector norm during Jacobi iterative steps:"<<endl;
    for(int i=0;i<solutionsJacobi.size();i++)
    {
        cout<<solutionsJacobi[i]<<endl;
        jacobiStepsFile<<solutionsJacobi[i]<<endl;
    }
    cout<<endl;

    Matrix dif = xGauss-xJacobi;
    cout<<"Gauss-Seidel vs Jacobi difference:"<<endl<<dif<<endl<<"Norm: "<<dif.norm()<<endl<<endl;
}



void directMethod(Matrix& A, Matrix& b)
{
    Matrix x = lu(A, b);
    cout<<"LU factorisation solution:"<<endl<<x<<endl<<endl;
    Matrix r = A*x;
    r = r - b;
    cout<<"LU solution residuum vector norm: "<<r.norm()<<endl<<endl;
}

void speedComparison(int skip, string filename)
{
    vector<int>timeGauss;
    vector<int>timeJacobi;
    vector<int>timeLU;
    vector<double>sols;
    for(int i=0;i<10;i++)
    {
        int N = speedTestMatrixSizes[i];
        Matrix A = generateAMatrix(N, 5+E, -1, -1);
        Matrix b = generateBMatrix(N, F);
        auto time = Time::now();
        Matrix x = gauss(A, b, DESIRED_NORM, sols, skip);
        auto now = Time::now();
        timeGauss.push_back(chrono::duration_cast<chrono::microseconds>(now-time).count());
        time = Time::now();
        x = jacobi(A, b, DESIRED_NORM, sols, skip);
        now = Time::now();
        timeJacobi.push_back(chrono::duration_cast<chrono::microseconds>(now-time).count());
        time = Time::now();
        x = lu(A, b);
        now = Time::now();
        timeLU.push_back(chrono::duration_cast<chrono::microseconds>(now-time).count());
    }
    ofstream file("graph_data/"+filename);
    cout<<endl<<"Skipping: "<<skip<<endl;
    cout<<"Matrix size, Time Gauss-Seidel, Time Jacobi, Time LU"<<endl;
    file<<"Matrix size, Time Gauss-Seidel, Time Jacobi, Time LU"<<endl;
    for(int i=0;i<timeLU.size();i++)
    {
        cout<<speedTestMatrixSizes[i]<<", "<<timeGauss[i]<<", "<<timeJacobi[i]<<", "<<timeLU[i]<<endl;
        file<<speedTestMatrixSizes[i]<<", "<<timeGauss[i]<<", "<<timeJacobi[i]<<", "<<timeLU[i]<<endl;
    }
}

void speedTest()
{
    cout<<endl<<endl<<"===== SPEED TEST ====="<<endl;
    speedComparison(30, "speed_test_skip_30.txt");
    speedComparison(1, "speed_test_skip_1");

}

void secondSystem()
{
    cout<<endl<<endl<<"===== SECOND SYSTEM OF EQUATIONS ====="<<endl;
    int N = 9*C*D;
    Matrix A = generateAMatrix(N, 3, -1,-1);
    Matrix b = generateBMatrix(N, F);
    iterativeMethods(A, b, "system2");
    directMethod(A, b);
}

void firstSystem()
{
    cout<<endl<<endl<<"===== FIRST SYSTEM OF EQUATIONS ====="<<endl;
    int N = 9*C*D;
    Matrix A = generateAMatrix(N, 5+E, -1, -1);
    Matrix b = generateBMatrix(N, F);
    iterativeMethods(A, b, "system1");
}

int main() {
    firstSystem();
    secondSystem();
    speedTest();
    return 0;
}
