//
// Created by player402 on 02.05.24.
//

#ifndef MN_GAUSS_H
#define MN_GAUSS_H

#endif //MN_GAUSS_H

Matrix gauss(Matrix A, Matrix b, double desired_norm, vector<double>&solutions)
{
    Matrix x(b.getRows(), 1, 1.0);
    Matrix d = A.diag();
    Matrix l = A.lower();
    Matrix u = A.upper();
    Matrix dl = d+l;
    Matrix element2 = dl.substituteForward(b);
    double current_norm = x.norm();
    solutions.push_back(current_norm);
    while(current_norm>desired_norm)
    {
        Matrix ux = u*x;
        Matrix element1 = dl.substituteForward(ux);
        element1 = element1*(-1);
        x = element1 + element2;
        Matrix r= A*x;
        r = r - b;
        current_norm = r.norm();
        solutions.push_back(current_norm);
    }
    return x;
}