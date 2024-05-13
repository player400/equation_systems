//
// Created by player402 on 02.05.24.
//

#ifndef MN_GAUSS_H
#define MN_GAUSS_H

#endif //MN_GAUSS_H

Matrix jacobi(Matrix& A, Matrix& b, double desired_norm, vector<double>&solutions)
{
    Matrix x(b.getRows(), 1, 1.0);
    Matrix d = A.diag();
    Matrix l = A.lower();
    Matrix u = A.upper();
    Matrix lu = l+u;
    Matrix invd = d;
    invd.elementWiseInverse();
    Matrix element2 = invd*b;
    Matrix revinvd = invd*(-1);
    Matrix factor = revinvd*lu;

    double current_norm = x.norm();
    solutions.push_back(current_norm);
    while(current_norm>desired_norm)
    {
        if(solutions.size()>=1000)
        {
            break;
        }
        Matrix element1 = factor*x;
        x = element1 + element2;
        if(solutions.size()%1 == 0)
        {
            Matrix r = A*x;
            r = r - b;
            current_norm = r.norm();
        }

        solutions.push_back(current_norm);
    }
    return x;
}