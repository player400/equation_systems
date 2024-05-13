//
// Created by player402 on 02.05.24.
//

#ifndef MN_GAUSS_H
#define MN_GAUSS_H

#endif //MN_GAUSS_H

Matrix gauss(Matrix& A, Matrix& b, double desired_norm, vector<double>&solutions, int skipAssesments = 1)
{
    bool doSkipping = skipAssesments < 2 ? false : true;
    Matrix x(b.getRows(), 1, 1.0);
    Matrix r = A*x;
    r = r - b;
    double current_norm = r.norm();
    solutions.push_back(current_norm);
    while(current_norm>desired_norm)
    {
        if(solutions.size()>=2000)
        {
            break;
        }
        for(int i=0;i<x.getRows();i++)
        {
            double element = b.getElement(0, i);
            for(int j=0;j<x.getRows();j++)
            {
                if(j==i)
                {
                    continue;
                }
                element -= A.getElement(j, i)*x.getElement(0, j);
            }
            element /= A.getElement(i, i);
            x.setElement(element,0, i);
        }

        if(!doSkipping || solutions.size()%skipAssesments==0)
        {
            r = A*x;
            r = r - b;
            current_norm = r.norm();
        }
        solutions.push_back(current_norm);
    }
    return x;
}