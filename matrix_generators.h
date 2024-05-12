//
// Created by player402 on 12.05.24.
//

#ifndef MN_MATRIX_GENERATORS_H
#define MN_MATRIX_GENERATORS_H

#endif //MN_MATRIX_GENERATORS_H

Matrix generateAMatrix(int N, double a1, double a2, double a3)
{
    Matrix result(N, N, 0.0);
    for(int i = 0;i < N;i++)
    {
        for(int j = 0;j < N;j++)
        {
            if(i==j)
            {
                result.setElement(a1, i, j);
            }
            else if(abs(i-j)==1)
            {
                result.setElement(a2, i, j);
            }
            else if(abs(i-j)==2)
            {
                result.setElement(a3, i, j);
            }
        }
    }
    return result;
}

Matrix generateBMatrix(int N, int f)
{
    Matrix result(N, 1);
    for(int i = 1;i<=N;i++)
    {
        double element = sin(i*(f+1));
        result.setElement(element, 0, i-1);
    }
    return result;
}
