//
// Created by player402 on 12.05.24.
//

#ifndef MN_LU_H
#define MN_LU_H

#endif //MN_LU_H

Matrix lu(Matrix& A, Matrix& b)
{
    int m = A.getCols();
    Matrix U(A);
    Matrix L(m);
    for(int i=2;i<=m;i++)
    {
        for(int j=1;j<i;j++)
        {
            L.setElement(U.getElement(j-1, i-1)/U.getElement(j-1, j-1), j-1, i-1);
            U.insertSubMatrix(U,0, j-1, m-1, j-1, 0, i-1,Matrix::InsertMode::SUBTRACT, L.getElement(j-1, i-1));
        }
    }
    Matrix y = L.substituteForward(b);
    Matrix x = U.substituteBackwards(y);
    return x;
}