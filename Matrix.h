//
// Created by player402 on 02.05.24.
//

#ifndef MN_MATRIX_H
#define MN_MATRIX_H


#pragma once

#include <vector>
#include <cstdlib>
#include <ostream>

using namespace  std;

class Matrix
{
private:
    vector<double>elements;
    int rows;
    int cols;


public:

    enum InsertMode
    {
        ADD,
        SUBTRACT,
        MULTIPLY,
        DIVIDE,
        REPLACE,
    };

    void setElement(double element, int col, int row);

    double getElement(int col, int row) const;

    void insertSubMatrix(Matrix& source, int colUpperLeftElement, int rowUpperLeftElement, int colLowerRightElement, int rowLowerRightElement, int col, int row, InsertMode mode = InsertMode::REPLACE, double multipliedByScalar = 1.0);

    int getRows() const;

    int getCols() const;

    double norm() const;

    Matrix lower(int scope = 1) const;

    Matrix diag(int scope = 0) const;

    Matrix upper(int scope = 1) const;

    void elementWiseInverse();

    Matrix substituteForward(Matrix vector);

    Matrix substituteBackwards(Matrix vector);

    Matrix(int rows, int cols);

    Matrix(int rows, int cols, double elements []);

    Matrix(int rows, int cols, double elements);

    Matrix(int size, double elements=1.0);

    Matrix(Matrix& old);

    Matrix(Matrix&& old);

    Matrix& operator=(const Matrix& right);

    friend Matrix operator*(Matrix& left, Matrix& right);

    friend Matrix operator*(Matrix& left, double right);

    friend Matrix operator+(Matrix& left, Matrix& right);

    friend Matrix operator-(Matrix& left, Matrix& right);

    friend ostream& operator<<(ostream& os, Matrix& right);

};



#endif //MN_MATRIX_H
