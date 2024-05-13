//
// Created by player402 on 02.05.24.
//

#include <valarray>
#include "Matrix.h"

void Matrix::setElement(double element, int col, int row)
{
    elements[row * cols + col] = element;
}

void Matrix::insertSubMatrix(Matrix& source, int colUpperLeftElement, int rowUpperLeftElement, int colLowerRightElement, int rowLowerRightElement, int col, int row, InsertMode mode,  double multipliedByScalar)
{
    for(int i=rowUpperLeftElement, currentRow = row;i<=rowLowerRightElement;i++, currentRow++)
    {
        for(int j = colUpperLeftElement, currentCol = col; j<=colLowerRightElement; j++, currentCol++)
        {
            double toInsert=source.getElement(j, i)*multipliedByScalar;
            double element;
            if(mode!=REPLACE)
            {
                element = getElement(currentCol, currentRow);
                if(mode==ADD)
                {
                    element+=toInsert;
                }
                if(mode==SUBTRACT)
                {
                    element-=toInsert;
                }
                if(mode==MULTIPLY)
                {
                    element*=toInsert;
                }
                if(mode==DIVIDE)
                {
                    element/=toInsert;
                }
            }
            else
            {
                element = toInsert;
            }
            setElement(element,currentCol, currentRow);
        }
    }
}

double Matrix::getElement(int col, int row) const
{
    return elements[row * cols + col];
}

int Matrix::getRows() const
{
    return rows;
}

int Matrix::getCols() const
{
    return cols;
}

double Matrix::norm() const
{
    double sum=0;
    for(int i = 0;i < elements.size();i++)
    {
        sum+=elements[i]*elements[i];
    }
    return sqrt(sum);
}

Matrix::Matrix(int rows, int cols) :elements(cols*rows)
{
    this->rows = rows;
    this->cols = cols;
}

Matrix::Matrix(int rows, int cols, double elements []): Matrix(rows, cols)
{
    for(int i = 0;i < this->elements.size();i++)
    {
        this->elements[i] = elements[i];
    }
}


Matrix::Matrix(int size, double elements):Matrix(size, size) {
    for(int i=0 ;i<size;i++)
    {
        setElement(elements, i, i);
    }
}

Matrix::Matrix(int rows, int cols, double elements ): Matrix(rows, cols)
{
    for(int i = 0;i < this->elements.size();i++)
    {
        this->elements[i] = elements;
    }
}

Matrix::Matrix(Matrix& old)
{
    cols = old.getCols();
    rows = old.getRows();
    elements.resize(cols * rows);
    for (int i = 0; i < elements.size(); i++)
    {
        elements[i] = old.elements[i];
    }
}

Matrix::Matrix(Matrix&& old)
{
    cols = old.getCols();
    rows = old.getRows();
    elements.resize(cols * rows);
    for (int i = 0; i < elements.size(); i++)
    {
        elements[i] = old.elements[i];
    }
}

Matrix& Matrix::operator=(const Matrix& right)
{
    cols = right.getCols();
    rows = right.getRows();
    elements.resize(cols * rows);
    for (int i = 0; i < elements.size(); i++)
    {
        elements[i] = right.elements[i];
    }
    return *this;
}

Matrix operator*(Matrix& left, Matrix& right)
{
    int new_rows = left.getRows();
    int new_cols = right.getCols();
    int subvector_number = right.getRows();
    int offset = 0;
    Matrix result(new_rows, new_cols);
    for (int i = 0; i < new_cols; i++)
    {
        for (int k = 0; k < new_rows; k++)
        {
            double sum =0.0;
            for (int j = 0; j < subvector_number; j++)
            {
                double scalar = right.getElement(i, j);
                sum+=scalar * left.getElement(j, k);

            }
            result.setElement(sum,i, k );
        }
    }
    return result;
}

ostream& operator<<(ostream& os, Matrix& right)
{
    os<<"[";
    for(int i = 0;i<right.rows;i++)
    {
        for(int j = 0;j<right.cols;j++)
        {
            os<<right.getElement(j, i);
            if(j != right.cols-1)
            {
                os<<", ";
            }
        }
        if(i != right.rows-1)
        {
            os<<";"<<endl;
        }
    }
    os<<"]";
    return os;
}

Matrix operator+(Matrix &left, Matrix &right) {
    Matrix result(left.rows, left.cols);
    for(int i=0;i<result.elements.size();i++)
    {
        result.elements[i] = left.elements[i]+right.elements[i];
    }
    return result;
}

Matrix Matrix::lower(int scope) const {
    Matrix result(rows, cols, 0.0);
    for(int i=0;i<rows;i++)
    {
        for(int j=0;j<cols;j++)
        {
            if(i-j >= scope)
            {
                result.setElement(getElement(j, i), j, i);
            }
        }
    }
    return result;
}

Matrix Matrix::diag(int scope) const {
    Matrix result(rows, cols, 0.0);
    for(int i=0;i<rows;i++)
    {
        for(int j=0;j<cols;j++)
        {
            if(abs(i-j) <= scope)
            {
                result.setElement(getElement(j, i), j, i);
            }
        }
    }
    return result;
}

Matrix Matrix::upper(int scope) const {
    Matrix result(rows, cols, 0.0);
    for(int i=0;i<rows;i++)
    {
        for(int j=0;j<cols;j++)
        {
            if(j-i >= scope)
            {
                result.setElement(getElement(j, i), j, i);
            }
        }
    }
    return result;
}

Matrix Matrix::substituteForward(Matrix vector) {
    Matrix copy(*this);

    for(int i=0; i<cols;i++)
    {
        for(int j=i+1; j<rows;j++)
        {
            double ratio = copy.getElement(i, j)/copy.getElement(i, i);
            copy.setElement(copy.getElement(i, j)-copy.getElement(i, i)*ratio,i, j );
            vector.setElement(vector.getElement(0, j)-vector.getElement(0, i)*ratio, 0, j);
        }
    }
    for(int i=0; i<cols;i++)
    {
        vector.setElement(vector.getElement(0, i)/copy.getElement(i, i),0, i);
    }
    return vector;
}

Matrix Matrix::substituteBackwards(Matrix vector) {
    Matrix copy(*this);

    for(int i=cols-1; i>0;i--)
    {
        for(int j=i-1; j>=0;j--)
        {
            double ratio = copy.getElement(i, j)/copy.getElement(i, i);
            copy.setElement(copy.getElement(i, j)-copy.getElement(i, i)*ratio,i, j );
            vector.setElement(vector.getElement(0, j)-vector.getElement(0, i)*ratio, 0, j);
        }
    }
    for(int i=0; i<cols;i++)
    {
        vector.setElement(vector.getElement(0, i)/copy.getElement(i, i),0, i);
    }
    return vector;
}

void Matrix::elementWiseInverse() {
    for(int i=0;i<elements.size();i++)
    {
        if(elements[i]!=0)
        {
            elements[i] = 1/elements[i];
        }
    }
}

Matrix operator-(Matrix &left, Matrix &right) {
    Matrix result(left.rows, left.cols);
    for(int i=0;i<result.elements.size();i++)
    {
        result.elements[i] = left.elements[i]-right.elements[i];
    }
    return result;
}

Matrix operator*(Matrix& left, double right) {
    Matrix result(left);
    for(int i=0;i<result.elements.size();i++)
    {
        result.elements[i] = result.elements[i]*right;
    }
    return result;
}

Matrix operator*(double right, Matrix& left) {
    return left*right;
}