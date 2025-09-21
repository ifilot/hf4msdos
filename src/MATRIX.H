#ifndef _MATRIX_H
#define _MATRIX_H

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <iostream.h>

/*
 * Custom matrix object; assumes matrices are smaller than 64 KiB by which
 * a near memory model can be used for them.
 */
class Matrix {

private:
    int m_rows; // number of rows
    int m_cols; // number of columns
    double *m_data; // pointer to first matrix element

public:
    /*
     * Default copy constructor
     */
    Matrix();

    /*
     * Build zero matrix
     */
    Matrix(int rows, int cols);

    /*
     * Copy constructor
     */
    Matrix(const Matrix& other);

    /*
     * Assignment constructor
     */
    Matrix& operator=(const Matrix& other);

    /*
     * Destructor
     */
    ~Matrix();

    /*
     * Access matrix element method
     */
    double& at(int r, int c);

    /*
     * Expose underlying data object
     */
    inline double* data() { return this->m_data; }

    /*
     * Alternative matrix access element method
     */
    double& operator()(int r, int c);

    /*
     * Non-altering matrix access element method
     */
    const double& operator()(int r, int c) const;

    /*
     * Output matrix
     */
    void print(const char* format) const;

    /*
     * Get number of rows
     */
    inline int rows() const {return this->m_rows;}

    /*
     * Get number of columns
     */
    inline int cols() const {return this->m_cols;}

    /*
     * Standard matrix multiplication
     */
    friend Matrix operator*(const Matrix &A, const Matrix &B);

    /*
     * Element-wise matrix addition
     */
    friend Matrix operator+(const Matrix &A, const Matrix &B);

    /*
     * Create matrix transpose
     */
    Matrix transpose() const;

    /*
     * Create zero matrix from this one
     */
    Matrix zeros_like() const;

    /*
     * Self-transpose the matrix
     */
    void self_transpose();

private:

};

#endif _MATRIX_H