#include "matrix.h"

Matrix::Matrix():m_rows(0), m_cols(0), m_data(0) {}

Matrix::Matrix(int rows, int cols):m_rows(rows), m_cols(cols), m_data(0) {
    this->m_data = (double*)calloc(m_rows * m_cols, sizeof(double));
}

Matrix::Matrix(const Matrix& other):m_rows(other.m_rows), m_cols(other.m_cols), m_data(0) {
    unsigned nrel = other.m_rows * other.m_cols;
    if(nrel) {
        this->m_data = (double*)malloc(nrel * sizeof(double));
        if(this->m_data) {
            memcpy(m_data, other.m_data, nrel * sizeof(double));
        }
    }
}

Matrix& Matrix::operator=(const Matrix& other) {
    if(this == &other) return *this;

    unsigned bytes_new = other.m_rows * other.m_cols * sizeof(double);
    unsigned bytes_cur = this->m_rows * this->m_cols * sizeof(double);
    if(bytes_new != bytes_cur) {
        if(this->m_data) {
            free(this->m_data);
        }
        this->m_data = 0;
        if(bytes_new) {
            this->m_data = (double*)malloc(bytes_new);
        }
    }
    this->m_rows = other.m_rows;
    this->m_cols = other.m_cols;
    if(this->m_data && other.m_data && bytes_new) {
        memcpy(this->m_data, other.m_data, bytes_new);
    }
    return *this;
}

Matrix::~Matrix() {
    if(this->m_data) {
        free(this->m_data);
    }
}

double& Matrix::at(int r, int c) {
    return this->m_data[r * this->m_cols + c];
}

double& Matrix::operator()(int r, int c) {
    return this->m_data[r * this->m_cols + c];
}

const double& Matrix::operator()(int r, int c) const {
    return this->m_data[r * this->m_cols + c];
}

void Matrix::print(const char* format) const {
    for(int i=0; i<this->m_rows; ++i) {
        for(int j=0; j<this->m_cols; ++j) {
            printf(format, this->at(i,j));
        }
        printf("\n");
    }
}

Matrix operator*(const Matrix &A, const Matrix &B) {
    int m = A.rows();
    int n = A.cols();
    int p = B.cols();

    if(n != B.rows()) {
        cerr << "Invalid matrix multiplication!" << endl;
        return Matrix();
    }

    Matrix C(m, p);
    int i,j,k;

    for(i=0; i<m; ++i) {
        for(j=0; j<p; ++j) {
            double sum = 0.0;
            for(k=0; k<n; ++k) {
                sum += A(i,k) * B(k,j);
            }
            C(i,j) = sum;
        }
    }

    return C;
}

Matrix operator+(const Matrix &A, const Matrix &B) {
    if(A.rows() != B.rows() || A.cols() != B.cols()) {
        cerr << "Invalid matrix addition!" << endl;
        return Matrix();
    }

    Matrix C(A.rows(), A.cols());
    int i,j;

    for(i=0; i<A.rows(); ++i) {
        for(j=0; j<A.cols(); ++j) {
            C(i,j) = A(i,j) + B(i,j);
        }
    }

    return C;
}

/*
 * Create matrix transpose
 */
Matrix Matrix::transpose() const {
    Matrix C(this->m_cols, this->m_rows);
    int i,j;

    for(i=0; i<this->m_cols; ++i) {
        for(j=0; j<this->m_rows; ++j) {
            C(i,j) = this->at(j,i);
        }
    }

    return C;
}

/*
 * Create zero matrix from this one
 */
Matrix Matrix::zeros_like() const {
    return Matrix(this->m_cols, this->m_rows);
}

/*
 * Self-transpose the matrix
 */
void Matrix::self_transpose() {
    double v;
    const int cols = this->m_cols;
    const int rows = this->m_rows;
    int i,j;

    // swap row and col setting
    this->m_cols = rows;
    this->m_rows = cols;

    // swap values
    for(i=0; i<rows; ++i) {
        for(j=i+1; j<cols; ++j) {
            v = this->at(i,j);
            this->at(i,j) = this->at(j,i);
            this->at(j,i) = v;
        }
    }
}