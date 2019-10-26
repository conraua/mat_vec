#include "Base.h"
#include "Matrix.h"
#include "Vector.h"
#include <tuple>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <limits>
#include <cmath>

namespace mat_vec {
    Matrix::Matrix(size_t size, double value) : matrix_ptr(new double *[size]), rows_(size), cols_(size) {
        for (size_t i = 0; i < size; i++) {
            this->matrix_ptr[i] = new double[size];
        }
        for (size_t i = 0; i < size; i++) {
            for (size_t j = 0; j < size; j++) {
                this->matrix_ptr[i][j] = value;
            }
        }
    }

    Matrix Matrix::eye(size_t size) {
        Matrix newMatrix(size);
        for (size_t i = 0; i < size; i++) {
            for (size_t j = 0; j < size; j++) {
                if (i == j) {
                    newMatrix.matrix_ptr[i][j] = 1.0;
                } else {
                    newMatrix.matrix_ptr[i][j] = 0.0;
                }
            }
        }
        return newMatrix;
    }

    Matrix::Matrix(size_t rows, size_t cols, double value) : matrix_ptr(new double *[rows]), rows_(rows), cols_(cols) {
        for (size_t i = 0; i < this->rows_; i++) {
            this->matrix_ptr[i] = new double[this->cols_];
        }
        for (size_t i = 0; i < this->rows_; i++) {
            for (size_t j = 0; j < this->cols_; j++) {
                this->matrix_ptr[i][j] = value;
            }
        }
    }

    Matrix::Matrix(const Matrix &src) : matrix_ptr(new double *[src.rows_]), rows_(src.rows_), cols_(src.cols_) {
        for (size_t i = 0; i < this->rows_; i++) {
            this->matrix_ptr[i] = new double[this->cols_];
        }
        for (size_t i = 0; i < this->rows_; i++) {
            for (size_t j = 0; j < this->cols_; j++) {
                this->matrix_ptr[i][j] = src.matrix_ptr[i][j];
            }
        }
    }

    Matrix &Matrix::operator=(const Matrix &rhs) {
        for (size_t i = 0; i < this->rows_; ++i) {
            delete[] this->matrix_ptr[i];
        }
        delete[] this->matrix_ptr;
        this->rows_ = rhs.rows_;
        this->cols_ = rhs.cols_;
        this->matrix_ptr = new double*[this->rows_];
        for (size_t i = 0; i < this->rows_; i++) {
            this->matrix_ptr[i] = new double[this->cols_];
        }
        for (size_t i = 0; i < this->rows_; ++i) {
            for (size_t j = 0; j < this->cols_; ++j) {
                this->matrix_ptr[i][j] = rhs.matrix_ptr[i][j];
            }
        }
        return *this;
    }

    Matrix::~Matrix() {
        for (size_t i = 0; i < this->rows_; i++) {
            delete[] matrix_ptr[i];
        }
        delete[] matrix_ptr;
    }

    void Matrix::reshape(size_t rows, size_t cols) {
        Matrix tempMatrix(rows, cols, 0);
        double *temp = new double[this->rows_ * this->cols_];
        int number = 0;
        for (size_t i = 0; i < this->rows_; ++i) {
            for (size_t j = 0; j < this->cols_; ++j) {
                temp[number] = this->matrix_ptr[i][j];
                number++;
            }
        }
        number = 0;
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                tempMatrix.matrix_ptr[i][j] = temp[number];
                number++;
            }
        }
        delete[] temp;
        *this = tempMatrix;
    }

    std::pair<size_t, size_t> Matrix::shape() const {
        return std::pair<size_t, size_t>(this->rows_, this->cols_);
    }

    double Matrix::get(size_t row, size_t col) const {
        return (this->matrix_ptr[row][col]);
    }

    Matrix Matrix::operator+(const Matrix &rhs) const {
        if (this->rows_ != rhs.rows_ || this->cols_ != rhs.cols_) {
            Matrix newMatrix(this->rows_, this->cols_, 0);
            for (size_t i = 0; i < newMatrix.rows_; i++) {
                for (size_t j = 0; j < newMatrix.cols_; j++) {
                    newMatrix.matrix_ptr[i][j] = std::numeric_limits<double>::infinity();
                }
            }
            return newMatrix;
        } else {
            Matrix newMatrix(this->rows_, this->cols_, 0);
            for (size_t i = 0; i < newMatrix.rows_; i++) {
                for (size_t j = 0; j < newMatrix.cols_; j++) {
                    newMatrix.matrix_ptr[i][j] = this->matrix_ptr[i][j] + rhs.matrix_ptr[i][j];
                }
            }
            return newMatrix;
        }
    }

    Matrix &Matrix::operator+=(const Matrix &rhs) {
        if (this->rows_ != rhs.rows_ || this->cols_ != rhs.cols_) {
            for (size_t i = 0; i < this->rows_; i++) {
                for (size_t j = 0; j < this->cols_; j++) {
                    this->matrix_ptr[i][j] = std::numeric_limits<double>::infinity();
                }
            }
            return *this;
        } else {
            for (size_t i = 0; i < this->rows_; i++) {
                for (size_t j = 0; j < this->cols_; j++) {
                    this->matrix_ptr[i][j] += rhs.matrix_ptr[i][j];
                }
            }
            return *this;
        }
    }

    Matrix Matrix::operator-(const Matrix &rhs) const {
        if (this->rows_ != rhs.rows_ || this->cols_ != rhs.cols_) {
            Matrix newMatrix(this->rows_, this->cols_, 0);
            for (size_t i = 0; i < newMatrix.rows_; i++) {
                for (size_t j = 0; j < newMatrix.cols_; j++) {
                    newMatrix.matrix_ptr[i][j] = std::numeric_limits<double>::infinity();
                }
            }
            return newMatrix;
        } else {
            Matrix newMatrix(this->rows_, this->cols_, 0);
            for (size_t i = 0; i < newMatrix.rows_; i++) {
                for (size_t j = 0; j < newMatrix.cols_; j++) {
                    newMatrix.matrix_ptr[i][j] = this->matrix_ptr[i][j] - rhs.matrix_ptr[i][j];
                }
            }
            return newMatrix;
        }
    }

    Matrix &Matrix::operator-=(const Matrix &rhs) {
        if (this->rows_ != rhs.rows_ || this->cols_ != rhs.cols_) {
            for (size_t i = 0; i < this->rows_; i++) {
                for (size_t j = 0; j < this->cols_; j++) {
                    this->matrix_ptr[i][j] = std::numeric_limits<double>::infinity();
                }
            }
            return *this;
        } else {
            for (size_t i = 0; i < this->rows_; i++) {
                for (size_t j = 0; j < this->cols_; j++) {
                    this->matrix_ptr[i][j] -= rhs.matrix_ptr[i][j];
                }
            }
            return *this;
        }
    }

    Matrix Matrix::operator*(const Matrix &rhs) const {
        if (this->cols_ != rhs.rows_) {
            Matrix newMatrix(this->rows_, this->cols_, 0);
            for (size_t i = 0; i < newMatrix.rows_; i++) {
                for (size_t j = 0; j < newMatrix.cols_; j++) {
                    newMatrix.matrix_ptr[i][j] = std::numeric_limits<double>::infinity();
                }
            }
            return newMatrix;
        } else {
            Matrix newMatrix(this->rows_, rhs.cols_, 0);
            for (size_t i = 0; i < newMatrix.rows_; i++) {
                for (size_t j = 0; j < newMatrix.cols_; j++) {
                    for (size_t k = 0; k < this->cols_; k++)
                        newMatrix.matrix_ptr[i][j] += this->matrix_ptr[i][k] * rhs.matrix_ptr[k][j];
                }
            }
            return newMatrix;
        }
    }


    Matrix &Matrix::operator*=(const Matrix &rhs) {
        if (this->cols_ != rhs.rows_) {
            for (size_t i = 0; i < this->rows_; i++) {
                for (size_t j = 0; j < this->cols_; j++) {
                    this->matrix_ptr[i][j] = std::numeric_limits<double>::infinity();
                }
            }
            return *this;
        } else {
            Matrix tempMatrix(*this);
            reshape(this->rows_, rhs.cols_);
            for (size_t i = 0; i < this->rows_; i++) {
                for (size_t j = 0; j < this->cols_; j++) {
                    for (size_t k = 0; k < tempMatrix.cols_; k++)
                        this->matrix_ptr[i][j] += tempMatrix.matrix_ptr[i][k] * rhs.matrix_ptr[k][j];
                }
            }
            return *this;
        }
    }

    Matrix Matrix::operator*(double k) const {
        Matrix newMatrix(this->rows_, this->cols_, 0);
        for (size_t i = 0; i < newMatrix.rows_; i++) {
            for (size_t j = 0; j < newMatrix.cols_; j++) {
                newMatrix.matrix_ptr[i][j] = this->matrix_ptr[i][j] * k;
            }
        }
        return newMatrix;
    }

    Matrix &Matrix::operator*=(double k) {
        for (size_t i = 0; i < this->rows_; i++) {
            for (size_t j = 0; j < this->cols_; j++) {
                this->matrix_ptr[i][j] *= k;
            }
        }
        return *this;
    }

    Matrix Matrix::operator/(double k) const {
        if (k != 0) {
            Matrix newMatrix(this->rows_, this->cols_, 0);
            for (size_t i = 0; i < newMatrix.rows_; i++) {
                for (size_t j = 0; j < newMatrix.cols_; j++) {
                    newMatrix.matrix_ptr[i][j] = this->matrix_ptr[i][j] / k;
                }
            }
            return newMatrix;
        } else {
            Matrix newMatrix(this->rows_, this->cols_, 0);
            for (size_t i = 0; i < newMatrix.rows_; i++) {
                for (size_t j = 0; j < newMatrix.cols_; j++) {
                    newMatrix.matrix_ptr[i][j] = std::numeric_limits<double>::infinity();
                }
            }
            return newMatrix;
        }
    }

    Matrix &Matrix::operator/=(double k) {
        if (k != 0) {
            for (size_t i = 0; i < this->rows_; i++) {
                for (size_t j = 0; j < this->cols_; j++) {
                    this->matrix_ptr[i][j] /= k;
                }
            }
            return *this;
        } else {
            for (size_t i = 0; i < this->rows_; i++) {
                for (size_t j = 0; j < this->cols_; j++) {
                    this->matrix_ptr[i][j] = std::numeric_limits<double>::infinity();
                }
            }
            return *this;
        }
    }

    Matrix Matrix::transposed() const {
        Matrix newMatrix(this->cols_, this->rows_, 0);
        for (size_t i = 0; i < newMatrix.rows_; i++) {
            for (size_t j = 0; j < newMatrix.cols_; j++) {
                newMatrix.matrix_ptr[i][j] = this->matrix_ptr[j][i];
            }
        }
        return newMatrix;
    }

    void Matrix::transpose() {
        Matrix tempMatrix(*this);
        reshape(this->cols_, this->rows_);
        for (size_t i = 0; i < this->rows_; i++) {
            for (size_t j = 0; j < this->cols_; j++) {
                this->matrix_ptr[i][j] = tempMatrix.matrix_ptr[j][i];
            }
        }
    }

    double Matrix::det() const {
        double EPS = 0.00000000001;
        double det_ = 1;
        Matrix tempMatrix(*this);
        for (size_t i = 0; i < tempMatrix.rows_; i++) {
            size_t k = i;
            for (size_t j = i + 1; j < tempMatrix.rows_; j++)
                if (std::abs(tempMatrix.matrix_ptr[j][i]) > std::abs(tempMatrix.matrix_ptr[k][i]))
                    k = j;
            std::swap(tempMatrix.matrix_ptr[i], tempMatrix.matrix_ptr[k]);
            if (i != k)
                det_ = -det_;
            det_ *= tempMatrix.matrix_ptr[i][i];
            for (size_t j = i + 1; j < tempMatrix.rows_; j++)
                tempMatrix.matrix_ptr[i][j] /= tempMatrix.matrix_ptr[i][i];
            for (size_t j = 0; j < tempMatrix.rows_; j++)
                if (j != i && std::abs(tempMatrix.matrix_ptr[j][i]) > EPS)
                    for (size_t k = i + 1; k < tempMatrix.rows_; k++)
                        tempMatrix.matrix_ptr[j][k] -= tempMatrix.matrix_ptr[i][k] * tempMatrix.matrix_ptr[j][i];
        }
        return det_;
    }

    Matrix Matrix::inv() const {
        Matrix inverse(rows_);
        double det_ = getDet(matrix_ptr, rows_);
        if (rows_ == 1) {
            inverse.matrix_ptr[0][0] = 1 / det_;
            return inverse;
        }
        int sign = 1;
        double **minor__;
        minor__ = new double *[rows_];
        for (size_t i = 0; i < rows_; ++i) {
            minor__[i] = new double[rows_];
        }
        for (size_t i = 0; i < rows_; i++) {
            for (size_t j = 0; j < rows_; j++) {
                getMinor(rows_, matrix_ptr, minor__, i, j);
                sign = ((i + j) % 2 == 0) ? 1 : -1;
                inverse.matrix_ptr[j][i] = sign * (getDet(minor__, rows_ - 1)) / det_;
            }
        }
        for (size_t i = 0; i < rows_; ++i) {
            delete[] minor__[i];
        }
        delete[] minor__;
        return inverse;
    }

    Vector Matrix::operator*(const Vector &vec) const {
        Vector newVector(vec.size());
        for (size_t i = 0; i < newVector.size(); i++) {
            for (size_t j = 0; j < newVector.size(); j++) {
                newVector[i] += this->matrix_ptr[i][j] * vec[j];
            }
        }
        return newVector;
    }

    bool Matrix::operator==(const Matrix &rhs) const {
        if ((this->rows_ != rhs.rows_) || (this->cols_ != rhs.cols_)) {
            return false;
        } else {
            for (size_t i = 0; i < this->rows_; i++) {
                for (size_t j = 0; j < this->cols_; j++) {
                    if (this->matrix_ptr[i][j] != rhs.matrix_ptr[i][j]) {
                        return false;
                    }
                }
            }
            return true;
        }
    }

    bool Matrix::operator!=(const Matrix &rhs) const {
        return(!(*this == rhs));
    }

    double& Matrix::operator()(size_t row, size_t col)
    {
        if ((row < this->rows_) && (col < this->cols_)) {
            return this->matrix_ptr[row][col];
        } else {
            return this->matrix_ptr[0][0];
        }
    }

    double Matrix::getDet(double **a, size_t size) const {
        double d = 0;
        if (size == 1)
            return a[0][0];
        double sign = 1;
        double **temp;
        temp = new double *[size];
        for (size_t i = 0; i < size; ++i) {
            temp[i] = new double[size];
        }
        for (size_t y = 0; y < size; y++) {
            getMinor(size, a, temp, 0, y);
            d += sign * a[0][y] * getDet(temp, size - 1);
            sign = -sign;
        }
        for (size_t i = 0; i < size; ++i) {
            delete[] temp[i];
        }
        delete[] temp;
        return d;
    }

    void Matrix::getMinor(size_t size, double **c, double **minor_, size_t this_row, size_t this_col) const {
        size_t i = 0;
        size_t j = 0;
        for (size_t row = 0; row < size; row++) {
            for (size_t col = 0; col < size; col++) {
                if (row != this_row && col != this_col) {
                    minor_[i][j++] = c[row][col];
                    if (j == size - 1) {
                        j = 0;
                        i++;
                    }
                }
            }
        }
    }
}
