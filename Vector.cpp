#include "Vector.h"
#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <limits>

namespace mat_vec {
    Vector operator*(double k, const Vector &v) {
        Vector newVector(v.size());
        for (size_t i = 0; i < newVector.size(); i++) {
            newVector[i] = v[i] * k;
        }
        return newVector;
    }

    Vector::Vector(size_t size, double value) : vector_ptr(new double[size]), size_(size) {
        for (size_t i = 0; i < this->size_; i++) {
            this->vector_ptr[i] = value;
        }
    }

    Vector::Vector(const Vector &src) : vector_ptr(new double[src.size()]), size_(src.size()) {
        for (size_t i = 0; i < this->size_; i++) {
            this->vector_ptr[i] = src.vector_ptr[i];
        }
    }

    Vector &Vector::operator=(const Vector &rhs) {
        delete[] this->vector_ptr;
        this->vector_ptr = new double[rhs.size()];
        this->size_ = rhs.size();
        for (size_t i = 0; i < this->size_; i++) {
            this->vector_ptr[i] = rhs.vector_ptr[i];
        }
        return *this;
    }

    Vector::~Vector() {
        delete[] this->vector_ptr;
    }

    size_t Vector::size() const {
        return (this->size_);
    }

    double Vector::operator[](size_t n) const {
        return (this->vector_ptr[n]);
    }

    double &Vector::operator[](size_t n) {
        return (this->vector_ptr[n]);
    }

    double Vector::norm() const {
        double L2Norm = 0;
        for (size_t i = 0; i < this->size_; i++) {
            L2Norm += (this->vector_ptr[i] * this->vector_ptr[i]);
        }
        return std::sqrt(L2Norm);
    }

    Vector Vector::normalized() const {
        double L2Norm = this->norm();
        Vector newVector(this->size_);
        for (size_t i = 0; i < newVector.size(); i++) {
            newVector[i] = (this->vector_ptr[i] / L2Norm);
        }
        return newVector;
    }

    void Vector::normalize() {
        double L2Norm = this->norm();
        for (size_t i = 0; i < this->size_; i++) {
            this->vector_ptr[i] /= L2Norm;
        }
    }

    Vector Vector::operator+(const Vector &rhs) const {
        if (this->size_ != rhs.size_) {
            Vector newVector(this->size_);
            for (size_t i = 0; i < newVector.size_; i++) {
                newVector.vector_ptr[i] = std::numeric_limits<double>::infinity();
            }
            return newVector;
        } else {
            Vector newVector(this->size_);
            for (size_t i = 0; i < this->size_; i++) {
                newVector[i] = this->vector_ptr[i] + rhs[i];
            }
            return newVector;
        }
    }

    Vector &Vector::operator+=(const Vector &rhs) {
        if (this->size_ != rhs.size_) {
            for (size_t i = 0; i < this->size_; i++) {
                this->vector_ptr[i] = std::numeric_limits<double>::infinity();
            }
            return *this;
        } else {
            for (size_t i = 0; i < this->size_; i++) {
                this->vector_ptr[i] += rhs[i];
            }
            return *this;
        }
    }

    Vector Vector::operator-(const Vector &rhs) const {
        if (this->size_ != rhs.size_) {
            Vector newVector(this->size_);
            for (size_t i = 0; i < this->size_; i++) {
                newVector.vector_ptr[i] = std::numeric_limits<double>::infinity();
            }
            return newVector;
        } else {
            Vector newVector(this->size_);
            for (size_t i = 0; i < this->size_; i++) {
                newVector[i] = this->vector_ptr[i] - rhs[i];
            }
            return newVector;
        }
    }

    Vector &Vector::operator-=(const Vector &rhs) {
        if (this->size_ != rhs.size_) {
            for (size_t i = 0; i < this->size_; i++) {
                this->vector_ptr[i] = std::numeric_limits<double>::infinity();
            }
            return *this;
        } else {
            for (size_t i = 0; i < this->size_; i++) {
                this->vector_ptr[i] -= rhs[i];
            }
            return *this;
        }
    }

    Vector Vector::operator^(const Vector &rhs) const {
        if (this->size_ != rhs.size_) {
            Vector newVector(this->size_);
            for (size_t i = 0; i < newVector.size_; i++) {
                newVector.vector_ptr[i] = std::numeric_limits<double>::infinity();
            }
            return newVector;
        } else {
            Vector newVector(this->size_);
            for (size_t i = 0; i < this->size_; i++) {
                newVector[i] = this->vector_ptr[i] * rhs[i];
            }
            return newVector;
        }
    }

    Vector &Vector::operator^=(const Vector &rhs) {
        if (this->size_ != rhs.size_) {
            for (size_t i = 0; i < this->size_; i++) {
                this->vector_ptr[i] = std::numeric_limits<double>::infinity();
            }
            return *this;
        } else {
            for (size_t i = 0; i < this->size_; i++) {
                this->vector_ptr[i] *= rhs[i];
            }
            return *this;
        }
    }

    double Vector::operator*(const Vector &rhs) const {
        if (this->size_ != rhs.size_) {
            return std::numeric_limits<double>::infinity();
        } else {
            double res = 0;
            for (size_t i = 0; i < this->size_; i++) {
                res += this->vector_ptr[i] * rhs[i];
            }
            return res;
        }
    }

    Vector Vector::operator*(double k) const {
        Vector newVector(this->size_);
        for (size_t i = 0; i < newVector.size(); i++) {
            newVector[i] = this->vector_ptr[i] * k;
        }
        return newVector;
    }

    Vector &Vector::operator*=(double k) {
        for (size_t i = 0; i < this->size_; i++) {
            this->vector_ptr[i] *= k;
        }
        return *this;
    }

    Vector Vector::operator/(double k) const {
        if (k != 0) {
            Vector newVector(this->size_);
            for (size_t i = 0; i < newVector.size(); i++) {
                newVector[i] = this->vector_ptr[i] / k;
            }
            return newVector;
        } else {
            Vector newVector(this->size_);
            for (size_t i = 0; i < newVector.size(); i++) {
                newVector[i] = std::numeric_limits<double>::infinity();
            }
            return newVector;
        }
    }

    Vector &Vector::operator/=(double k) {
        if (k != 0) {
            for (size_t i = 0; i < this->size_; i++) {
                this->vector_ptr[i] /= k;
            }
            return *this;
        } else {
            for (size_t i = 0; i < this->size_; i++) {
                this->vector_ptr[i] = std::numeric_limits<double>::infinity();
            }
            return *this;
        }
    }

    Vector Vector::operator*(const Matrix &mat) const {
        Vector newVector(this->size_);
        for (size_t i = 0; i < newVector.size_; i++) {
            for (size_t j = 0; j < newVector.size_; j++) {
                newVector[i] += mat.get(i, j) * this->vector_ptr[j];
            }
        }
        return newVector;
    }

    Vector &Vector::operator*=(const Matrix &mat) {
        Vector tempVector(this->size_, 0.0);
        for (size_t i = 0; i < tempVector.size_; i++) {
            for (size_t j = 0; j < tempVector.size_; j++) {
                tempVector[i] += mat.get(i, j) * this->vector_ptr[j];
            }
        }
        *this = tempVector;
        return *this;
    }

    bool Vector::operator==(const Vector &rhs) const {
        if (this->size_ != rhs.size_) {
            return false;
        }
        double eps = 0.00000000001;
        for (size_t i = 0; i < this->size_; i++) {
            if (std::abs(this->vector_ptr[i] - rhs[i]) > eps) {
                return false;
            }
        }
        return true;
    }

    bool Vector::operator!=(const Vector &rhs) const {
        return (!(*this == rhs));
    }
}
