#define CATCH_CONFIG_MAIN
#include "Matrix.h"
#include "Vector.h"
#include "Base.h"
#include "catch.hpp"
#include <limits>

using namespace mat_vec;

const double INF = std::numeric_limits<double>::infinity();

TEST_CASE("Vector::operator=()") {
    Vector v1(5);
    Vector v2(1, 1.1);
    Vector v3(1, 1.1);
    v1 = v2;
    REQUIRE(v1 == v2);
}

TEST_CASE("Vector::size()") {
    Vector v1(5, 1.0);
     size_t size_ = v1.size();
     REQUIRE(size_ == 5);
}

TEST_CASE("Vector::operator[]()") {
    Vector v1(5, 5.0);
    double a = v1[1];
    double &b = v1[3];
    b = 0.0;
    REQUIRE(a != v1[3]);
}

TEST_CASE("Vector::norm()") {
    Vector v1 (4, 4.0);
    double length = v1.norm();
    REQUIRE(length == 8);
}

TEST_CASE("Vector::normalize()") {
    Vector v1 (4, 4.0);
    Vector v2(v1);
    v2 = v1.normalized();
    v1.normalize();
    REQUIRE(v1 == v2);
}

TEST_CASE("Vector::arithmetic_operators()") {
    Vector v1(5, 0.1);
    Vector v2(4, 0.0);
    Vector v3(5, 5.0);
    Vector v4 = v1 + v2;
    Vector v5(v1);
    REQUIRE(v4[0] == INF);
    v4 = v3 + v1;
    REQUIRE(v4[0] == 5.1);
    v4 = v1 - v2;
    REQUIRE(v4[0] == INF);
    v4 = v3 - v1;
    REQUIRE(v4[0] == 4.9);
    v4 = v1 ^ v2;
    REQUIRE(v4[0] == INF);
    v4 = v3 ^ v1;
    REQUIRE(v4[0] == 0.5);
    v2 += v1;
    REQUIRE(v2[0] == INF);
    v1 += v3;
    REQUIRE(v1[0] == 5.1);
    v2 -= v1;
    REQUIRE(v2[0] == INF);
    v1 -= v3;
    REQUIRE(v1 == v5);
    v2 ^= v1;
    REQUIRE(v2[0] == INF);
    v1 ^= v3;
    REQUIRE(v1 != v5);
}

TEST_CASE("Vector::operator*()") {
    Vector v1(5, 2.0);
    Vector v2(4, 1.0);
    double k = 0.5;
    Vector v3 = v1 * k;
    REQUIRE(v3[0] == 1);
    v3 = k * v1;
    REQUIRE(v3[0] == 1);
    v3 *= k;
    REQUIRE(v3[0] == 0.5);
    v3[0] = v1 * v1;
    REQUIRE(v3[0] == 20);
    v3[0] = v1 * v2;
    REQUIRE(v3[0] == INF);
}

TEST_CASE("Vector::operator/()") {
    Vector v1(5, 5.0);
    double k1 = 5.0;
    double k2 = 0.0;
    Vector v4(v1);
    v4 /= k1;
    REQUIRE(v4[0] == 1.0);
    v4 /= k2;
    REQUIRE(v4[0] == INF);
    v4 = v1 / k1;
    REQUIRE(v4[0] == 1.0);
    v4 = v1 / k2;
    REQUIRE(v4[0] == INF);
}

TEST_CASE("Vector::operator*(Matrix)") {
    Matrix m1(2, 0.0);
    Vector v2(2, 0.0);
    Vector v1(2, 0.0);
    Vector v3(v1);
    m1(0, 0) = 1;
    m1(0, 1) = 2;
    m1(1, 0) = 3;
    m1(1, 1) = 4;
    v2[0] = 1;
    v2[1] = 2;
    v1[0] = 5;
    v1[1] = 11;
    v3 = m1 * v2;
    REQUIRE(v3 == v1);
    v3 = v2 * m1;
    REQUIRE(v3 == v1);
    v2 *= m1;
    REQUIRE(v2 == v1);
}

TEST_CASE("Vector::operator==()") {
    Vector v1(3, 1.0);
    Vector v2(3, 1.0);
    Vector v3(2, 1.0);
    Vector v4(3, 0.0);
    REQUIRE(v1 == v2);
    REQUIRE(!(v1 == v3));
    REQUIRE(!(v2 == v4));
}

TEST_CASE("Vector::operator!=()") {
    Vector v1(3, 1.0);
    Vector v2(3, 1.0);
    Vector v3(2, 1.0);
    Vector v4(3, 0.0);
    REQUIRE(!(v1 != v2));
    REQUIRE(v1 != v3);
    REQUIRE(v2 != v4);
}

TEST_CASE("Matrix::eye()") {
    Matrix m1(3, 0.0);
    m1(0, 0) = 1.0;
    m1(1, 1) = 1.0;
    m1(2, 2) = 1.0;
    REQUIRE(m1 == Matrix::eye(3));
}

TEST_CASE("Matrix::operator=") {
    Matrix m1(3, 5.0);
    Matrix m2(5, 1.3);
    Matrix m3(5, 1.3);
    m1 = m2;
    REQUIRE(m1 == m3);
}

TEST_CASE("Matrix::reshape()") {
    Matrix m1(3, 2, 1.3);
    Matrix m2(2, 3, 1.3);
    m2(1, 2) = 0.0;
    m1(2, 1) = 0.0;
    m1.reshape(2, 3);
    REQUIRE(m1 == m2);
}

TEST_CASE("Matrix::shape()") {
    Matrix m1(3, 2, 0.0);
    Matrix m2(3, 2, 1.0);
    std::pair<size_t, size_t> s1 = m1.shape();
    std::pair<size_t, size_t> s2 = m2.shape();
    REQUIRE(s1 == s2);
}

TEST_CASE("Matrix::get()") {
    Matrix m1(5);
    m1(2, 2) = 5.0;
    double a = m1.get(2, 2);
    REQUIRE(a == 5.0);
}

TEST_CASE("Matrix::arithmetic_operators()") {
    Matrix m1(5, 2.0);
    Matrix m2(4, 2.0);
    Matrix m3(5, 1.0);
    Matrix m4 = m1 + m2;
    Matrix m5(m1);
    REQUIRE(m4(0, 0) == INF);
    m4 = m3 + m1;
    REQUIRE(m4(0, 0) == 3.0);
    m4 = m1 - m2;
    REQUIRE(m4(0, 0) == INF);
    m4 = m3 - m1;
    REQUIRE(m4(0, 0) == -1.0);
    m4 = m1 * m2;
    REQUIRE(m4(0, 0) == INF);
    m4 = m3 * m1;
    REQUIRE(m4(0, 0) == 10.0);
    m2 += m1;
    REQUIRE(m2(0, 0) == INF);
    m1 += m3;
    REQUIRE(m1(0, 0) == 3.0);
    m2 -= m1;
    REQUIRE(m2(0, 0) == INF);
    m1 -= m3;
    REQUIRE(m1 == m5);
    m2 *= m1;
    REQUIRE(m2(0, 0) == INF);
    m1 *= m3;
    REQUIRE(m1 != m5);
}

TEST_CASE("Matrix::operator*() & Matrix::operator/()") {
    Matrix m1(5, 1.0);
    Matrix m2(m1);
    double k1 = 0.5;
    double k2 = 0.0;
    m2 = m1 * k1;
    REQUIRE(m2(0, 0) == 0.5);
    m2 = m1 / k1;
    REQUIRE(m2(0, 0) == 2.0);
    m1 *= k1;
    REQUIRE(m1(0,0) == 0.5);
    m1 /= k1;
    REQUIRE(m1(0,0) == 1.0);
    m1 /= k2;
    REQUIRE(m1(0, 0) == INF);
    m2 = m1 / k2;
    REQUIRE(m2(0, 0) == INF);
}

TEST_CASE("Matrix::transposed()") {
    Matrix m1(5, 0.0);
    Matrix m2(5, 0.0);
    Matrix m3(3, 0.0);
    m1(0, 1) = 2;
    m1(2, 1) = 2;
    m2(1, 2) = 2;
    m2(1, 0) = 2;
    m3 = m1.transposed();
    REQUIRE(m3 == m2);
    m1.transpose();
    REQUIRE(m1 == m2);
}

TEST_CASE("Matrix::det() & Matrix::inv()") {
    Matrix m1(2, 0.0);
    Matrix m2(2, 0.0);
    Matrix m3(2, 0.0);
    m1(0, 0) = 1;
    m1(0, 1) = 2;
    m1(1, 0) = 3;
    m1(1, 1) = 4;
    m2(0, 0) = 10;
    m2(0, 1) = 20;
    m2(1, 0) = 30;
    m2(1, 1) = 40;
    m3(0, 0) = 11;
    m3(0, 1) = 22;
    m3(1, 0) = 33;
    m3(1, 1) = 44;
    REQUIRE((m1.det() == -2));
    REQUIRE((m1 * m1.inv() == Matrix::eye(2)));
    Matrix m4(1, 5.0);
    REQUIRE((m4 * m4.inv() == Matrix::eye(1)));
}

TEST_CASE("Matrix::operator==()") {
    Matrix m1(3, 1.0);
    Matrix m2(3, 1.0);
    Matrix m3(2, 1.0);
    Matrix m4(3, 0.0);
    REQUIRE(m1 == m2);
    REQUIRE(!(m1 == m3));
    REQUIRE(!(m2 == m4));
}

TEST_CASE("Matrix::operator!=()") {
    Matrix m1(3, 1.0);
    Matrix m2(3, 1.0);
    Matrix m3(2, 1.0);
    Matrix m4(3, 0.0);
    REQUIRE(!(m1 != m2));
    REQUIRE(m1 != m3);
    REQUIRE(m2 != m4);
}

TEST_CASE("Matrix::operator()()") {
    Matrix m1(5, 5.0);
    m1(0, 0) = 0.0;
    double a = m1(9, 9);
    REQUIRE(a == m1(0, 0));
}