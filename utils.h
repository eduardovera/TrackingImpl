#ifndef UTILS
#define UTILS

#include <dlib/matrix.h>
#include <Point.h>
#include <limits>

using namespace std;
using namespace dlib;

enum TYPE_NULL_SPACE_VECTOR {
    LEFT_NULL_SPACE_VECTOR,
    RIGHT_NULL_SPACE_VECTOR
};

matrix<double> getNullSpaceVector(const matrix<double> &m, TYPE_NULL_SPACE_VECTOR type) {
    matrix<double> u, s, v;
    matrix<double> null_space_vector;
    svd(m, u, s, v);
    double min_autovalue = numeric_limits<double>::max();
    int min_index = -1;
    for (int i = 0; i < s.nr(); i++) {
        if (s(i, i) < min_autovalue) {
            min_autovalue = s(i, i);
            min_index = i;
        }
    }
    if (type == RIGHT_NULL_SPACE_VECTOR) {
        null_space_vector = colm(v, min_index);
    } else if (type == LEFT_NULL_SPACE_VECTOR) {
        null_space_vector = colm(u, min_index);
    }
    return null_space_vector;
}

void force_restriction(matrix<double, 3, 3> &F) {
    matrix<double, 3, 3> u, s, v;
    svd(F, u, s, v);
    double min_autovalue = numeric_limits<double>::max();
    int min_index = -1;
    for (int i = 0; i < s.nr(); i++) {
        if (s(i, i) < min_autovalue) {
            min_autovalue = s(i, i);
            min_index = i;
        }
    }
    s(min_index, min_index) = 0;
    F = u * (s * trans(v));
}

double distanceBetweenPoints(const Point &p1, const Point &p2) {
    return ( ( (p1(0) - p2(0)) * (p1(0) - p2(0)) ) + ( (p1(1) - p2(1)) * (p1(1) - p2(1))) );
}

double getGeometricDistance(const matrix<double, 3, 1> &x, const matrix<double, 3, 1> &x_, const matrix<double, 3, 3> &F) {
    matrix<double> l, l_;
    l = F * x;
    l_ = trans(trans(x_) * F);

    double d1, d2;
    d1 = std::abs(dot(l, x_));
    d2 = std::abs(dot(l_, x));

    return d1 + d2;
}

matrix<double, 3, 3> crossMatrixForm(const matrix<double, 3, 1> &p) {
    matrix<double, 3, 3> crossMatrix;
    crossMatrix(0, 1) = -p(2);
    crossMatrix(1, 0) = p(2);

    crossMatrix(0, 2) = p(1);
    crossMatrix(2, 0) = -p(1);

    crossMatrix(1, 2) = -p(0);
    crossMatrix(2, 1) = p(0);
    return crossMatrix;
}

pair<matrix<double, 3, 4>, matrix<double, 3, 4>> getProjectionMatrixes(matrix<double, 3 ,3> F) {
    matrix<double, 3, 4> P;
    set_colm(P, range(0, 2)) = identity_matrix<double, 3>();
    set_colm(P, 3) = 0;

    matrix<double, 3, 4> P_;
    matrix<double, 3, 1> right_e = getNullSpaceVector(F, RIGHT_NULL_SPACE_VECTOR);

    set_colm(P_, range(0, 2)) = crossMatrixForm(right_e) * F;
    set_colm(P_, 3) = right_e;

    return pair<matrix<double, 3, 4>, matrix<double, 3, 4>>(P, P_);
}

matrix<double, 4, 1> DLT_Triangulation(const Point &x, const Point &x_, const pair<matrix<double, 3, 4>, const matrix<double, 3, 4>> &projectionMatrixes) {
    matrix<double, 4, 4> A;
    const matrix<double, 3, 4> P = projectionMatrixes.first;
    const matrix<double, 3, 4> P_ = projectionMatrixes.second;

    set_rowm(A, 0) = (x(0) * rowm(P, 2)) - rowm(P, 0);
    set_rowm(A, 1) = (x(1) * rowm(P, 2)) - rowm(P, 1);
    set_rowm(A, 2) = (x_(0) * rowm(P_, 2)) - rowm(P_, 0);
    set_rowm(A, 3) = (x_(1) * rowm(P_, 2)) - rowm(P_, 1);
    matrix<double, 4, 1> X = getNullSpaceVector(A, RIGHT_NULL_SPACE_VECTOR);

//    X /= X(3);

//    matrix<double> k = P_ * X;

//    cout << x_ << endl;
//    cout << X << endl;
//    cout << k / k(2) << endl;
//    getchar();

    return X;

}

double getQuadraticDistance(const std::vector<pair<Point, Point>> &matches, const pair<matrix<double, 3, 4>, matrix<double, 3, 4>> &projectionMatrixes) {
    double error = 0.0;
    matrix<double, 4, 1> X;
    matrix<double, 3, 1> p1_est;
    matrix<double, 3, 1> p2_est;
    for (const pair<Point, Point> &p : matches) {
        X = DLT_Triangulation(p.first, p.second, projectionMatrixes);
        p1_est = projectionMatrixes.first * X;
        p2_est = projectionMatrixes.second * X;
        error += distanceBetweenPoints(p.first, (Point)p1_est) + distanceBetweenPoints(p.second, (Point)p2_est);
    }
    return error;
}


#endif // UTILS

