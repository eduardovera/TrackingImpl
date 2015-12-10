#ifndef POINT_H
#define POINT_H

#include <dlib/matrix.h>

using namespace dlib;

class Point : public matrix<double, 1, 3> {

    public:

        Point() {

        }

        Point(matrix<double> m) {
            // TO DO: asserts
            double inv_m2 = 1.0 / m(2);

            this->operator ()(0) = m(0) * inv_m2;
            this->operator ()(1) = m(1) * inv_m2;
            this->operator ()(2) = 1.0;
        }

        Point(matrix<double, 1, 3> &m) {

            double inv_m2 = 1.0 / m(2);

            this->operator ()(0) = m(0) * inv_m2;
            this->operator ()(1) = m(1) * inv_m2;
            this->operator ()(2) = 1.0;

        }

        Point(double x, double y) {
            this->operator ()(0) = x;
            this->operator ()(1) = y;
            this->operator ()(2) = 1.0;
        }

        matrix<double, 3, 1> innerMatrix() {
            matrix<double, 3, 1> p;
            p(0) = this->operator ()(0);
            p(1) = this->operator ()(1);
            p(2) = this->operator ()(2);
            return p;
        }


};

Point operator*(const matrix<double> &T, const Point &p) {
    matrix<double, 1, 3> p_ = T * (matrix<double, 1, 3>) p;
    return (Point)(p_);
}

#endif // POINT_H

