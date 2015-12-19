#ifndef SOLVER
#define SOLVER
using namespace std;

const double PI = acos(-1);

int solve_deg2(double a, double b, double c, double & x1, double & x2) {
    double delta = b * b - 4 * a * c;

    if (delta < 0) return 0;

    double inv_2a = 0.5 / a;

    if (delta == 0) {
        x1 = -b * inv_2a;
        x2 = x1;
        return 1;
    }

    double sqrt_delta = sqrt(delta);
    x1 = (-b + sqrt_delta) * inv_2a;
    x2 = (-b - sqrt_delta) * inv_2a;
    return 2;
}

int solve_deg3(double a, double b, double c, double d, double & x0, double & x1, double & x2) {
    if (a == 0) {
        // Solve second order sytem
        if (b == 0) {
            // Solve first order system
            if (c == 0)
                return 0;

            x0 = -d / c;
            return 1;
        }

        x2 = 0;
        return solve_deg2(b, c, d, x0, x1);
    }

    // Calculate the normalized form x^3 + a2 * x^2 + a1 * x + a0 = 0
    double inv_a = 1. / a;
    double b_a = inv_a * b, b_a2 = b_a * b_a;
    double c_a = inv_a * c;
    double d_a = inv_a * d;

    // Solve the cubic equation
    double Q = (3 * c_a - b_a2) / 9;
    double R = (9 * b_a * c_a - 27 * d_a - 2 * b_a * b_a2) / 54;
    double Q3 = Q * Q * Q;
    double D = Q3 + R * R;
    double b_a_3 = (1. / 3.) * b_a;

    if (Q == 0) {
        if (R == 0) {
            x0 = x1 = x2 = -b_a_3;
            return 3;
        } else {
            x0 = pow(2 * R, 1 / 3.0) - b_a_3;
            return 1;
        }
    }

    if (D <= 0) {
        // Three real roots
        double theta = acos(R / sqrt(-Q3));
        double sqrt_Q = sqrt(-Q);
        x0 = 2 * sqrt_Q * cos(theta / 3.0) - b_a_3;
        x1 = 2 * sqrt_Q * cos((theta + 2 * PI) / 3.0) - b_a_3;
        x2 = 2 * sqrt_Q * cos((theta + 4 * PI) / 3.0) - b_a_3;

        return 3;
    }

    // D > 0, only one real root
    double AD = pow(fabs(R) + sqrt(D), 1.0 / 3.0) * (R > 0 ? 1 : (R < 0 ? -1 : 0));
    double BD = (AD == 0) ? 0 : -Q / AD;

    // Calculate the only real root
    x0 = AD + BD - b_a_3;

    return 1;
}

std::vector<matrix<double, 3, 3>> solve(const matrix<double, 3, 3> &F1, const matrix<double, 3, 3> &F2) {
    double f11 = F1(0, 0);
    double f12 = F1(0, 1);
    double f13 = F1(1, 2);
    double f21 = F1(1, 0);
    double f22 = F1(1, 1);
    double f23 = F1(1, 2);
    double f31 = F1(2, 0);
    double f32 = F1(2, 1);
    double f33 = F1(2, 2);
    double g11 = F2(0, 0);
    double g12 = F2(0, 1);
    double g13 = F2(1, 2);
    double g21 = F2(1, 0);
    double g22 = F2(1, 1);
    double g23 = F2(1, 2);
    double g31 = F2(2, 0);
    double g32 = F2(2, 1);
    double g33 = F2(2, 2);

    double a, b, c, d;

    d = -(g13*g22*g31) + (g12*g23*g31) + (g13*g21*g32) - (g11*g23*g32) - (g12*g21*g33) + (g11*g22*g33);

    c = -(f33*g12*g21) + (f32*g13*g21) + (f33*g11*g22) - (f31*g13*g22) - (f32*g11*g23) + (f31*g12*g23) +
         (f23*g12*g31) - (f22*g13*g31) - (f13*g22*g31) + (3*g13*g22*g31) + (f12*g23*g31) - (3*g12*g23*g31) -
         (f23*g11*g32) + (f21*g13*g32) + (f13*g21*g32) - (3*g13*g21*g32) - (f11*g23*g32) + (3*g11*g23*g32) +
         (f22*g11*g33) - (f21*g12*g33) - (f12*g21*g33) + (3*g12*g21*g33) + (f11*g22*g33) - (3*g11*g22*g33);

    b = -(f23*f32*g11) + (f22*f33*g11) + (f23*f31*g12) - (f21*f33*g12) - (f22*f31*g13) + (f21*f32*g13) +
         (f13*f32*g21) - (f12*f33*g21) + (2*f33*g12*g21) - (2*f32*g13*g21) - (f13*f31*g22) + (f11*f33*g22) -
         (2*f33*g11*g22) + (2*f31*g13*g22) + (f12*f31*g23) - (f11*f32*g23) + (2*f32*g11*g23) - (2*f31*g12*g23) -
         (f13*f22*g31) + (f12*f23*g31) - (2*f23*g12*g31) + (2*f22*g13*g31) + (2*f13*g22*g31) - (3*g13*g22*g31) -
         (2*f12*g23*g31) + (3*g12*g23*g31) + (f13*f21*g32) - (f11*f23*g32) + (2*f23*g11*g32) - (2*f21*g13*g32) -
         (2*f13*g21*g32) + (3*g13*g21*g32) + (2*f11*g23*g32) - (3*g11*g23*g32) - (f12*f21*g33) + (f11*f22*g33) -
         (2*f22*g11*g33) + (2*f21*g12*g33) + (2*f12*g21*g33) - (3*g12*g21*g33) - (2*f11*g22*g33) + (3*g11*g22*g33);

    a = -(f13*f22*f31) + (f12*f23*f31) + (f13*f21*f32) - (f11*f23*f32) - (f12*f21*f33) + (f11*f22*f33) + (f23*f32*g11) -
         (f22*f33*g11) - (f23*f31*g12) + (f21*f33*g12) + (f22*f31*g13) - (f21*f32*g13) - (f13*f32*g21) + (f12*f33*g21) -
         (f33*g12*g21) + (f32*g13*g21) + (f13*f31*g22) - (f11*f33*g22) + (f33*g11*g22) - (f31*g13*g22) - (f12*f31*g23) +
         (f11*f32*g23) - (f32*g11*g23) + (f31*g12*g23) + (f13*f22*g31) - (f12*f23*g31) + (f23*g12*g31) - (f22*g13*g31) -
         (f13*g22*g31) + (g13*g22*g31) + (f12*g23*g31) - (g12*g23*g31) - (f13*f21*g32) + (f11*f23*g32) - (f23*g11*g32) +
         (f21*g13*g32) + (f13*g21*g32) - (g13*g21*g32) - (f11*g23*g32) + (g11*g23*g32) + (f12*f21*g33) - (f11*f22*g33) +
         (f22*g11*g33) - (f21*g12*g33) - (f12*g21*g33) + (g12*g21*g33) + (f11*g22*g33) - (g11*g22*g33);

    std::vector<matrix<double, 3, 3>> F_vec;
    double x0, x1, x2;
    int roots_size = solve_deg3(a, b, c, d, x0, x1, x2);
    if (roots_size == 1) {
        F_vec.push_back((x0 * F1) + ((1 - x0) * F2));
    }
    if (roots_size == 2) {
        F_vec.push_back((x1 * F1) + ((1 - x1) * F2));
    }
    if (roots_size == 3) {
        F_vec.push_back((x2 * F1) + ((1 - x2) * F2));
    }
    return F_vec;
}

#endif // SOLVER

