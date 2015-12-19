#include <iostream>
#include <CImg.h>
#include <iostream>
#include <fstream>
#include <settings.h>
#include <utils.h>
#include <dlib/matrix.h>
#include <Point.h>
#include <chrono>
#include <solver.h>

using namespace std;

struct eigenvalue_sorter {

    public:
        double value;
        int index;

        eigenvalue_sorter(const double value, const int index) {
            this->index = index;
            this->value = value;
        }

        inline bool operator < (eigenvalue_sorter e) const {
            return value < e.value;
        }

};

matrix<double> getFundamental7Points(Tracking_settings &settings, const int &iterations = 1) {
    matrix<double> U, S, V;
    matrix<double, 9, 1> f1, f2;
    matrix<double, 3, 3> F1, F2;
    matrix<double> A(settings.matches.size(), 9);
    double index_min_error;
    std::vector<matrix<double, 3, 3>> fundamentals;

    int itr = 0;
    while (itr != iterations) {
        for (long j = 0; j < A.nr(); j++) {
            A(j, 0) = settings.matches[j].second(0) * settings.matches[j].first(0);
            A(j, 1) = settings.matches[j].second(0) * settings.matches[j].first(1);
            A(j, 2) = settings.matches[j].second(0);
            A(j, 3) = settings.matches[j].second(1) * settings.matches[j].first(0);
            A(j, 4) = settings.matches[j].second(1) * settings.matches[j].first(1);
            A(j, 5) = settings.matches[j].second(1);
            A(j, 6) = settings.matches[j].first(0);
            A(j, 7) = settings.matches[j].first(1);
            A(j, 8) = 1;
        }
        svd(A, U, S, V);

        std::vector<eigenvalue_sorter> eigenvalues;
        for (int i = 0; i < S.nc(); i++) {
            eigenvalues.push_back(eigenvalue_sorter(S(i, i), i));
         }

        sort(eigenvalues.begin(), eigenvalues.end());

        for (int i = 0; i < V.nc(); i++) {
            f1(i) = V(i, eigenvalues[1].index);
            f2(i) = V(i, eigenvalues[0].index);
        }

        F1 = reshape(f1, 3, 3);
        F2 = reshape(f2, 3, 3);

        fundamentals = solve(F1, F2);

        double min_error = numeric_limits<double>::max();
        double error;
        for (size_t i = 0; i < fundamentals.size(); i++) {
            error = getQuadraticDistance(settings.matches, getProjectionMatrixes(fundamentals[i]));
            if (error < min_error) {
                min_error = error;
                index_min_error = i;
            }
        }
        matrix<double, 4, 1> X;
        matrix<double, 3, 1> p1_est, p2_est;
        pair<matrix<double, 3, 4>, matrix<double, 3, 4>> projectionMatrixes = getProjectionMatrixes(fundamentals[index_min_error]);

        for (pair<Point, Point> &p : settings.matches) {
            X = DLT_Triangulation(p.first, p.second, projectionMatrixes);
            p1_est = projectionMatrixes.first * X;
            p2_est = projectionMatrixes.second * X;
            p = make_pair<Point, Point>((Point) p1_est, (Point) p2_est);
        }

        itr++;

    }

    return fundamentals[index_min_error];
}

matrix<double> getFundamental8Points(Tracking_settings &settings, const int &iterations = 1) {

    matrix<double> A(settings.matches.size(), 9);
    matrix<double, 4, 1> X;
    matrix<double, 3, 3> F;
    matrix<double, 3, 1> p1_est, p2_est;

    int i = 0;
    while (i != iterations) {
        for (long j = 0; j < A.nr(); j++) {
            A(j, 0) = settings.matches[j].second(0) * settings.matches[j].first(0);
            A(j, 1) = settings.matches[j].second(0) * settings.matches[j].first(1);
            A(j, 2) = settings.matches[j].second(0);
            A(j, 3) = settings.matches[j].second(1) * settings.matches[j].first(0);
            A(j, 4) = settings.matches[j].second(1) * settings.matches[j].first(1);
            A(j, 5) = settings.matches[j].second(1);
            A(j, 6) = settings.matches[j].first(0);
            A(j, 7) = settings.matches[j].first(1);
            A(j, 8) = 1;
        }

        F = reshape(getNullSpaceVector(A, RIGHT_NULL_SPACE_VECTOR), 3, 3);
        force_restriction(F);
        pair<matrix<double, 3, 4>, matrix<double, 3, 4>> projectionMatrixes = getProjectionMatrixes(F);

        for (pair<Point, Point> &p : settings.matches) {
            X = DLT_Triangulation(p.first, p.second, projectionMatrixes);
            p1_est = projectionMatrixes.first * X;
            p2_est = projectionMatrixes.second * X;

            p = make_pair<Point, Point>((Point) p1_est, (Point) p2_est);
        }
        i++;
    }

    return F;

}

int main() {
    auto t_start = std::chrono::high_resolution_clock::now();

    Tracking_settings settings;

    matrix<double, 3, 3> u, s, v, W, Tx, R, E_;

    W = zeros_matrix(W);
    W(0, 0) = 1;
    W(1, 2) = -1;
    W(2, 1) = 1;

//    W(2, 2) = 1;
//    W(0, 1) = 1;
//    W(1, 0) = -1;

    matrix<double, 3, 3> F = getFundamental8Points(settings, 5);

    pair<matrix<double, 3, 4>, matrix<double, 3, 4>> projectionMatrixes = getProjectionMatrixes(F);
    cout << "ERROR: " << getQuadraticDistance(settings.matches, projectionMatrixes) << endl;

    matrix<double, 3, 3> E = trans(settings.instrinsic) * (F * settings.instrinsic);
    cout << "Essential Matrix: " << endl << E / E(2, 2) << endl;

    svd(E, u, s, v);
    cout << s << endl;

    Tx = u * W * s * trans(u);
    cout << "[T]x: " << endl << Tx / Tx(0, 1) << endl;

    R = u * inv(W) * trans(v);
    cout << "R: " << endl << R / R(2, 2)<< endl;

    E_ = Tx * R;
    cout << "Rebuilt Essential Matrix: " << endl << E_ / E_(2, 2)  << endl;

    auto t_elasped = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t_start).count()/1000.0;
    std::cout << "Tracking:                Total ["  << t_elasped << "]   "  << 1.0/t_elasped << " fps\n";

    return 0;
}

