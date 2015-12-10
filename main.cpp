#include <iostream>
#include <CImg.h>
#include <iostream>
#include <fstream>
#include <settings.h>
#include <utils.h>
#include <dlib/matrix.h>
#include <Point.h>
#include <chrono>

using namespace std;

pair<matrix<double, 3, 4>, matrix<double, 3, 4>> getProjectionMatrixesGroundTruth(Tracking_settings &settings) {
    matrix<double, 3, 4> extrinsic, P, P_;
    set_colm(extrinsic, range(0, 2)) = settings.r_1;
    set_colm(extrinsic, 3) = settings.t_1;
    P = settings.instrinsic * extrinsic;

    set_colm(extrinsic, range(0, 2)) = settings.r_2;
    set_colm(extrinsic, 3) = settings.t_2;
    P_ = settings.instrinsic * extrinsic;


    return pair<matrix<double, 3, 4>, matrix<double, 3, 4>>(P, P_);

}

matrix<double, 3, 3> getFundamentalGroundTruth(Tracking_settings &settings) {
    pair<matrix<double, 3, 4>, matrix<double, 3, 4>> projectionMatrixes = getProjectionMatrixesGroundTruth(settings);
    matrix<double, 3, 4> P = projectionMatrixes.first;
    matrix<double, 3, 4> P_ = projectionMatrixes.second;

    matrix<double> C = getNullSpaceVector(P, RIGHT_NULL_SPACE_VECTOR);

    matrix<double, 3, 1> e_ = P_ * C;

    matrix<double, 3, 3> F = crossMatrixForm(e_) * P_ * pinv(P);

    return F;
}


matrix<double> getFundamentalMatrix(Tracking_settings &settings, const int &iterations = 1) {

    matrix<double> A(settings.matches.size(), 9);
    matrix<double, 4, 1> X;
    matrix<double, 3, 3> F;
    matrix<double, 3, 1> p1_est, p2_est;
    double error = 0.0, error2 = 0.0;

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
//        cout << F << endl;
//        F = getFundamentalGroundTruth(settings);
//        cout << getFundamentalGroundTruth(settings) << endl;
//        getchar();
        pair<matrix<double, 3, 4>, matrix<double, 3, 4>> projectionMatrixes = getProjectionMatrixes(F);

        error = 0.0;
        error2 = 0.0;
        for (pair<Point, Point> &p : settings.matches) {
            X = DLT_Triangulation(p.first, p.second, projectionMatrixes);
            p1_est = projectionMatrixes.first * X;
            p2_est = projectionMatrixes.second * X;
            error += getQuadraticDistance(p.first, (Point)p1_est) + getQuadraticDistance(p.second, (Point)p2_est);
            error2 += getGeometricDistance(p.first.innerMatrix(), p.second.innerMatrix(), F);
//            cout << (p.first) << " " << (p.second) << endl;
            p = make_pair<Point, Point>((Point)p1_est, (Point)p2_est);
//            cout << trans(p1_est / p1_est(2)) << " " << trans(p2_est / p2_est(2)) << endl;
//            getchar();
        }
        i++;
    }
    cout << "Error: " << error << endl;
    cout << "Error GEOM: " << error2 << endl;

    return F;

}

int main() {

    auto t_start = std::chrono::high_resolution_clock::now();

    Tracking_settings settings;

//    getFundamentalGroundTruth(settings);

//    cout << settings.instrinsic << endl;
//    cout << settings.image_A.width() << " x " << settings.image_A.height() << endl;
//    cout << settings.image_B.width() << " x " << settings.image_B.height() << endl;

    matrix<double, 3, 3> F = getFundamentalMatrix(settings);
    matrix<double, 3, 1> right_e = getNullSpaceVector(F, RIGHT_NULL_SPACE_VECTOR);
//    cout << right_e / right_e(2) << endl;

    matrix<double, 3, 3> E = trans(settings.instrinsic) * F * settings.instrinsic;
//    cout << F << endl;

    E /= E(2, 2);

    cout << E << endl;

    matrix<double, 3, 3> u, s, v, W, Tx, R, Z;

    W = zeros_matrix(W);
    W(0, 0) = 1;
    W(1, 2) = -1;
    W(2, 1) = 1;

    Z = zeros_matrix(Z);
    Z(0, 1) = 1;
    Z(1, 0) = -1;

    svd(E, u, s, v);

    cout << s << endl;

    Tx = u * W * s * trans(u);
//    cout << Tx << endl;

    R = u * inv(W) * trans(v);
//    R = u * trans(W) * vt;
//    cout << R << endl;

    Z = Tx * R;

    Z /= Z(2, 2);

    cout << Z << endl;

    auto t_elasped = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t_start).count()/1000.0;
    std::cout << "Tracking:                Total ["  << t_elasped << "]   "  << 1.0/t_elasped << " fps\n";

    return 0;
}

