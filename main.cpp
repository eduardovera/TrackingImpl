#include <iostream>
#include <CImg.h>
#include <iostream>
#include <fstream>
#include <settings.h>
#include <dlib/matrix/matrix.h>
#include <limits>
#include <chrono>

using namespace std;

matrix<double> getFundamentalMatrix(Tracking_settings &settings, const int &iterations = 1) {
    std::vector<pair<Point, Point>> points;

    matrix<double, 3, 3> T = identity_matrix<double, 3>();
//    T(0, 0) = settings.image_A.width() + settings.image_A.height();
//    T(0, 2) = settings.image_A.width();
//    T(1, 1) = T(0, 0);
//    T(1, 2) = settings.image_A.height();

    for (pair<Point, Point> match : settings.matches) {
        Point p = T * match.first;
        Point p_ = T * match.second;
        points.push_back(pair<Point, Point>(p, p_));
    }

    matrix<double> A(settings.matches.size(), 9);

    int i = 0;
    while (i != iterations) {
        matrix<double> u, s, vt;
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

        matrix<double, 9, 1> f;

        svd(A, u, s, vt);
        double minS = numeric_limits<double>::max();
        for (int i = 0; i < s.nr(); i++) {
            if (s(i, i) < minS) {
                f = colm(vt, i);
                minS = s(i, i);
            }
        }

        matrix<double, 3, 3> F = reshape(f, 3, 3);

//        F = inv(T) * F * T;

        cout << F << endl;

        double error = 0.0;
        for (pair<Point, Point> p : points) {
            error += p.second * F * trans(p.first);
        }
        cout << "Error: " << error << endl;

        i++;
    }




    return T;

}

int main() {

    auto t_start = std::chrono::high_resolution_clock::now();

    Tracking_settings settings;

    cout << settings.instrinsic << endl;
    cout << settings.image_A.width() << " x " << settings.image_A.height() << endl;
    cout << settings.image_B.width() << " x " << settings.image_B.height() << endl;

    cout << settings.matches[0].first << endl;
    cout << settings.matches[0].second << endl;

    getFundamentalMatrix(settings);

    auto t_elasped = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t_start).count()/1000.0;
    std::cout << "Tracking:                Total ["  << t_elasped << "]   "  << 1.0/t_elasped << " fps\n";

    return 0;
}

