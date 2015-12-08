#include <iostream>
#include <CImg.h>
#include <iostream>
#include <fstream>
#include <settings.h>

using namespace std;

matrix<double> getFundamentalMatrix(Tracking_settings &settings) {
    std::vector<pair<Point, Point>> points;

    matrix<double, 3, 3> T = identity_matrix<double, 3>();
    T(0, 0) = settings.image_A.width() + settings.image_A.height();
    T(0, 2) = settings.image_A.width();
    T(1, 1) = T(0, 0);
    T(1, 2) = settings.image_A.height();

    for (pair<Point, Point> match : settings.matches) {
        Point p = match.first.applyTransform(T);
        Point p_ = match.second.applyTransform(T);
        points.push_back(pair<Point, Point>(p, p_));
    }


    return T;

}

int main() {

    Tracking_settings settings;

    cout << settings.instrinsic << endl;
    cout << settings.image_A.width() << " x " << settings.image_A.height() << endl;
    cout << settings.image_B.width() << " x " << settings.image_B.height() << endl;

    cout << settings.matches[0].first << endl;
    cout << settings.matches[0].second << endl;

    getFundamentalMatrix(settings);

    return 0;
}

