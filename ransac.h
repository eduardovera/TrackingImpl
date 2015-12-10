#ifndef RANSAC
#define RANSAC

#include <dlib/matrix.h>
#include <Point.h>
#include <iostream>
#include <time.h>

using namespace std;
using namespace cv;

class Ransac {

private:
    static void randomIndexes(const int maxIndex, vector<int> &pointIndexes, const vector<pair<Point, Point>> &params) {
        pointIndexes.clear();
        srand(time(NULL));
        while (pointIndexes.size() != 8) {
            int r = rand() % maxIndex;
            if (!contains(r, pointIndexes)) {
                pointIndexes.push_back(r);
            }
        }
    }

    static bool contains(const int &value, const vector<int> &v) {
        for (int i = 0; i < v.size(); i++) {
            if (v[i] == value) {
                return true;
            }
        }
        return false;
    }

    static void computeInliers(const vector<pair<Point, Point>> &params, const matrix<double, 3, 3> &F, vector<pair<Point, Point>> &inliers_per_iteration) {
        const double T = 0.97; //sqrt[3.84*(0.5)²];

        double distance = 0.0;
        inliers_per_iteration.clear();
        for (int i = 0; i < params.size(); i++) {
            distance = getGeometricDistance(params[i], F, NULL);
            if (distance < T) {
                inliers_per_iteration.push_back(params[i]);
            }
        }
    }

public:

    static void run(const Mat1i &imageA, const Mat1i &imageB, const vector<pair<Mat1d, Mat1d>> &params, vector<pair<Mat1d, Mat1d>> &inliers, const int &distance_type) {
        vector<pair<Mat1d, Mat1d>> inliers_per_iteration;
        vector<Mat1d> F;
        int max_inliers_size = 0;
        int N = INT_MAX;
        double distance = 0;
        int k = 1;
        vector<pair<Mat1d, Mat1d>> selected_params;
        vector<int> indexes;
        inliers.clear();
        while (N > k) {
            selected_params.clear();
            randomIndexes(params.size(), indexes, params);

            for (int i = 0; i < indexes.size(); i++) {
                selected_params.push_back(params[indexes[i]]);
            }

            fundamentalMatrix(selected_params, imageA, imageB, F, true, "7-POINTS");

            for (int i = 0; i < F.size(); i++) {
                computeInliers(params, F[i], inliers_per_iteration);
                if (inliers_per_iteration.size() > max_inliers_size) {
                    inliers = inliers_per_iteration;
                    max_inliers_size = inliers.size();
                    double e = 1.0 - (static_cast<double>(static_cast<double>(inliers_per_iteration.size()) / params.size()));
                    N = log(1.0 - 0.95) / (log(1 - pow(1 - e, 7)));
                }
            }
            k++;
        }
        /* GERANDO ARQUIVO COM INLIERS PARA GARANTIR OS MESMOS INLIERS PARA TODOS OS MÉTODOS DURANTE OS TESTES */
        string fname = "inliers.dat";
        ofstream f;
        f.open(fname);
        f << inliers.size() << endl;
        for (int i = 0; i < inliers.size(); i++) {
            f << inliers[i].first(0) << " " << inliers[i].first(1) << " " << inliers[i].second(0) << " " << inliers[i].second(1) << endl;
        }
        f.close();
    }
};

#endif // RANSAC

