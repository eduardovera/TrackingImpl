#ifndef SETTINGS
#define SETTINGS

#include <iostream>
#include <fstream>
#include <vector>
#include <dlib/matrix.h>
#include <dlib/string.h>
#include <cassert>
#include <CImg.h>

using namespace std;
using namespace dlib;
using namespace cimg_library;

class Point : public matrix<double, 1, 3> {

    public:

        Point() {

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

        Point& operator * (const matrix<double, 3, 3> &T, const Point &p) {
            matrix<double, 1, 3> p_temp;
            p_temp(0) = p(0);
            p_temp(1) = p(1);
            p_temp(2) = p(2);

            return (T * p_temp);
        }

};

class Tracking_settings {

    private:
        void load_instrinsic(string filename) {
            this->instrinsic = identity_matrix<double, 3>();
            ifstream file(filename);
            string line;
            if (file.is_open()) {
                getline(file, line);
                instrinsic(0) = stod(line);
                getline(file, line);
                instrinsic(4) = stod(line);
                getline(file, line);
                instrinsic(2) = stod(line);
                getline(file, line);
                instrinsic(5) = stod(line);
            }

        }

        void load_images(string filenameA, string filenameB) {
            image_A.load(filenameA.c_str());
            image_B.load(filenameB.c_str());
        }

        void load_matches(string filename) {
            ifstream file(filename);
            std::vector<string> tokenized_line;
            string line;
            if (file.is_open()) {
                getline(file, line);
                int num_matches = stoi(line);
                while (num_matches != 0) {
                    getline(file, line);
                    tokenized_line = dlib::split(line, " ");
                    matches.push_back(make_pair<Point, Point>(Point(stod(tokenized_line[0]), stod(tokenized_line[1])), Point(stod(tokenized_line[2]), stod(tokenized_line[3]))));
                    num_matches--;
                }
            }
        }

    public:
        matrix<double> instrinsic;
        CImg<uint8_t> image_A;
        CImg<uint8_t> image_B;
        std::vector<pair<Point, Point>> matches;

        Tracking_settings() {
            load_instrinsic("input/intrinsic.dat");
            load_images("input/templeR0035.png", "input/templeR0039.png");
            load_matches("input/matches.txt");
        }


};


#endif // SETTINGS

