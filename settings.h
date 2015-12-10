#ifndef SETTINGS
#define SETTINGS

#include <iostream>
#include <fstream>
#include <vector>
#include <dlib/matrix.h>
#include <dlib/string.h>
#include <Point.h>
#include <cassert>
#include <CImg.h>

using namespace std;
using namespace dlib;
using namespace cimg_library;

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

        void load_extrinsic(string filename, matrix<double, 3, 3> &R, matrix<double, 3, 1> &t) {
            ifstream file(filename);
            string line;
            if (file.is_open()) {
                for (int j = 0; j < 3; j++) {
                    for (int i = 0; i < 3; i++) {
                        getline(file, line);
                        R(j, i) = stod(line);
                    }
                }
                int i = 0;
                while (i != 3) {
                    getline(file, line);
                    t(i) = stod(line);
                    i++;
                }
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
        matrix<double, 3, 3> r_1;
        matrix<double, 3, 1> t_1;

        matrix<double, 3, 3> r_2;
        matrix<double, 3, 1> t_2;

        CImg<uint8_t> image_A;
        CImg<uint8_t> image_B;
        std::vector<pair<Point, Point>> matches;


        Tracking_settings() {
            load_instrinsic("input/intrinsic-temple.dat");
            load_extrinsic("input/extrinsic-temple-R0039.dat", r_1, t_1);
            load_extrinsic("input/extrinsic-temple-R0035.dat", r_2, t_2);
//            load_images("input/templeR0039.png", "input/templeR0035.png");
            load_matches("input/matches.txt");
        }


};


#endif // SETTINGS

