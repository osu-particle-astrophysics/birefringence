#include "antenna_measurements.hh"
#include <fstream>
#include <iostream>
#include <cmath>
#include <string>

void load_antenna_measurements(const Config& cfg, const Geometry& geom, StationData& data) {
    std::string sfile;
    if (geom.WHICHPOL == 0) {
        sfile = "dave_data/day359_pol0.log";
    } else {
        sfile = "dave_data/day359_pol1.log";
    }

    std::ifstream myfile(sfile.c_str());
    std::ifstream davea5file("dave_data/a5_amyformat.txt");

    int NSHOTS = 640;
    double running_rms = 0.0;
    double running_mean = 0.0;

    if (myfile.is_open()) {
        for (int i = 0; i < NSHOTS; i++) {
            int this_station, this_day, this_pol;
            double this_depth, this_snrmax;

            myfile >> this_station >> this_day >> this_pol >> this_depth >> this_snrmax;

            this_snrmax =
                this_snrmax *
                std::sqrt((this_depth - geom.station_depths[this_station - 1]) * (this_depth - geom.station_depths[this_station - 1]) +
                          geom.horizontal_distances[this_station - 1] * geom.horizontal_distances[this_station - 1]) /
                std::sqrt((-1000.0 - geom.station_depths[0]) * (-1000.0 - geom.station_depths[0]) +
                          geom.horizontal_distances[0] * geom.horizontal_distances[0]);

            if (this_snrmax > 0.0 && !(this_station == 5 && this_pol == 0)) {
                data.vdepth_data[this_station - 1].push_back(this_depth);
                data.vsnrmax[this_station - 1].push_back(this_snrmax);
                data.vtotal_distances[this_station - 1].push_back(
                    std::sqrt((this_depth - geom.station_depths[this_station - 1]) * (this_depth - geom.station_depths[this_station - 1]) +
                              geom.horizontal_distances[this_station - 1] * geom.horizontal_distances[this_station - 1]));

                if (i > 2) {
                    running_mean = 0.0;
                    running_rms = 0.0;

                    for (int j = 0; j < 3; j++) {
                        running_mean += data.vsnrmax[this_station - 1][data.vsnrmax[this_station - 1].size() - j - 1];
                    }
                    running_mean /= 3.0;

                    for (int j = 0; j < 3; j++) {
                        double diff = data.vsnrmax[this_station - 1][data.vsnrmax[this_station - 1].size() - j - 1] - running_mean;
                        running_rms += diff * diff;
                    }
                    running_rms = std::sqrt(running_rms / 2.0);
                }

                data.vsnrmax_err[this_station - 1].push_back(running_rms);
                data.vdepth_data_err[this_station - 1].push_back(0.0);
                data.vtotal_distances_err[this_station - 1].push_back(0.0);
            }
        }
    }

    NSHOTS = 33;
    if (geom.WHICHPOL == 0 && davea5file.is_open()) {
        std::cout << "i'm reading dave's a5 file.\n";

        for (int i = 0; i < NSHOTS; i++) {
            int this_station = 5;
            int this_channel;
            double this_depth, this_snrmax;
            std::string stemp;

            davea5file >> stemp >> this_channel >> this_depth >> this_snrmax;

            if (this_snrmax > 0.0) {
                data.vdepth_data[this_station - 1].push_back(this_depth);
                data.vsnrmax[this_station - 1].push_back(this_snrmax);
                data.vtotal_distances[this_station - 1].push_back(
                    std::sqrt((this_depth - geom.station_depths[this_station - 1]) * (this_depth - geom.station_depths[this_station - 1]) +
                              geom.horizontal_distances[this_station - 1] * geom.horizontal_distances[this_station - 1]));

                if (i > 2) {
                    running_mean = 0.0;
                    running_rms = 0.0;

                    for (int j = 0; j < 3; j++) {
                        running_mean += data.vsnrmax[this_station - 1][data.vsnrmax[this_station - 1].size() - j - 1];
                    }
                    running_mean /= 3.0;

                    for (int j = 0; j < 3; j++) {
                        double diff = data.vsnrmax[this_station - 1][data.vsnrmax[this_station - 1].size() - j - 1] - running_mean;
                        running_rms += diff * diff;
                    }
                    running_rms = std::sqrt(running_rms / 2.0);
                }

                data.vsnrmax_err[this_station - 1].push_back(running_rms);
                data.vdepth_data_err[this_station - 1].push_back(0.0);
                data.vtotal_distances_err[this_station - 1].push_back(0.0);
            }
        }
    }
}
