#ifndef R_CONFIG_H_
#define R_CONFIG_H_

#include <json/writer.h>
#include <vector>
#include <iostream>

class Rconf {
    public:
        std::string output;
        bool use_sector;
        double loa;
        double roa;
        bool use_curv;
        bool use_radar_model;
        bool useAMSL;
        int range;
        int antennaHeight;
        float minAngle;
        int target;
        float maxAngle;
        int peak_power; int antenna_gain; double wavelength;
        int crossSection; double detachable_threshold;
        double srcLat;
        double srcLon;

    static bool is_available(std::optional <const std::string> path = std::nullopt);
    static Rconf load(std::optional <const std::string> path = std::nullopt);

};

#endif
