#ifndef GA_CONFIG_H_
#define GA_CONFIG_H_

#include <vector>
#include <iostream>
#include <json/json.h>
struct Weights {
    double w1;
    double w2;
    double w3;
    double w4;
    double x;
};

struct Parameters {
    double crossoverRate;
    double mutationRate;
    int numGenerations;
    int numPopulation;
    Weights weights;
    double overlapFactor;
    bool save;
    int displaySolutions;
    double varianceThreshold;
    int forgiveFactor;
    int saveSolutions;
    std::string output;
    bool use_curv;
};

struct RadarConf {
    int range;
    int antennaHeight;
    float minAngle;
    int target;
    float maxAngle;
    int quantity;
};

struct HomogenousConf {
    int chromosomeLength;
    bool is_case;
    RadarConf radar;
};

struct HeterogenousConf {
    int chromosomeLength;
    bool is_case;
    int target;
    std::vector<RadarConf> radars;
};

struct autoGenConf {
    int minRadars;
    int maxRadars;
    int target;
    bool is_case;
    std::vector<RadarConf> radars;
};

class GaConf {
    public:
    Parameters parameters;
    HomogenousConf homo;
    HeterogenousConf hetero;
    autoGenConf autoG;

    static bool is_available(std::optional <const std::string> path = std::nullopt);
    static GaConf load(std::optional <const std::string> path = std::nullopt);
};

#endif
