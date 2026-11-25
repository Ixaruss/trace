#include "config.h"
#include <fstream>
#include <filesystem>
#include <json/reader.h>
using namespace std;

bool Rconf::is_available(optional<const string> path) {
    return filesystem::exists(path.value_or("./config/radar.json"));
}


Rconf from_jsonn(const Json::Value& root) {
    Rconf config;
    config.output = root["output"].asString();
    config.use_curv = root["use_earth_curvature"].asBool();
    config.use_sector = root["sector_flag"].asBool();
    config.loa = root["loa"].asDouble();
    config.roa = root["roa"].asDouble();
    config.use_radar_model = root["use_radar_model"].asBool();
    config.useAMSL = root["useAMSL"].asBool();
    config.range = root["range"].asInt();
    config.antennaHeight = root["antenna"].asInt();
    config.target = root["target"].asInt();
    config.minAngle = root["min-angle"].asDouble();
    config.maxAngle = root["max-angle"].asDouble();
    config.peak_power = root["peak_power"].asInt();
    config.antenna_gain = root["antenna_gain"].asInt();
    config.wavelength = root["wavelength"].asDouble();
    config.crossSection = root["crossSection"].asInt();
    config.detachable_threshold = root["detachable_threshold"].asDouble();
    config.srcLat = root["srcLat"].asDouble();
    config.srcLon = root["srcLon"].asDouble();
    return config;
}

Rconf Rconf::load(optional<const string> path) {
    string p = path.value_or("./config/radar.json");

    ifstream f(p);
    if (!f.is_open())
        throw runtime_error("Cannot open config file: " + p);

    Json::Value root;
    Json::CharReaderBuilder reader;
    string errors;
    if (!Json::parseFromStream(reader, f, &root, &errors))
        throw runtime_error("Failed to parse config JSON: " + errors);

    return from_jsonn(root);
}
