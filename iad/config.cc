#include "config.h"
using namespace std;

bool IadConf::is_available(optional<const string> path) {
    return filesystem::exists(path.value_or("./config/iadconf.json"));
}

void parseHeteroRadarDetails(const Json::Value& node, IadRadarConf& r) {
    r.range = node["range"].asDouble();
    r.antennaHeight = node["antenna"].asDouble();
    r.minAngle = node["min-angle"].asDouble();
    r.maxAngle = node["max-angle"].asDouble();
    r.lat = node["lat"].asDouble();
    r.lon = node["lon"].asDouble();
    r.label = node["label"].asString();
}

void parseRadarArray(const Json::Value& arrayNode, std::vector<IadRadarConf>& list) {
    for (const auto& item : arrayNode) {
        IadRadarConf details;
        parseHeteroRadarDetails(item["radar"], details);
        list.push_back(details);
    }
}


IadConf IadConf::load(optional<const string> path) {
    string p = path.value_or("./config/iadconf.json");

    ifstream f(p);
    if (!f.is_open())
        throw runtime_error("Cannot open config file: " + p);

    Json::Value root;
    Json::CharReaderBuilder reader;
    string errors;
    if (!Json::parseFromStream(reader, f, &root, &errors))
        throw runtime_error("Failed to parse config JSON: " + errors);

    IadConf config;

    config.target = root["target"].asInt();
    config.use_curv = root["use_earth_curvature"].asBool();
    config.output_path = root["output_file"].asString();
    parseRadarArray(root["radars"],config.radars);

    return config;
}
