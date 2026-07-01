#include "config.h"
#include <fstream>
#include <filesystem>
#include <stdexcept>
using namespace std;

bool GaConf::is_available(optional<const string> path) {
    return filesystem::exists(path.value_or("./gaconf.json"));
}

void parseHeteroRadarDetails(const Json::Value& node, RadarConf& r) {
    r.range = node["range"].asDouble();
    r.antennaHeight = node["antenna"].asDouble();
    r.minAngle = node["min-angle"].asDouble();
    r.maxAngle = node["max-angle"].asDouble();
    r.quantity = node["quantity"].asInt();
}

void parseRadarArray(const Json::Value& arrayNode, std::vector<RadarConf>& list) {
    for (const auto& item : arrayNode) {
        RadarConf details;
        parseHeteroRadarDetails(item["radar"], details);
        list.push_back(details);
    }
}

GaConf from_json(const Json::Value& root) {
    GaConf config;

    config.parameters.crossoverRate = root["crossover-rate"].asDouble();
    config.parameters.mutationRate = root["mutation-rate"].asDouble();
    config.parameters.numGenerations = root["num-generation"].asInt();
    config.parameters.numPopulation = root["num-population"].asInt();
    config.parameters.displaySolutions = root["display_solutions"].asInt();
    config.parameters.saveSolutions = root["save_solutions"].asInt();
    config.parameters.save = root["save"].asBool();
    config.parameters.output = root["output"].asString();
    config.parameters.scatterOutput = root["scatter-output"].asString();

    // --- Weights Sub-Object ---
    Json::Value wNode = root["weights"];
    config.parameters.weights.w1 = wNode["w1"].asDouble();
    config.parameters.weights.w2 = wNode["w2"].asDouble();
    config.parameters.weights.w3 = wNode["w3"].asDouble();
    config.parameters.weights.w4 = wNode["w4"].asDouble();
    config.parameters.weights.x = wNode["x"].asInt();

    // --- Homogenous Settings ---
    Json::Value homoNode = root["homogenous"];
    config.homo.chromosomeLength = homoNode["chromosome-length"].asInt(); // spelling matches JSON
    config.homo.is_case = homoNode["case"].asBool();

    Json::Value hRadarNode = homoNode["radar"];
    config.homo.radar.range = hRadarNode["range"].asInt();
    config.homo.radar.antennaHeight = hRadarNode["antenna"].asInt();
    config.homo.radar.target = hRadarNode["target"].asInt();
    config.homo.radar.minAngle = hRadarNode["min-angle"].asDouble();
    config.homo.radar.maxAngle = hRadarNode["max-angle"].asDouble();

    // --- Heterogenous Settings ---
    Json::Value heteroNode = root["heterogenous"];
    config.hetero.chromosomeLength = heteroNode["chromosome-length"].asInt();
    config.hetero.is_case = heteroNode["case"].asBool();
    config.hetero.target = heteroNode["target"].asInt();
     parseRadarArray(heteroNode["radars"],config.hetero.radars);
     for (size_t i = 0; i < config.hetero.radars.size(); ++i) {
         config.hetero.radars[i].target =config.hetero.target;
     }
    // --- Auto Settings ---
    Json::Value autoNode = root["auto"];
    config.autoG.minRadars = autoNode["min-radar"].asInt();
    config.autoG.maxRadars = autoNode["max-radar"].asInt();
    config.autoG.is_case = autoNode["case"].asBool();
    parseRadarArray(autoNode["radars"],config.autoG.radars);


    return config;
}

GaConf GaConf::load(optional<const string> path) {
    string p = path.value_or("./gaconf.json");

    ifstream f(p);
    if (!f.is_open())
        throw runtime_error("Cannot open config file: " + p);

    Json::Value root;
    Json::CharReaderBuilder reader;
    string errors;
    if (!Json::parseFromStream(reader, f, &root, &errors))
        throw runtime_error("Failed to parse config JSON: " + errors);

    return from_json(root);
}
