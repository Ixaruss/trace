#include "config.h"
#include <fstream>
#include <filesystem>
#include <stdexcept>
using namespace std;

bool GaConf::is_available(optional<const string> path) {
    return filesystem::exists(path.value_or("./gaconf.json"));
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

GaConf GaConf::from_json(const Json::Value& root) {
    GaConf config;

    config.parameters.crossoverRate = root["crossover-rate"].asDouble();
    config.mutation_rate = root["mutation-rate"].asDouble();
    config.num_generation = root["num-generation"].asInt();
    config.num_population = root["num-population"].asInt();
    config.display_solutions = root["display_solutions"].asInt();
    config.save_solutions = root["save_solutions"].asInt();
    config.output = root["output"].asString();
    config.scatter_output = root["scatter-output"].asString();

    // --- Weights Sub-Object ---
    Json::Value wNode = root["weights"];
    config.weights.w1 = wNode["w1"].asDouble();
    config.weights.w2 = wNode["w2"].asDouble();
    config.weights.w3 = wNode["w3"].asDouble();
    config.weights.w4 = wNode["w4"].asDouble();
    config.weights.x = wNode["x"].asInt();

    // --- Homogenous Settings ---
    Json::Value homoNode = root["homogenous"];
    config.homogenous.chromosome_length = homoNode["chromosome-lentgh"].asInt(); // spelling matches JSON
    config.homogenous.is_case = homoNode["case"].asBool();

    Json::Value hRadarNode = homoNode["radar"];
    config.homogenous.radar.range = hRadarNode["range"].asInt();
    config.homogenous.radar.antenna = hRadarNode["antenna"].asInt();
    config.homogenous.radar.target = hRadarNode["target"].asInt();
    config.homogenous.radar.min_angle = hRadarNode["min-angle"].asDouble();
    config.homogenous.radar.max_angle = hRadarNode["max-angle"].asDouble();

    // --- Heterogenous Settings ---
    Json::Value heteroNode = root["heterogenous"];
    config.heterogenous.chromosome_length = heteroNode["chromosome-lentgh"].asInt();
    config.heterogenous.is_case = heteroNode["case"].asBool();
    config.heterogenous.target = heteroNode["target"].asInt();
    config.heterogenous.radars = parseRadarArray(heteroNode["radars"]);

    // --- Auto Settings ---
    Json::Value autoNode = root["auto"];
    config.auto_settings.min_radar = autoNode["min-radar"].asInt();
    config.auto_settings.max_radar = autoNode["max-radar"].asInt();
    config.auto_settings.radars = parseRadarArray(autoNode["radars"]);



        return conf;
}
