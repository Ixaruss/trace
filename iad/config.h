#ifndef CONFIG_IAD_
#define CONFIG_IAD_

#include <vector>
#include <string>
#include <optional>
#include <json/json.h>
#include <filesystem>
#include <fstream>

struct IadRadarConf {
    int range;
    int antennaHeight;
    float minAngle;
    float maxAngle;
    double lat;
    double lon;
    std::string label;
};

class IadConf {
    public:
    std::vector<IadRadarConf> radars;
    int target;
    bool use_curv;
    std::string output_path;

    static bool is_available(std::optional <const std::string> path = std::nullopt);
    static IadConf load(std::optional <const std::string> path = std::nullopt);

};
#endif
