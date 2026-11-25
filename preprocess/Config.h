#ifndef CONFIG_H_
#define CONFIG_H_

#include <gdal.h>
#include <iostream>
#include <string>
#include <optional>
#include <iostream>
#include <json/json.h>
#include <map>


//json parser enums and class
enum Target {
    BUILD,
    SERVER
};
struct Input {
    std::string road;
    std::string rail;
    std::string water_area;
    std::string water_lines;
    std::string ib;
    std::string elevation;
};

struct Output {
    std::string matrix_file;
    std::string elevation_file;
    std::string feasible;
    std::string slope_file;
};

enum Layer_vals {
    ROADS,
    RAILS,
    WATER_LINES,
    WATER_AREAS,
    IB,
    SLOPE,
    ALT
};
struct Case {
    bool roads;
    bool rails;
    bool waterbodies;
    bool waterlines;
    bool ib;
    bool slope;
};
struct Grid_Conf {
  int resolution = 0;
  float origin_x = 0;
  float origin_y = 0;
  double min_x = 0, max_x = 0;
  double min_y = 0, max_y = 0;
  int   rows = 0;
  int   cols = 0;
  // 0 = single flat bins (elev.bin/slope.bin/map.bin); >0 = tiled folders.
  // Must default to 0: if the JSON has no "block_size" key this used to be
  // whatever garbage was on the stack, i.e. tiling switched itself on or off.
  int block_size = 0;

};

struct Distance {
    int from_road;
    int from_rail;
    int from_water;
    int from_ib;
};

class Config {
    public: Target target;
    public: Input input;
    public: Output output;
    public: Grid_Conf grid;
    public: Case fcase;
    public: Distance distance;
    public: bool agl_flag;

    public: std::vector<std::pair<double,double>> area_op; // can be a polygon or retangle, if not rectangle find highest breath and hightest height and assign the values to above variable to form a rectangle;

    static Config load(std::optional <const std::string> path = std::nullopt);
    void save(std::string path = "./config/config.json") const;
    static Config from_json(const Json::Value& j);
    Json::Value to_json() const;
    static bool is_available(std::optional <const std::string> path = std::nullopt);
};

//helper functions
std::string _extract_part(const std::string& file_path, unsigned short with_extension = 0);
#endif
