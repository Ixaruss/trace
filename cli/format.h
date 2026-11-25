#pragma once
// cli/output_format.h
//
// Every subcommand that can write tiff/svg/png repeated the same
// if(tiff){...}else if(svg){...}else{...} triage — four times across the
// old main.cc (twice inside viewshed alone, since its sector/non-sector
// branches were otherwise identical). Collapsed to one enum + one small
// dispatch function per engine type.

#include <string>
#include "../operations.h"   // ViewshedEngine, RadarEngine, OutputCrs, BandLayout

enum class OutputFormat { Png, Svg, Tiff };

// tiff/svg/png are the three mutually-exclusive flags a subcommand's
// output-format option_group collected. Priority order (tiff > svg > png)
// matches every existing if(tiff){}else if(svg){}else{} chain exactly —
// this does not change behavior, only names it.
inline OutputFormat resolve_output_format(bool tiff, bool svg, bool /*png*/) {
    if (tiff) return OutputFormat::Tiff;
    if (svg)  return OutputFormat::Svg;
    return OutputFormat::Png;
}

inline const char* extension_for(OutputFormat fmt) {
    switch (fmt) {
        case OutputFormat::Tiff: return ".tif";
        case OutputFormat::Svg:  return ".svg";
        case OutputFormat::Png:  return ".png";
    }
    return ".png";
}

// ViewshedEngine and RadarEngine name their three writers differently
// (writeViewshedTo* vs write*To*) and RadarEngine's take an extra
// RadarCoverageResult, so one dispatch function per engine type rather
// than forcing both through one template — clearer than the generic
// version would be for two call sites.
inline void write_viewshed_output(ViewshedEngine& engine, const std::string& path,
                                  OutputFormat fmt, int reductionPercent,
                                  OutputCrs crs, BandLayout band) {
    switch (fmt) {
        case OutputFormat::Tiff: engine.writeViewshedToGeoTIFF(path, reductionPercent, crs, band); break;
        case OutputFormat::Svg:  engine.writeViewshedToSVG(path); break;
        case OutputFormat::Png:  engine.writeViewshedToPNG(path); break;
    }
}

inline void write_radar_output(RadarEngine& engine, const RadarCoverageResult& res,
                               const std::string& path, OutputFormat fmt,
                               int reductionPercent, OutputCrs crs, BandLayout band) {
    switch (fmt) {
        case OutputFormat::Tiff: engine.writeCoverageToGeoTIFF(res, path, reductionPercent, crs, band); break;
        case OutputFormat::Svg:  engine.writeToSVG(res, path); break;
        case OutputFormat::Png:  engine.writeToPNG(res, path); break;
    }
}
