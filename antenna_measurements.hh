#pragma once

#include "config.hh"
#include "geometry.hh"
#include "station_data.hh"

void load_antenna_measurements(const Config& cfg, const Geometry& geom, StationData& data);

