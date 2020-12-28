/*
 * Copyright (c) 2013-2014, yinqiwen <yinqiwen@gmail.com>
 * Copyright (c) 2014-2020, Matt Stancliff <matt@genges.com>.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

#ifndef GEOHASH_INT_GEOHASH_HELPER_H
#define GEOHASH_INT_GEOHASH_HELPER_H

#include <math.h>
#include "geohash.h"
#include <stdbool.h>
__BEGIN_DECLS

#define GZERO(s)                                                               \
    do {                                                                       \
        (s).bits = (s).step = 0;                                               \
    } while (0)
#define GISZERO(s) (!(s).bits && !(s).step)
#define GISNOTZERO(s) ((s).bits || (s).step)

typedef uint64_t GeoHashFix52Bits;
typedef uint64_t GeoHashVarBits;

typedef struct GeoHashRadius {
    GeoHashBits hash;
    GeoHashArea area;
    GeoHashNeighbors neighbors;
} GeoHashRadius;

int GeoHashBitsComparator(const GeoHashBits *a, const GeoHashBits *b);
uint8_t geohashEstimateStepsByRadius(double rangeMeters);
bool geohashBoundingBox(double lat, double lng, double radiusMeters,
                        double *bounds);
GeoHashRadius geohashGetAreasByRadius(uint8_t coordType, double lat, double lng,
                                      double radiusMeters);
GeoHashRadius geohashGetAreasByRadiusWGS84(double lat, double lng,
                                           double radiusMeters);
GeoHashRadius geohashGetAreasByRadiusMercator(double lat, double lng,
                                              double radiusMeters);
GeoHashFix52Bits geohashAlign52Bits(const GeoHashBits hash);
double geohashGetXMercator(double lng);
double geohashGetYMercator(double lat);
double geohashGetXWGS84(double x);
double geohashGetYWGS84(double y);
bool geohashVerifyCoordinates(uint8_t coordType, double x, double y);
bool geohashGetDistanceIfInRadiusWGS84(double x1, double y1, double x2,
                                       double y2, double radius,
                                       double *distance);
bool geohashGetDistanceIfInRadiusMercator(double x1, double y1, double x2,
                                          double y2, double radius,
                                          double *distance);

double geohashCoordDistanceEarth(double lat1, double lng1, double lat2,
                                 double lng2);
double geohashCoordDistanceEarthFast(double lat1, double lng1, double lat2,
                                     double lng2);

__END_DECLS
#endif /* GEOHASH_INT_GEOHASH_HELPER_H */
