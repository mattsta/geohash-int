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

/* This started out as a C++ to C conversion from the ardb project.
 * This file originated as:
 * https://github.com/yinqiwen/ardb/blob/d42503/src/geo/geohash_helper.cpp
 * It has been through many modifications since.
 */

#include "geohashHelper.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923132169163975144
#endif

#define D_R (M_PI / 180.0)
#define R_MAJOR 6378137.0
#define R_MINOR 6356752.3142
#define RATIO (R_MINOR / R_MAJOR)
#define ECCENT (sqrt(1.0 - (RATIO * RATIO)))
#define COM (0.5 * ECCENT)

/* degrees to radians constant pi/1080 */
const double DEG_TO_RAD = 0.017453292519943295769236907684886;

/* Earth's quatratic mean radius for WGS-84 */
const double EARTH_RADIUS_IN_METERS = 6372797.560856;

const double MERCATOR_MAX = 20037726.37;
const double MERCATOR_MIN = -20037726.37;

static inline double degRad(double ang) {
    return ang * D_R;
}
static inline double radDeg(double ang) {
    return ang / D_R;
}

double mercatorY(double lat) {
    lat = fmin(89.5, fmax(lat, -89.5));
    double phi = degRad(lat);
    double sinphi = sin(phi);
    double con = ECCENT * sinphi;
    con = pow((1.0 - con) / (1.0 + con), COM);
    double ts = tan(0.5 * (M_PI * 0.5 - phi)) / con;
    return 0 - R_MAJOR * log(ts);
}

double mercatorX(double lon) {
    return R_MAJOR * degRad(lon);
}

double mercLon(double x) {
    return radDeg(x) / R_MAJOR;
}

double mercLat(double y) {
    double ts = exp(-y / R_MAJOR);
    double phi = M_PI_2 - 2 * atan(ts);
    double dphi = 1.0;
    int i;
    for (i = 0; fabs(dphi) > 0.000000001 && i < 15; i++) {
        double con = ECCENT * sin(phi);
        dphi =
            M_PI_2 - 2 * atan(ts * pow((1.0 - con) / (1.0 + con), COM)) - phi;
        phi += dphi;
    }

    return radDeg(phi);
}

/* You must *ONLY* estimate steps when you are encoding.
 * If you are decoding, always decode to GEO_STEP_MAX (26). */
uint8_t geohashEstimateStepsByRadius(double rangeMeters) {
    uint8_t step = 1;
    while (rangeMeters > 0 && rangeMeters < MERCATOR_MAX) {
        rangeMeters *= 2;
        step++;
    }

    step--;
    if (!step) {
        step = GEO_STEP_MAX; /* if range = 0, give max resolution */
    }

    return step > GEO_STEP_MAX ? GEO_STEP_MAX : step;
}

double geohashGetXWGS84(double x) {
    return mercLon(x);
}

double geohashGetYWGS84(double y) {
    return mercLat(y);
}

double geohashGetXMercator(double lng) {
    if (lng > 180 || lng < -180) {
        return lng;
    }

    return mercatorX(lng);
}

double geohashGetYMercator(double lat) {
    if (lat > 90 || lat < -90) {
        return lat;
    }

    return mercatorY(lat);
}

int geohashBitsComparator(const GeoHashBits *a, const GeoHashBits *b) {
    /* If step not equal, compare on step.  else, compare on bits. */
    return a->step != b->step ? a->step - b->step : a->bits - b->bits;
}

bool geohashBoundingBox(double lat, double lng, double radiusMeters,
                        double *bounds) {
    if (!bounds) {
        return false;
    }

    double latr, lonr;
    latr = degRad(lat);
    lonr = degRad(lng);

    double distance = radiusMeters / EARTH_RADIUS_IN_METERS;
    double minLat = latr - distance;
    double maxLat = latr + distance;

    /* Note: we're being lazy and not accounting for coordinates near poles */
    double minLng, maxLng;
    double differenceLng = asin(sin(distance) / cos(latr));
    minLng = lonr - differenceLng;
    maxLng = lonr + differenceLng;

    bounds[0] = radDeg(minLat);
    bounds[1] = radDeg(minLng);
    bounds[2] = radDeg(maxLat);
    bounds[3] = radDeg(maxLng);

    return true;
}

GeoHashRadius geohashGetAreasByRadius(uint8_t coordType, double lat, double lng,
                                      double radiusMeters) {
    GeoHashRange range = {{0}};
    GeoHashRadius radius = {{0}};
    GeoHashBits hash = {0};
    GeoHashNeighbors neighbors = {{0}};
    GeoHashArea area = {{0}};
    double deltaLng, deltaLat;
    double minLat, maxLat, minLon, maxLon;
    int steps;

    if (coordType == GEO_WGS84_TYPE) {
        double bounds[4];
        geohashBoundingBox(lat, lng, radiusMeters, bounds);
        minLat = bounds[0];
        minLon = bounds[1];
        maxLat = bounds[2];
        maxLon = bounds[3];
    } else {
        deltaLat = deltaLng = radiusMeters;
        minLat = lat - deltaLat;
        maxLat = lat + deltaLat;
        minLon = lng - deltaLng;
        maxLon = lng + deltaLng;
    }

    steps = geohashEstimateStepsByRadius(radiusMeters);

    geohashGetCoordRange(coordType, &range);
    geohashEncode(&range, lat, lng, steps, &hash);
    geohashNeighbors(&hash, &neighbors);
    geohashDecode(&range, hash, &area);

    if (area.range.lat.min < minLat) {
        GZERO(neighbors.south);
        GZERO(neighbors.southWest);
        GZERO(neighbors.southEast);
    }

    if (area.range.lat.max > maxLat) {
        GZERO(neighbors.north);
        GZERO(neighbors.northEast);
        GZERO(neighbors.northWest);
    }

    if (area.range.lng.min < minLon) {
        GZERO(neighbors.west);
        GZERO(neighbors.southWest);
        GZERO(neighbors.northWest);
    }

    if (area.range.lng.max > maxLon) {
        GZERO(neighbors.east);
        GZERO(neighbors.southEast);
        GZERO(neighbors.northEast);
    }

    radius.hash = hash;
    radius.neighbors = neighbors;
    radius.area = area;
    return radius;
}

GeoHashRadius geohashGetAreasByRadiusWGS84(double lat, double lng,
                                           double radiusMeters) {
    return geohashGetAreasByRadius(GEO_WGS84_TYPE, lat, lng, radiusMeters);
}

GeoHashRadius geohashGetAreasByRadiusMercator(double lat, double lng,
                                              double radiusMeters) {
    return geohashGetAreasByRadius(GEO_MERCATOR_TYPE, lat, lng, radiusMeters);
}

GeoHashFix52Bits geohashAlign52Bits(const GeoHashBits hash) {
    uint64_t bits = hash.bits;
    bits <<= (GEOHASH_BIT_WIDTH - hash.step * 2);
    return bits;
}

/* calculate distance using haversin great circle distance formula */
double geohashCoordDistanceEarth(double lat1, double lng1, double lat2,
                                 double lng2) {
    double lat1Radins = degRad(lat1);
    double lng1Radins = degRad(lng1);
    double lat2Radins = degRad(lat2);
    double lng2Radins = degRad(lng2);
    double u = sin(0.5 * (lat2Radins - lat1Radins));
    double v = sin(0.5 * (lng2Radins - lng1Radins));
    return 2.0 * EARTH_RADIUS_IN_METERS *
           asin(sqrt(u * u + cos(lat1Radins) * cos(lat2Radins) * v * v));
}

/* calculate distance using linear(ish) distance metrics. */
double geohashCoordDistanceEarthFast(double lat1, double lng1, double lat2,
                                     double lng2) {
    /* Also see: https://gis.stackexchange.com/questions/58653/
     *      and: http://www.thekompf.com/gps/distcalc.html */
    double l = (lat2 + lat1) * 0.5 * DEG_TO_RAD;
    double dx = 111.3 * cos(l) * (lng2 - lng1);
    double dy = 111.3 * (lat2 - lat1);

    /* If you want to get *really* fancy, replace sqrt() with something from:
     * wikipedia article:
     *  Methods_of_computing_square_roots
     *    #Approximations_that_depend_on_the_floating_point_representation */
    return sqrt(dx * dx + dy * dy);
}

static bool _geohashGetDistanceIfInRadius(uint8_t coordType, double x1,
                                          double y1, double x2, double y2,
                                          double radius, double *distance,
                                          double (*distCalc)(double, double,
                                                             double, double)) {
    if (coordType == GEO_WGS84_TYPE) {
        *distance = distCalc(y1, x1, y2, x2);
    } else {
        /* else, we have a square grid, use p'theorem */
        *distance = geohashCoordDistanceEarthFast(y1, x1, y2, x2);
    }

    return (*distance <= radius);
}

bool geohashGetDistanceIfInRadiusWGS84(double x1, double y1, double x2,
                                       double y2, double radius,
                                       double *distance) {
    return _geohashGetDistanceIfInRadius(GEO_WGS84_TYPE, x1, y1, x2, y2, radius,
                                         distance, geohashCoordDistanceEarth);
}

bool geohashGetDistanceIfInRadiusWGS84Fast(double x1, double y1, double x2,
                                           double y2, double radius,
                                           double *distance) {
    return _geohashGetDistanceIfInRadius(GEO_WGS84_TYPE, x1, y1, x2, y2, radius,
                                         distance,
                                         geohashCoordDistanceEarthFast);
}

bool geohashGetDistanceIfInRadiusMercator(double x1, double y1, double x2,
                                          double y2, double radius,
                                          double *distance) {
    return _geohashGetDistanceIfInRadius(GEO_MERCATOR_TYPE, x1, y1, x2, y2,
                                         radius, distance, NULL);
}

bool geohashVerifyCoordinates(uint8_t coordType, double x, double y) {
    GeoHashRange range;
    geohashGetCoordRange(coordType, &range);

    if (x < range.lng.min || x > range.lng.max || y < range.lat.min ||
        y > range.lat.max) {
        return false;
    }

    return true;
}
