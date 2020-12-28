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

#include "geohash.h"

/**
 * Hashing works like this:
 * Divide the world into 4 buckets.  Label each one as such:
 *  -----------------
 *  |       |       |
 *  |       |       |
 *  | 0,1   | 1,1   |
 *  -----------------
 *  |       |       |
 *  |       |       |
 *  | 0,0   | 1,0   |
 *  -----------------
 */

/* These are constraints from EPSG:900913 / EPSG:3785 / OSGEO:41001 */
/* We can't geocode at the north/south pole. */
static const GeoHashRange _GeoRangeWGS84 = {
    .lat = {.min = -85.05112878, .max = 85.05112878},
    .lng = {.min = -180.0, .max = 180.0}};

static const GeoHashRange _GeoRangeMercator = {
    .lat = {.min = -20037726.37, .max = 20037726.37},
    .lng = {.min = -20037726.37, .max = 20037726.37}};

void geohashRangeSetWGS84(GeoHashRange *range) {
    *range = _GeoRangeWGS84;
}
const GeoHashRange *geohashRangeGetWGS84(void) {
    return &_GeoRangeWGS84;
}

void geohashRangeSetMercator(GeoHashRange *range) {
    *range = _GeoRangeMercator;
}
const GeoHashRange *geohashRangeGetMercator(void) {
    return &_GeoRangeMercator;
}

bool geohashGetCoordRange(GeoType coordType, GeoHashRange *range) {
    switch (coordType) {
    case GEO_WGS84_TYPE:
        geohashRangeSetWGS84(range);
        break;
    case GEO_MERCATOR_TYPE:
        geohashRangeSetMercator(range);
        break;
    default:
        return false;
    }

    return true;
}

/* If we're compiled with -mbmi2 and running on a Haswell platform (2013+),
 * we can compute the z-order curve in ~2.5 CPU cycles per dimension.
 * Below, both intrinsic {de,}interleave64 functions take less than 6 CPU
 * cycles (as compared with ~18 CPU cycles for the non-intrinsic version). */
#ifdef __BMI2__
#include <immintrin.h>
static const bool _geohashBuiltForBMI2 = true;
static const uint64_t __X = 0x5555555555555555ULL; /* 0x5 = 0101 */
static const uint64_t __Y = 0xAAAAAAAAAAAAAAAAULL; /* 0xA = 1010 */
static uint64_t interleave64(uint32_t x, uint32_t y) {
    /* pdep(a, b) = parallel deposit from 'a' according to selector mask 'b' */
    return _pdep_u64(x, __X) | _pdep_u64(y, __Y);
}

static void deinterleave64(uint64_t combined, uint32_t *x, uint32_t *y) {
    /* pext(a, b) = parallel extract from 'a' according to selector mask 'b' */
    *x = (uint32_t)_pext_u64(combined, __X);
    *y = (uint32_t)_pext_u64(combined, __Y);
}
#else
static const bool _geohashBuiltForBMI2 = false;
/* Interleave lower bits of x and y, so the bits of x
 * are in the even positions and bits from y in the odd;
 * x and y must initially be less than 2^32 (65536).
 * From:  https://graphics.stanford.edu/~seander/bithacks.html#InterleaveBMN
 */
static const uint64_t B[] = {0x5555555555555555ULL, 0x3333333333333333ULL,
                             0x0F0F0F0F0F0F0F0FULL, 0x00FF00FF00FF00FFULL,
                             0x0000FFFF0000FFFFULL, 0x00000000FFFFFFFFULL};
static const uint8_t S[] = {0, 1, 2, 4, 8, 16};
static inline uint64_t interleave64(uint32_t X, uint32_t Y) {
    uint64_t x = X;
    uint64_t y = Y;

    x = (x | (x << S[5])) & B[4];
    y = (y | (y << S[5])) & B[4];

    x = (x | (x << S[4])) & B[3];
    y = (y | (y << S[4])) & B[3];

    x = (x | (x << S[3])) & B[2];
    y = (y | (y << S[3])) & B[2];

    x = (x | (x << S[2])) & B[1];
    y = (y | (y << S[2])) & B[1];

    x = (x | (x << S[1])) & B[0];
    y = (y | (y << S[1])) & B[0];

    return x | (y << 1);
}

/* reverse the interleave process
 * derived from http://stackoverflow.com/questions/4909263
 */
static inline void deinterleave64(uint64_t interleaved, uint32_t *X,
                                  uint32_t *Y) {
    uint64_t x = interleaved;
    uint64_t y = interleaved >> 1;

    x = (x | (x >> S[0])) & B[0];
    y = (y | (y >> S[0])) & B[0];

    x = (x | (x >> S[1])) & B[1];
    y = (y | (y >> S[1])) & B[1];

    x = (x | (x >> S[2])) & B[2];
    y = (y | (y >> S[2])) & B[2];

    x = (x | (x >> S[3])) & B[3];
    y = (y | (y >> S[3])) & B[3];

    x = (x | (x >> S[4])) & B[4];
    y = (y | (y >> S[4])) & B[4];

    x = (x | (x >> S[5])) & B[5];
    y = (y | (y >> S[5])) & B[5];

    *X = x;
    *Y = y;
}
#endif

bool geohashEncode(const GeoHashRange *range, double lat, double lng,
                   uint8_t step, GeoHashBits *hash) {
    if (NULL == hash || step > GEO_STEP_MAX || !step || RANGEISZERO(range)) {
        return false;
    }

    hash->bits = 0;
    hash->step = step;

    if (lat < range->lat.min || lat > range->lat.max || lng < range->lng.min ||
        lng > range->lng.max) {
        return false;
    }

    double latOffset =
        (lat - range->lat.min) / (range->lat.max - range->lat.min);
    double lngOffset =
        (lng - range->lng.min) / (range->lng.max - range->lng.min);

    /* convert to fixed point based on the step size */
    latOffset *= (1 << step);
    lngOffset *= (1 << step);

    uint32_t latOffsetInt = (uint32_t)latOffset;
    uint32_t lngOffsetInt = (uint32_t)lngOffset;

    /* Lat is Y and Lng is X.
     * lat/lng pairs are (y, x).
     * We encode our hash as (x, y) and that's (lng, lat) */
    hash->bits = interleave64(lngOffsetInt, latOffsetInt);
    return true;
}

bool geohashEncodeType(uint8_t coordType, double lat, double lng, uint8_t step,
                       GeoHashBits *hash) {
    GeoHashRange range = {{0}};
    geohashGetCoordRange(coordType, &range);
    return geohashEncode(&range, lat, lng, step, hash);
}

bool geohashEncodeWGS84(double lat, double lng, uint8_t step,
                        GeoHashBits *hash) {
    return geohashEncodeType(GEO_WGS84_TYPE, lat, lng, step, hash);
}

bool geohashEncodeMercator(double lat, double lng, uint8_t step,
                           GeoHashBits *hash) {
    return geohashEncodeType(GEO_MERCATOR_TYPE, lat, lng, step, hash);
}

bool geohashDecode(const GeoHashRange *range, const GeoHashBits hash,
                   GeoHashArea *area) {
    if (HASHISZERO(hash) || NULL == area || RANGEISZERO(range)) {
        return false;
    }

    area->hash = hash;
    uint8_t step = hash.step;

    uint32_t gridLng; /* X, lng, east-west */
    uint32_t gridLat; /* Y, lat, north-south */
    deinterleave64(hash.bits, &gridLng, &gridLat);

    double latScale = range->lat.max - range->lat.min;
    double lngScale = range->lng.max - range->lng.min;

    /* - divide by 2^step.
     * - for 0-1 coordinate
     *    - multiply by scale
     *    - add to the min to get the absolute coordinate. */
    area->range.lat.min =
        range->lat.min + (gridLat * 1.0 / (1ULL << step)) * latScale;
    area->range.lat.max =
        range->lat.min + ((gridLat + 1) * 1.0 / (1ULL << step)) * latScale;
    area->range.lng.min =
        range->lng.min + (gridLng * 1.0 / (1ULL << step)) * lngScale;
    area->range.lng.max =
        range->lng.min + ((gridLng + 1) * 1.0 / (1ULL << step)) * lngScale;

    return true;
}

bool geohashDecodeType(uint8_t coordType, const GeoHashBits hash,
                       GeoHashArea *area) {
    GeoHashRange r = {{0}};
    geohashGetCoordRange(coordType, &r);
    return geohashDecode(&r, hash, area);
}

bool geohashDecodeWGS84(const GeoHashBits hash, GeoHashArea *area) {
    return geohashDecodeType(GEO_WGS84_TYPE, hash, area);
}

bool geohashDecodeMercator(const GeoHashBits hash, GeoHashArea *area) {
    return geohashDecodeType(GEO_MERCATOR_TYPE, hash, area);
}

bool geohashDecodeAreaToLatLng(const GeoHashArea *area, double *lat,
                               double *lng) {
    if (!lat || !lng) {
        return false;
    }

    double y = (0.5 * (area->range.lat.min + area->range.lat.max));
    double x = (0.5 * (area->range.lng.min + area->range.lng.max));

    *lat = y;
    *lng = x;
    return true;
}

bool geohashDecodeToLatLngType(uint8_t coordType, const GeoHashBits hash,
                               double *lat, double *lng) {
    GeoHashArea area = {{0}};
    if (!lat || !lng || !geohashDecodeType(coordType, hash, &area)) {
        return false;
    }

    return geohashDecodeAreaToLatLng(&area, lat, lng);
}

bool geohashDecodeToLatLngWGS84(const GeoHashBits hash, double *lat,
                                double *lng) {
    return geohashDecodeToLatLngType(GEO_WGS84_TYPE, hash, lat, lng);
}

bool geohashDecodeToLatLngMercator(const GeoHashBits hash, double *lat,
                                   double *lng) {
    return geohashDecodeToLatLngType(GEO_MERCATOR_TYPE, hash, lat, lng);
}

typedef enum moveWhich { MOVE_X = 3, MOVE_Y } moveWhich;

static void geohashMove(GeoHashBits *hash, int8_t direction, moveWhich which) {
    if (direction == 0) {
        return;
    }

    uint64_t x = hash->bits & 0xAAAAAAAAAAAAAAAALL;
    uint64_t y = hash->bits & 0x5555555555555555LL;

    uint64_t move = which == MOVE_X ? x : y;

    uint64_t zzMask =
        which == MOVE_X ? 0x5555555555555555LL : 0xAAAAAAAAAAAAAAAALL;
    uint64_t recombineMask =
        which == MOVE_X ? 0xAAAAAAAAAAAAAAAALL : 0x5555555555555555LL;

    uint64_t zz = zzMask >> (sizeof(zz) - hash->step * 2);
    if (direction > 0) {
        move = move + (zz + 1);
    } else {
        move = move | zz;
        move = move - (zz + 1);
    }

    move &= (recombineMask >> (sizeof(move) - hash->step * 2));

    hash->bits = which == MOVE_X ? (move | y) : (x | move);
}

void geohashNeighbors(const GeoHashBits *hash, GeoHashNeighbors *neighbors) {
    neighbors->east = *hash;
    neighbors->west = *hash;
    neighbors->north = *hash;
    neighbors->south = *hash;
    neighbors->southEast = *hash;
    neighbors->southWest = *hash;
    neighbors->northEast = *hash;
    neighbors->northWest = *hash;

    /* Solution grid:
     *   We have our center 'hash', we want every surrounding hash too.
     *   So, for [CENTER], get all corners.
     *          +--------------------------------------------+
     *          |              |               |             |
     *          |      NW      |       N       |      NE     |
     *          |    (-1, +1)  |    (0, +1)    |   (+1, +1)  |
     *          |              |               |             |
     *          |              |               |             |
     *          +--------------------------------------------+
     *          |              |               |             |
     *          |              |               |             |
     *          |      W       |               |      E      |
     *          |    (-1, 0)   |   Original    |   (+1, 0)   |
     *          |              |               |             |
     *          |              |               |             |
     *          +--------------------------------------------+
     *          |              |               |             |
     *          |       SW     |       S       |     SE      |
     *          |    (-1, -1)  |    (0, -1)    |   (+1, -1)  |
     *          |              |               |             |
     *          |              |               |             |
     *          +--------------------------------------------+
     */

    /* X goes Right */
    geohashMove(&neighbors->east, 1, MOVE_X);
    geohashMove(&neighbors->east, 0, MOVE_Y);

    /* X goes Left */
    geohashMove(&neighbors->west, -1, MOVE_X);
    geohashMove(&neighbors->west, 0, MOVE_Y);

    /* Y goes Down */
    geohashMove(&neighbors->south, 0, MOVE_X);
    geohashMove(&neighbors->south, -1, MOVE_Y);

    /* Y goes Up */
    geohashMove(&neighbors->north, 0, MOVE_X);
    geohashMove(&neighbors->north, 1, MOVE_Y);

    /* X goes Left, Y goes Up */
    geohashMove(&neighbors->northWest, -1, MOVE_X);
    geohashMove(&neighbors->northWest, 1, MOVE_Y);

    /* X goes Up, Y goes Up */
    geohashMove(&neighbors->northEast, 1, MOVE_X);
    geohashMove(&neighbors->northEast, 1, MOVE_Y);

    /* X goes Up, Y goes Down */
    geohashMove(&neighbors->southEast, 1, MOVE_X);
    geohashMove(&neighbors->southEast, -1, MOVE_Y);

    /* X goes Down, Y goes Down */
    geohashMove(&neighbors->southWest, -1, MOVE_X);
    geohashMove(&neighbors->southWest, -1, MOVE_Y);
}

bool geohashBuiltForBMI2(void) {
    return _geohashBuiltForBMI2;
}
