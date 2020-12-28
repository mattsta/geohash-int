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

#ifndef GEOHASH_INT_GEOHASH_H
#define GEOHASH_INT_GEOHASH_H

#include <stddef.h>
#include <stdbool.h>
#include <stdint.h>
__BEGIN_DECLS

#define HASHISZERO(r) (!(r).bits && !(r).step)
#define _RANGEISZERO(r) (!(r).max && !(r).min)
#define RANGEISZERO(r) (_RANGEISZERO((r)->lng) || _RANGEISZERO((r)->lat))

typedef enum GeoType {
    GEO_WGS84_TYPE = 16,
    GEO_MERCATOR_TYPE,
} GeoType;

#define GEOHASH_BIT_WIDTH 52
#define GEO_STEP_MAX 26

typedef enum GeoDirection {
    GEOHASH_NORTH = 7,
    GEOHASH_EAST,
    GEOHASH_WEST,
    GEOHASH_SOUTH,
    GEOHASH_SOUTH_WEST,
    GEOHASH_SOUTH_EAST,
    GEOHASH_NORT_WEST,
    GEOHASH_NORT_EAST
} GeoDirection;

typedef struct GeoHashBits {
    /* We currently use a maximum of 52 bits for the hash data.  That's because
     * we limit geocoding to the first 52 bits of a 'double'.
     * We can use the next 5 bits (l(26)/l(2)) to store the step size (max: 26).
     * Now we have 7 unused bits for other purposes. */

    /* If we store 'step' as an individual uint8_t, that causes GeoHashBits
     * to become 16 bytes with struct padding.  By using struct bitfields,
     * we keep the size of GeoHashBits to 8 bytes.
     * GeoHashNeighbors includes GeoHashBits 8 times and takes up 128 bytes
     * with excessive padding.  With the bitfield approach, GeoHashNeighbors
     * only uses 64 bytes.
     * GeoHashRadius uses 192 bytes with padded GeoHashBits and only
     * only 112 bytes using the struct bitfield version. */

    /* If we want to use full width, 64-bit unsigned geohashes in the future,
     * we may need to revisit this alignment.  Though, we could take six
     * unused bits below and extend 'bits' to 52+6 = 58 bits.
     *  52 bits is enough to geolocate down to a radius of 0.6 meters (~2 ft).
     *  58 bits is enough to geolocate down to a radius of 7.5 centimeters.
     *  64 bits is enough to geolocate down to a radius of 9.3 millimeters.
     * 128 bits is enough to geolocate down to a radius of 0.0021 nanometers.
     * (for comparison: DNA is 2 nanometers wide) */

    /* To calculate radius table for geohash bit length:
           radius = 20037726 # radius at 2 bits of precision; larger than earth
       for i in xrange(4, 66, 2):
                radius = radius / 2.0;
            print "{:2} bits (step {:2}) radius {:.4f} m".format(i, i/2, radius)
     */

    uint64_t bits : GEOHASH_BIT_WIDTH; /* currently 52 */
    uint64_t step : 5;
    uint64_t unused : 7;
} GeoHashBits;

typedef struct GeoHashRange {
    struct {
        double min;
        double max;
    } lng;
    struct {
        double min;
        double max;
    } lat;
} GeoHashRange;

typedef struct GeoHashArea {
    GeoHashBits hash;
    GeoHashRange range;
} GeoHashArea;

typedef struct GeoHashNeighbors {
    GeoHashBits north;
    GeoHashBits east;
    GeoHashBits west;
    GeoHashBits south;
    GeoHashBits northEast;
    GeoHashBits southEast;
    GeoHashBits northWest;
    GeoHashBits southWest;
} GeoHashNeighbors;

const GeoHashRange *geohashRangeGetWGS84(void);
const GeoHashRange *geohashRangeGetMercator(void);
void geohashRangeSetWGS84(GeoHashRange *range);
void geohashRangeSetMercator(GeoHashRange *range);
bool geohashGetCoordRange(GeoType coordType, GeoHashRange *range);
bool geohashEncode(const GeoHashRange *range, double lat, double lng,
                   uint8_t step, GeoHashBits *hash);
bool geohashEncodeType(uint8_t coordType, double lat, double lng, uint8_t step,
                       GeoHashBits *hash);
bool geohashEncodeMercator(double lat, double lng, uint8_t step,
                           GeoHashBits *hash);
bool geohashEncodeWGS84(double lat, double lng, uint8_t step,
                        GeoHashBits *hash);
bool geohashDecode(const GeoHashRange *range, const GeoHashBits hash,
                   GeoHashArea *area);
bool geohashDecodeType(uint8_t coordType, const GeoHashBits hash,
                       GeoHashArea *area);
bool geohashDecodeMercator(const GeoHashBits hash, GeoHashArea *area);
bool geohashDecodeWGS84(const GeoHashBits hash, GeoHashArea *area);
bool geohashDecodeAreaToLatLng(const GeoHashArea *area, double *lat,
                               double *lng);
bool geohashDecodeToLatLngType(uint8_t coordType, const GeoHashBits hash,
                               double *lat, double *lng);
bool geohashDecodeToLatLngWGS84(const GeoHashBits hash, double *lat,
                                double *lng);
bool geohashDecodeToLatLngMercator(const GeoHashBits hash, double *lat,
                                   double *lng);
void geohashNeighbors(const GeoHashBits *hash, GeoHashNeighbors *neighbors);
bool geohashBuiltForBMI2(void);

__END_DECLS
#endif /* GEOHASH_INT_GEOHASH_H */
