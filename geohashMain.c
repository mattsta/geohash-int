#include <stdio.h>
#include <inttypes.h> /* PRIu64 */
#include "geohash.h"
#include "geohashHelper.h"

/* this is a sad mess.  we need globally relevant scans and bounding boxes */
int main() {
    GeoHashBits hash;
    hash.bits = 3;
    hash.step = 2;
    GeoHashNeighbors neighbors;
    geohashNeighbors(&hash, &neighbors);
    printf("%" PRIu64 "\n", hash.bits);
    printf("%" PRIu64 "\n", neighbors.east.bits);
    printf("%" PRIu64 "\n", neighbors.west.bits);
    printf("%" PRIu64 "\n", neighbors.south.bits);
    printf("%" PRIu64 "\n", neighbors.north.bits);
    printf("%" PRIu64 "\n", neighbors.northWest.bits);
    printf("%" PRIu64 "\n", neighbors.northEast.bits);
    printf("%" PRIu64 "\n", neighbors.southEast.bits);
    printf("%" PRIu64 "\n", neighbors.southWest.bits);

    const GeoHashRange *r = geohashRangeGetWGS84();
    printf("Encoding lat -32.1, lng 120.3\n");
    geohashEncode(r, -32.1, 120.3, 8, &hash);
    printf("Encoded using %d steps: %" PRIu64 "\n", hash.step, hash.bits);
    GeoHashArea area;
    geohashDecode(r, hash, &area);
    printf("Lat decoded: %.2f %.2f\n", area.range.lat.min, area.range.lat.max);
    printf("Lng decoded: %.2f %.2f\n", area.range.lng.min, area.range.lng.max);

    return 0;
}
