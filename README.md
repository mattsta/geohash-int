Feature Hacking Branch
======================

Check out the fancy SSE4/BMI2 single instruction bit interleaving for geohash creation: https://github.com/mattsta/geohash-int/blob/599c7532839a136525cc96bf8c653cad9cb117d3/geohash.c#L73-L91

Also includes other refactorings against the original `geohash-int` release from 2014, but code has diverged somewhat since then.

## Building:

```bash
mkdir build
cd build
cmake ..
make -j12
```

geohash-int
======

A fast C99 geohash library which only provide int64 hash result. 
[GeoHash](http://en.wikipedia.org/wiki/Geohash).

## Embedding geohash-int

Just copy all files  into your project. 




