# Instructions

## Install NFLlib with AVX2 flag enabled

First, download and install [NFLlib](https://github.com/quarkslab/NFLlib). You need cmake, GMP and Mpfr, as well as a C++11 compiler to build NFLLib.
```sh
git clone https://github.com/quarkslab/NFLlib
mkdir deps
cd deps
cmake ../NFLlib -DCMAKE_BUILD_TYPE=Release -DNFL_OPTIMIZED=ON -DNTT_AVX2
make
make test
cd ../
```


## Compile & Run 
Run following commands to test verifiable encryption script:
```sh
mkdir build
cd build
cmake ..
make vericrypt
./vericrypt
```