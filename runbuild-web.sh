set -x
set -e

source ../emsdk/emsdk_env.sh

emcmake cmake -B build -DUSE_EMSCRIPTEN=1 -DPLATFORM=Web
cmake --build build -j$(sysctl -n hw.ncpu) --verbose