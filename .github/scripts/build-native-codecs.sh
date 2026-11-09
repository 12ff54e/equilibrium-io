#!/usr/bin/env bash
set -euo pipefail

# Usage: bash build-native-codecs.sh BUILD_DIRECTORY ON|OFF
# Match the Wasm dependency pins. Build NetCDF against this HDF5 when enabled.
mkdir -p "$1"
codec_root=$(cd "$1" && pwd)
codec_netcdf=$2
codec_prefix="$codec_root/install"
mkdir -p "$codec_root/downloads" "$codec_root/src"

download() {
    local name=$1 url=$2 checksum=$3
    local archive="$codec_root/downloads/$name.tar.gz"
    if [[ ! -f "$archive" ]]; then
        curl --fail --location --retry 3 "$url" --output "$archive"
    fi
    printf '%s  %s\n' "$checksum" "$archive" | sha256sum --check
    mkdir -p "$codec_root/src/$name"
    tar -xzf "$archive" --strip-components=1 -C "$codec_root/src/$name"
}

download hdf5 \
    https://codeload.github.com/HDFGroup/hdf5/tar.gz/refs/tags/hdf5_1.14.6 \
    09ee1c671a87401a5201c06106650f62badeea5a3b3941e9b1e2e1e08317357f
cmake -S "$codec_root/src/hdf5" -B "$codec_root/build-hdf5" -G Ninja \
    -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$codec_prefix" \
    -DCMAKE_INSTALL_LIBDIR=lib -DCMAKE_INSTALL_RPATH="$codec_prefix/lib" \
    -DBUILD_SHARED_LIBS=ON -DBUILD_STATIC_LIBS=OFF -DBUILD_TESTING=OFF \
    -DHDF5_BUILD_HL_LIB=ON -DHDF5_BUILD_CPP_LIB=OFF \
    -DHDF5_BUILD_FORTRAN=OFF -DHDF5_BUILD_JAVA=OFF \
    -DHDF5_BUILD_TOOLS=OFF -DHDF5_BUILD_UTILS=OFF \
    -DHDF5_BUILD_EXAMPLES=OFF -DHDF5_BUILD_DOC=OFF \
    -DHDF5_ENABLE_PARALLEL=OFF -DHDF5_ENABLE_SZIP_SUPPORT=OFF \
    -DHDF5_ENABLE_Z_LIB_SUPPORT=ON -DHDF5_ENABLE_PLUGIN_SUPPORT=OFF
cmake --build "$codec_root/build-hdf5" --parallel 2
cmake --install "$codec_root/build-hdf5"

if [[ "$codec_netcdf" == ON ]]; then
    download netcdf \
        https://codeload.github.com/Unidata/netcdf-c/tar.gz/refs/tags/v4.9.3 \
        990f46d49525d6ab5dc4249f8684c6deeaf54de6fec63a187e9fb382cc0ffdff
    cmake -S "$codec_root/src/netcdf" -B "$codec_root/build-netcdf" -G Ninja \
        -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$codec_prefix" \
        -DCMAKE_INSTALL_LIBDIR=lib -DCMAKE_INSTALL_RPATH="$codec_prefix/lib" \
        -DCMAKE_PREFIX_PATH="$codec_prefix" -DHDF5_ROOT="$codec_prefix" \
        -DBUILD_SHARED_LIBS=ON -DNETCDF_FIND_SHARED_LIBS=ON \
        -DNETCDF_ENABLE_HDF5=ON -DNETCDF_ENABLE_DAP=OFF \
        -DNETCDF_ENABLE_NCZARR=OFF -DNETCDF_ENABLE_REMOTE_FUNCTIONALITY=OFF \
        -DNETCDF_ENABLE_PLUGINS=OFF -DNETCDF_ENABLE_LIBXML2=OFF \
        -DNETCDF_BUILD_UTILITIES=OFF -DNETCDF_ENABLE_EXAMPLES=OFF \
        -DNETCDF_ENABLE_TESTS=OFF -DNETCDF_ENABLE_DOXYGEN=OFF
    cmake --build "$codec_root/build-netcdf" --parallel 2
    cmake --install "$codec_root/build-netcdf"
fi
