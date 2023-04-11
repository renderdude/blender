#!/bin/bash

pushd ../build_xcode_cycles/bin/Debug/
ln -s ../../intern/cycles/kernel/osl/shaders shader
popd
cp ../lib/darwin/materialx/libraries/stdlib/genosl/include/*h ../build_xcode_cycles/bin/Debug/shader
cp ./intern/cycles/kernel/osl/shaders/stdcycles.h ../build_xcode_cycles/bin/Debug/shader
cp ../lib/darwin/osl/share/OSL/shaders/stdosl.h ../build_xcode_cycles/bin/Debug/shader
