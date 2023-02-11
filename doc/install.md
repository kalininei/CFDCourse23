@page install Сборка

Tested on gcc 7.5.0

External dependencies: boost, oneapi/tbb

Define path to the oneapi/tbb include directory in the root CMakeLists.txt with the THIRD_PARTY_INSTALL_PATH parameter.

Then
```
mkdir build
cd build
cmake ..
make -j
```

to run the example proceed with
```
cd bin
./tal_executor_test
```
