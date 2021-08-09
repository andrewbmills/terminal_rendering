# cmake version to be used
cmake_minimum_required( VERSION 3.0 )

# project name
project( terminal_rendering )

# flags

# files

# include
include_directories( ./include ./src ./libs/Quartic-master ./libs)

# target
add_executable( main ./src/main.cpp ./src/shapes.cpp ./libs/Quartic-master/quartic.cpp )