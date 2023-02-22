set(FIND_EIGEN3_PATHS ~/.local /usr /usr/local)
set(EIGEN3_DIR ~/.local)
find_path(EIGEN3_INCLUDE_DIRS
        NAMES eigen3
        PATH_SUFFIXES include
        PATHS ${FIND_EIGEN3_PATHS})

find_library(EIGEN3_LIBRARIES
        NAMES eigen3
        HINTS ${FIND_EIGEN3_PATHS})