set(FIND_EIGEN3_PATHS ~/.local /usr /usr/local /usr/share ${C_PATH})
find_path(EIGEN3_INCLUDE_DIRS
        NAMES eigen3
#        PATH_SUFFIXES include
        PATHS ${FIND_EIGEN3_PATHS})

#find_library(EIGEN3_LIBRARIES
#        NAMES eigen3
#        HINTS ${FIND_EIGEN3_PATHS})