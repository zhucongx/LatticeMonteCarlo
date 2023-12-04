set(FIND_TensorFlow_PATHS ~/TensorFlow /usr/ /usr/local)

find_path(TensorFlow_INCLUDE_DIRS
        NAMES tensorflow
        PATH_SUFFIXES include
        PATHS ${FIND_TensorFlow_PATHS})
find_library(TensorFlow_LIBRARIES
        NAMES tensorflow
        PATH_SUFFIXES lib
        HINTS ${FIND_TensorFlow_PATHS})