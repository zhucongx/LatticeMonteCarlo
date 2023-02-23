# Try to find Eigen http://eigen.tuxfamily.org/
#
# Searches for Eigen library even if it is NOT installed (i.e. NO Eigen3Config.cmake)
#
# This module supports requiring a minimum version, e.g. you can do
#   find_package(Eigen 3.1.2)
# to require version 3.1.2 or newer of Eigen. By deafault looks for 3.0.0
#
# You can also simply use find_package(Eigen) or find_package(Eigen REQUIRED) or find_package(Eigen 3.2 REQUIRED)
#
# Once done the following variables will be defined
#
#  EIGEN_FOUND              - True if Eigen was found on your system
#  EIGEN_USE_FILE           - The file making Eigen usable by simply doing include(${EIGEN_USE_FILE})
#  EIGEN_INCLUDE_DIRS       - List of directories of Eigen and it's dependencies
#
#  EIGEN_DEFINITIONS        - Definitions needed to build with Eigen
#  EIGEN_INCLUDE_DIR        - Directory where signature_of_Eigen_matrix_library can be found
#  EIGEN_VERSION     		- Full Eigen Package version
#  EIGEN_VERSION_MAJOR      - The major version of Eigen
#  EIGEN_VERSION_MINOR      - The minor version of Eigen
#  EIGEN_VERSION_PATCH      - The patch version of Eigen

if(NOT Eigen_FIND_VERSION)
    if(NOT Eigen_FIND_VERSION_MAJOR)
        set(Eigen_FIND_VERSION_MAJOR 3)
    endif()
    if(NOT Eigen_FIND_VERSION_MINOR)
        set(Eigen_FIND_VERSION_MINOR 0)
    endif()
    if(NOT Eigen_FIND_VERSION_PATCH)
        set(Eigen_FIND_VERSION_PATCH 0)
    endif()
    set(Eigen_FIND_VERSION "${Eigen_FIND_VERSION_MAJOR}.${Eigen_FIND_VERSION_MINOR}.${Eigen_FIND_VERSION_PATCH}")
    message( STATUS "Looking for default minimum Eigen version ${Eigen_FIND_VERSION}" )
else()
    message( STATUS "Looking for Eigen Library with minimum version ${Eigen_FIND_VERSION}" )
endif()


# macro to check version
macro(Eigen_Check_Version)
    file(READ "${EIGEN_INCLUDE_DIR}/Eigen/src/Core/util/Macros.h" _Eigen_version_header)

    string(REGEX MATCH "define[ \t]+EIGEN_WORLD_VERSION[ \t]+([0-9]+)" _Eigen_world_version_match "${_Eigen_version_header}")
    set(EIGEN_VERSION_MAJOR "${CMAKE_MATCH_1}")
    string(REGEX MATCH "define[ \t]+EIGEN_MAJOR_VERSION[ \t]+([0-9]+)" _Eigen_major_version_match "${_Eigen_version_header}")
    set(EIGEN_VERSION_MINOR "${CMAKE_MATCH_1}")
    string(REGEX MATCH "define[ \t]+EIGEN_MINOR_VERSION[ \t]+([0-9]+)" _Eigen_minor_version_match "${_Eigen_version_header}")
    set(EIGEN_VERSION_PATCH "${CMAKE_MATCH_1}")
    set(EIGEN_VERSION ${EIGEN_VERSION_MAJOR}.${EIGEN_VERSION_MINOR}.${EIGEN_VERSION_PATCH})

    if(${EIGEN_VERSION} VERSION_LESS ${Eigen_FIND_VERSION})
        set(EIGEN_VERSION_OK FALSE)
    else()
        set(EIGEN_VERSION_OK TRUE)
    endif()

    message(STATUS "Eigen version ${EIGEN_VERSION} found in ${EIGEN_INCLUDE_DIR}")

    if(NOT EIGEN_VERSION_OK)
        message(WARNING "Eigen version is less than requred version ${Eigen_FIND_VERSION}")
    endif()
endmacro(Eigen_Check_Version)

#include the Standard package handler
include(FindPackageHandleStandardArgs)

if (EIGEN_INCLUDE_DIR)

    # in User Provided (or Cached) location
    message(STATUS "Looking for Eigen via User Provided (or Cached) location")
    Eigen_Check_Version()
    find_package_handle_standard_args(Eigen DEFAULT_MSG EIGEN_INCLUDE_DIR EIGEN_VERSION_OK)
    set(EIGEN_INCLUDE_DIRS ${EIGEN_INCLUDE_DIR} ${EIGEN_INCLUDE_DIR}/unsupported)

else (EIGEN_INCLUDE_DIR)

    ## Check for Installed Eigen Package (Eigen3Config.cmake)
    find_package(Eigen3 QUIET)

    if(EIGEN3_FOUND OR EIGEN_FOUND)
        # found via Installed Eigen3Config.cmake
        message(STATUS "Looking for Eigen via Installed Eigen3Config.cmake")
        if(EIGEN3_INCLUDE_DIR)
            set (EIGEN_INCLUDE_DIR ${EIGEN3_INCLUDE_DIR})
        endif(EIGEN3_INCLUDE_DIR)
    else(EIGEN3_FOUND OR EIGEN_FOUND)
        # Look for uninstalled Eigen Packages
        message(STATUS "Looking for Eigen at Standard Locations")

        find_path(EIGEN_INCLUDE_DIR NAMES signature_of_eigen3_matrix_library
                PATHS
                ${CMAKE_INSTALL_PREFIX}/include
                ${KDE4_INCLUDE_DIR}
                PATH_SUFFIXES eigen3 eigen
                )
    endif(EIGEN3_FOUND OR EIGEN_FOUND)

    if(EIGEN_INCLUDE_DIR)
        Eigen_Check_Version()
    endif(EIGEN_INCLUDE_DIR)

    find_package_handle_standard_args(Eigen DEFAULT_MSG EIGEN_INCLUDE_DIR EIGEN_VERSION_OK)

    set(EIGEN_INCLUDE_DIRS ${EIGEN_INCLUDE_DIR} )
    list(APPEND EIGEN_INCLUDE_DIRS ${EIGEN_INCLUDE_DIR}/unsupported )

    set(EIGEN_DEFINITIONS  "" )

endif(EIGEN_INCLUDE_DIR)
