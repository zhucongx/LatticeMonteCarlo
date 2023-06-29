cmake_minimum_required(VERSION 3.14)
#################################################### SET  COMPILERS ####################################################
#set(CMAKE_C_COMPILER "gcc-13") # {gcc, gcc-13, mpicc, icc, mpiicc, clang}
#set(CMAKE_CXX_COMPILER "g++-13") # {g++, g++-13, mpicxx, icpc, mpiicpc, clang++}
#set(CMAKE_Fortran_COMPILER "gfortran-13") # {gfortran, gfortran-13, mpifort, ifort, mpiifort, flang}
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_C_STANDARD 17)

project(lmc VERSION 0.1.0)
set(CMAKE_VERBOSE_MAKEFILE ON)
###################################################### SET  PATHS ######################################################
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
set(EXECUTABLE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
###################################################### SET  FLAGS ######################################################
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}"
        CACHE INTERNAL "")

#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} \
set(CMAKE_CXX_FLAGS_DEBUG "\
-Wall -Wextra -Wshadow -Wnon-virtual-dtor -Wpointer-arith -Wold-style-cast \
-Wpedantic -Wcast-qual -Wunused -Woverloaded-virtual -Winline -Wundef \
-Wunreachable-code -Wcast-align -Weffc++ -Wmisleading-indentation \
-Wredundant-decls -pedantic-errors -Wnull-dereference -Wdouble-promotion \
-Wconversion -Wfatal-errors -Wextra-semi -O0 -g3"
        CACHE INTERNAL "")

#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} \
set(CMAKE_CXX_FLAGS_RELEASE "\
-O3 -finline-limit=1000000 -Wall"
        CACHE INTERNAL "")

message(STATUS "Detecting ${CMAKE_CXX_COMPILER_ID} compiler")
message(STATUS "Current build type is ${CMAKE_BUILD_TYPE}")
IF (CMAKE_BUILD_TYPE MATCHES Debug)
    SET(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS_DEBUG})
ELSEIF (CMAKE_BUILD_TYPE MATCHES Release)
    SET(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS_RELEASE})
ENDIF ()

#if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
#
#elseif (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
#
#elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
#
#elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
#
#else ()
###################################################### BUILD SRCS ######################################################
file(GLOB PROJECT_SRC ${PROJECT_SOURCE_DIR}/${PROJECT_NAME}/*/src/*.cpp)
file(GLOB PROJECT_SOURCE ${PROJECT_SOURCE_DIR}/${PROJECT_NAME}/*/src)
foreach (dir ${PROJECT_SOURCE})
    message(STATUS "Found source directory: ${dir}")
endforeach ()
file(GLOB PROJECT_INCLUDE ${PROJECT_SOURCE_DIR}/${PROJECT_NAME}/*/include)
foreach (dir ${PROJECT_INCLUDE})
    message(STATUS "Found include directory: ${dir}")
endforeach ()

add_executable(${PROJECT_NAME} ${PROJECT_NAME}/main.cpp ${PROJECT_SRC})
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_INCLUDE})
set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "${PROJECT_NAME}.exe")
###################################################### BUILD LIBS ######################################################
include(FetchContent)
####################################################### PKGCONFIG ######################################################
#find_package(PkgConfig)
######################################################## BOOST #########################################################
#set(Boost_USE_DEBUG_LIBS OFF) # ignore debug libs and
#set(Boost_USE_RELEASE_LIBS ON) # only find release libs
find_package(Boost COMPONENTS filesystem system iostreams REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE Boost::boost Boost::filesystem Boost::system Boost::iostreams)
#target_include_directories(${PROJECT_NAME} SYSTEM PRIVATE ${Boost_INCLUDE_DIR})
######################################################### MPI ##########################################################
find_package(MPI REQUIRED)
#target_compile_options(${PROJECT_NAME} PRIVATE ${MPI_CXX_COMPILE_FLAGS})
target_link_libraries(${PROJECT_NAME} PRIVATE MPI::MPI_CXX)
#target_include_directories(${PROJECT_NAME} PRIVATE ${MPI_CXX_INCLUDE_DIRS})
######################################################## OPENMP ########################################################
find_package(OpenMP REQUIRED)
#target_compile_options(${PROJECT_NAME} PRIVATE ${OpenMP_CXX_FLAGS})
target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)
#target_include_directories(${PROJECT_NAME} PRIVATE ${OpenMP_CXX_INCLUDE_DIRS})
######################################################## EIGEN #########################################################
find_package(MKL)
if (MKL_FOUND AND NOT (CMAKE_CXX_COMPILER_ID MATCHES "Clang"))
    message(STATUS "MKL found. Using MKL BLAS/LAPACK.")
    set(BLA_VENDOR Intel10_64lp)
    add_definitions(-DEIGEN_USE_MKL_ALL)
    target_compile_options(${PROJECT_NAME} PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
    target_link_libraries(${PROJECT_NAME} PUBLIC $<LINK_ONLY:MKL::MKL>)
    target_include_directories(${PROJECT_NAME} PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
else ()
    message(STATUS "MKL not found.")
    if (DEFINED ENV{MKLROOT})
        message(STATUS "But MKLROOT found, try MKL again: $ENV{MKLROOT}")
        set(BLA_VENDOR Intel10_64lp)
        add_definitions(-DEIGEN_USE_MKL_ALL)
        if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
            add_definitions(-mkl)
        endif ()
    else ()
        message(WARNING "MKLROOT not found. Using generic BLAS/LAPACK.")
        set(BLA_VENDOR Generic)
        add_definitions(-DEIGEN_USE_LAPACK -DEIGEN_USE_BLAS)
    endif ()
    find_package(BLAS)
    find_package(LAPACK)
    if (BLAS_FOUND AND LAPACK_FOUND)
        target_compile_options(${PROJECT_NAME} PUBLIC ${BLAS_LINKER_FLAGS})
        target_link_libraries(${PROJECT_NAME} PUBLIC ${BLAS_LIBRARIES})
        target_include_directories(${PROJECT_NAME} PUBLIC ${BLAS_INCLUDE_DIRS})
        target_compile_options(${PROJECT_NAME} PUBLIC ${LAPACK_LINKER_FLAGS})
        target_link_libraries(${PROJECT_NAME} PUBLIC ${LAPACK_LIBRARIES})
        target_include_directories(${PROJECT_NAME} PUBLIC ${LAPACK_INCLUDE_DIRS})
    else ()
        message(WARNING "BLAS/LAPACK not found. Using Eigen3 only.")
    endif ()
endif ()

find_package(Eigen3 3.3)
if (NOT EIGEN3_FOUND)
    message(WARNING "Eigen3 not found. Downloading Eigen3 from gitlab")
    FetchContent_Declare(
            Eigen3
            GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
            GIT_TAG 3.4.0
            GIT_PROGRESS TRUE
            GIT_SHALLOW TRUE
    )
    FetchContent_MakeAvailable(Eigen3)
    set(EIGEN3_INCLUDE_DIRS ${eigen3_SOURCE_DIR})
endif ()
target_include_directories(${PROJECT_NAME} SYSTEM PRIVATE ${EIGEN3_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} PRIVATE Eigen3::Eigen)
###################################################### ARMADILLO #######################################################
#find_package(ARMADILLO REQUIRED)
#target_include_directories(${PROJECT_NAME} PUBLIC
#        ${ARMADILLO_INCLUDE_DIRS})
#target_link_libraries(${PROJECT_NAME} PUBLIC
#        ${ARMADILLO_LIBRARIES})
###################################################### TENSORFLOW ######################################################
#find_package(TensorFlow REQUIRED)
#target_include_directories(${PROJECT_NAME} PUBLIC
#        ${TensorFlow_INCLUDE_DIRS})
#target_link_libraries(${PROJECT_NAME} PUBLIC
#        ${TensorFlow_LIBRARIES})
######################################################### JSON #########################################################
# Typically you don't care so much for a third party library's tests to be run from your own project's code.
# set(JSON_BuildTests OFF CACHE INTERNAL "")
# If you only include this third party in PRIVATE source files, you do not
# need to install it when your main project gets installed.
# set(JSON_Install OFF CACHE INTERNAL "")
# add_subdirectory(external/json)
FetchContent_Declare(
        json
        GIT_REPOSITORY https://github.com/ArthurSonzogni/nlohmann_json_cmake_fetchcontent
        GIT_TAG v3.11.2
        GIT_PROGRESS TRUE
        GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(json)
target_link_libraries(${PROJECT_NAME} PRIVATE nlohmann_json::nlohmann_json)
######################################################## SPGLIB ########################################################
#FetchContent_Declare(
#        Spglib
#        GIT_REPOSITORY https://github.com/spglib/spglib.git
#        GIT_TAG v2.1.0-rc2
#        GIT_PROGRESS TRUE
#        GIT_SHALLOW TRUE
#)
#set(SPGLIB_USE_OMP ON)
#FetchContent_MakeAvailable(Spglib)
##FetchContent_GetProperties(Spglib)
##if (NOT Spglib_POPULATED)
##    FetchContent_Populate(Spglib)
##    add_subdirectory(${spglib_SOURCE_DIR} ${spglib_BINARY_DIR} EXCLUDE_FROM_ALL)
##endif ()
#target_link_libraries(${PROJECT_NAME} PRIVATE Spglib::symspg)
###################################################### GOOGLETEST ######################################################
#set(TEST "on")
#IF(TEST)
#    message(STATUS "TEST mode on")
#    include(FetchContent)
#    FetchContent_Declare(
#            googletest
#            GIT_REPOSITORY https://github.com/google/googletest/
#            GIT_TAG release-1.11.0)
#    # For Windows: Prevent overriding the parent project's compiler/linker settings
#    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
#    FetchContent_MakeAvailable(googletest)
#
#    enable_testing()
#    include(GoogleTest)
##    macro(package_add_test TESTNAME)
##        add_executable(${TESTNAME} ${ARGN})
##        target_link_libraries(${TESTNAME} gtest gmock gtest_main)
##        gtest_discover_tests(${TESTNAME}
##                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test
##                PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/test"
##                )
##        set_target_properties(${TESTNAME} PROPERTIES FOLDER tests)
##    endmacro()
#
##    macro(package_add_test_with_libraries TESTNAME FILES LIBRARIES)
##        add_executable(${TESTNAME} ${FILES})
##        target_link_libraries(${TESTNAME} gtest gmock gtest_main ${LIBRARIES})
##        gtest_discover_tests(${TESTNAME}
##                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test
##                PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/test"
##                )
##        set_target_properties(${TESTNAME} PROPERTIES FOLDER tests)
##    endmacro()
#    set(TEST_NAME my_test)
#    add_executable(${TEST_NAME} gtest.cpp ${PROJECT_SRC})
#    set_target_properties(${TEST_NAME} PROPERTIES
#            OUTPUT_NAME "${TEST_NAME}.exe")
#    target_include_directories(${TEST_NAME} PUBLIC
#            ${PROJECT_INCLUDE}
#            ${Boost_INCLUDE_DIR})
#
#    target_link_libraries(${TEST_NAME} gtest gmock gtest_main )
#    gtest_discover_tests(${TEST_NAME})
#ELSE()
#    message(STATUS "TEST mode off")
#ENDIF()

##################################################### SHOW SUMMARY #####################################################
MESSAGE(STATUS "====================================================================================================")
MESSAGE(STATUS "Project Name:     ${PROJECT_NAME}")
MESSAGE(STATUS "C++ Compiler:     ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
MESSAGE(STATUS "C++ Standard:     C++ ${CMAKE_CXX_STANDARD}")
MESSAGE(STATUS "Build type:       ${CMAKE_BUILD_TYPE}")
MESSAGE(STATUS "Build type flags: ${CMAKE_CXX_FLAGS}")
MESSAGE(STATUS "Output Path:      ${EXECUTABLE_OUTPUT_DIRECTORY}")
MESSAGE(STATUS "====================================================================================================")
