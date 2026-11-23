cmake_minimum_required(VERSION 3.10)

project(dist)

option(CXX "enable C++ compilation" ON)
enable_language(CXX)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

file(GLOB MISC misc/*.hpp misc/*.cpp)
file(GLOB GRIMM GRIMM/*.h GRIMM/*.c)
file(GLOB DIST distance_algorithms/*.hpp distance_algorithms/*.cpp distance_algorithms/perm/*.c distance_algorithms/perm/*.h distance_algorithms/perm/*/*.c distance_algorithms/perm/*/*.h)
file(GLOB EXTER external/*.hpp external/*.cpp)
file(GLOB CYCLE cycle/*.hpp cycle/*.cpp)
file(GLOB HEUR heur/*.hpp heur/*.cpp)

set(CMAKE_CXX_STANDARD 14)
set(CXX_FILESYSTEM_LIBRARIES "stdc++fs")

#################################################################################################
# dist
#################################################################################################

find_package(OpenMP REQUIRED)
add_executable(${CMAKE_PROJECT_NAME} main_dist.cpp ${MISC} ${DIST} ${EXTER} ${CYCLE} ${GRIMM})
target_compile_options(dist PRIVATE -Wall PRIVATE "${OpenMP_CXX_FLAGS}")
target_link_libraries(dist PUBLIC ${CXX_FILESYSTEM_LIBRARIES} PRIVATE "${OpenMP_CXX_FLAGS}")

#################################################################################################
# dist_test
#################################################################################################

enable_testing()
find_package(OpenMP REQUIRED)
file(GLOB TEST quickcheck/quickcheck/*.hh)
add_executable(dist_test test.cpp ${MISC} ${DIST} ${TEST} ${EXTER} ${CYCLE} ${HEUR} ${GRIMM})
target_compile_options(dist_test PRIVATE -Wall PRIVATE "${OpenMP_CXX_FLAGS}")
target_link_libraries(dist_test PUBLIC ${CXX_FILESYSTEM_LIBRARIES} PRIVATE "${OpenMP_CXX_FLAGS}")

# One test per check, so that a failed assert only aborts its own check.
# Keep in sync with the checks list in test.cpp.
set(TEST_CHECKS
    PBounds3T PBounds4R PBounds45RT
    CheckGraphDegreesMap CheckGraphDecCorrespondenceMap
    CheckGraphDegreesBFS CheckGraphDecCorrespondenceBFS
    CheckGraphDegrees1C2C3C CheckGraphDecCorrespondence1C2C3C
    CheckGraphDecCorrespondenceGAMap CheckGraphDecCorrespondenceGABFS
    CheckGraphDecCorrespondenceGA1C2C3C
    GABoundsMap GABoundsBFS GABounds1C2C3C)
foreach(check ${TEST_CHECKS})
    add_test(NAME ${check} COMMAND dist_test ${check})
endforeach()

#################################################################################################
# decomposition
#################################################################################################

find_package(OpenMP REQUIRED)
add_executable(dec main_dec.cpp ${MISC} ${CYCLE} ${HEUR} ${GRIMM})
# target_compile_options(dec PRIVATE -Wall PRIVATE "")
# target_link_libraries(dec PUBLIC ${CXX_FILESYSTEM_LIBRARIES} PRIVATE "")
target_compile_options(dec PRIVATE -Wall PRIVATE "${OpenMP_CXX_FLAGS}")
target_link_libraries(dec PUBLIC ${CXX_FILESYSTEM_LIBRARIES} PRIVATE "${OpenMP_CXX_FLAGS}")

#################################################################################################
# To permutation
#################################################################################################

find_package(OpenMP REQUIRED)
add_executable(toperm main_toperm.cpp ${MISC} ${CYCLE})
target_compile_options(toperm PRIVATE -Wall PRIVATE "${OpenMP_CXX_FLAGS}")
target_link_libraries(toperm PUBLIC ${CXX_FILESYSTEM_LIBRARIES} PRIVATE "${OpenMP_CXX_FLAGS}")

#################################################################################################
# partition
#################################################################################################

find_package(OpenMP REQUIRED)
add_executable(part main_part.cpp ${MISC} ${CYCLE} ${HEUR} ${GRIMM})
target_compile_options(part PRIVATE -Wall PRIVATE "${OpenMP_CXX_FLAGS}")
target_link_libraries(part PUBLIC ${CXX_FILESYSTEM_LIBRARIES} PRIVATE "${OpenMP_CXX_FLAGS}")

#################################################################################################

#################################################################################################


if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
    include(FeatureSummary)
    feature_summary(WHAT ALL)
endif()
