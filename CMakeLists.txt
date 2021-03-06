cmake_minimum_required(VERSION 3.10)

project(voro++ VERSION 0.4.6 LANGUAGES CXX)
set(SOVERSION "0")

option(BUILD_VORO_EXECUTABLE "build voro executable" OFF)
option(BUILD_SHARED_LIBS "Build shared libs" ON)
option(TEST "test example" OFF)
include(GNUInstallDirs)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CXX_FLAGS)
  #release comes with -O3 by default
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CXX_FLAGS)

set (CMAKE_CXX_STANDARD 11)

file(GLOB LIB_VORO_SRC src/*.cc)
file(GLOB NOT_VORO_SOURCES src/v_base_wl.cc src/cmd_line.cc src/voro++.cc)
list(REMOVE_ITEM LIB_VORO_SRC ${NOT_VORO_SOURCES})


SET(
    EXEC_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/src/cmd_line.cc
)


IF(BUILD_VORO_EXECUTABLE)
    add_executable(voro++ ${EXEC_SRC})
    target_link_libraries(voro++ PUBLIC voro++)
endif(BUILD_VORO_EXECUTABLE)


add_library(voro++ ${LIB_VORO_SRC})
set_target_properties(voro++ PROPERTIES 
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/src"
  SOVERSION ${SOVERSION})
install(TARGETS voro++ EXPORT VORO_Targets LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
#for voto++.hh
target_include_directories(voro++ PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src> $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)


IF(TEST)
  find_package(ZLIB REQUIRED)
  # file(GLOB EXAMPLES LIST_DIRECTORIES true test/*)
  file(GLOB test test/*.cpp)
  add_executable(qw46 ${test})
  include_directories(test/)
  target_link_libraries(qw46 PUBLIC voro++ ${ZLIB_LIBRARIES})
  enable_testing()
  add_test(NAME qw46 COMMAND qw46 ${CMAKE_CURRENT_SOURCE_DIR}/test/fcc.xyz)
endif(TEST)

file(GLOB_RECURSE VORO_HEADERS src/*.hh)
install(FILES ${VORO_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}) 
install(EXPORT VORO_Targets FILE VOROConfig.cmake NAMESPACE VORO:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/VORO)
include(CMakePackageConfigHelpers)
write_basic_package_version_file("VOROConfigVersion.cmake" VERSION ${PROJECT_VERSION} COMPATIBILITY ExactVersion)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/VOROConfigVersion.cmake" DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/VORO)