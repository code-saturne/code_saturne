cmake_minimum_required(VERSION 3.10)
PROJECT(CoProcessingTest C CXX)

include(CMakeDependentOption)
find_package(ParaView REQUIRED)

if (NOT TARGET ParaView::PythonCatalyst)
  message(STATUS
    "${CMAKE_PROJECT_NAME} requires ParaView to be built with Catalyst and "
    "Python support enabled. Please rebuild ParaView (or point to a "
    "different build of ParaView) with PARAVIEW_ENABLE_CATALYST and "
    "PARAVIEW_ENABLE_PYTHON set to TRUE")
  return ()
endif()

if (NOT PARAVIEW_USE_MPI)
  message(STATUS
    "${CMAKE_PROJECT_NAME} requires ParaView to be built with MPI support "
    "enabled. Please rebuild ParaView (or point to a different build of "
    "ParaView) with PARAVIEW_USE_MPI set to TRUE")
  return ()
endif ()

add_executable(CoProcessingTest CoProcessingTest.cxx)
target_link_libraries(CoProcessingTest
  PRIVATE
    ParaView::PythonCatalyst
    VTK::CommonDataModel
    VTK::ParallelMPI)

