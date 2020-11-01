cmake_minimum_required(VERSION 3.10)
PROJECT(CoProcessingTest CXX C)

# include(CMakeDependentOption)
find_package(MPI COMPONENTS C CXX)
find_package(ParaView CONFIG COMPONENTS Catalyst PythonCatalyst)

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

# The entries after ParallelMPI are needed for the vtk histograms writer

add_executable(CoProcessingTest CoProcessingTest.cxx)
target_link_libraries(CoProcessingTest
  PRIVATE
    ParaView::PythonCatalyst
    VTK::CommonDataModel
    VTK::ParallelMPI
    VTK::IOImage
    VTK::ViewsContext2D
    VTK::RenderingContext2D
    VTK::ChartsCore
    VTK::png)

