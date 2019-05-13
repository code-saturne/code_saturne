cmake_minimum_required(VERSION 2.6)

PROJECT(CoProcessingTest)

FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})

ADD_EXECUTABLE(CoProcessingTest CoProcessingTest.cxx)
TARGET_LINK_LIBRARIES(CoProcessingTest vtkPVCatalyst vtkPVPythonCatalystPythonD)

