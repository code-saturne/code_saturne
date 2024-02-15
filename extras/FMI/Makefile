PROJECT_NAME   = FMU_CS
BUILD_DIR      = "binaries/linux64"
CXX            = g++
CXXFLAGS       = -m64 -pipe -std=c++11 -g -fPIC -Wall -W
LINK           = g++
LFLAGS         = -m64 -shared
INCPATH        = -I gen -I fmi -I.
MKDIR          = mkdir -p
ZIP            = zip
CHK_DIR_EXISTS = test -d

TARGET         = $(PROJECT_NAME).so

all: fmu

build:
	$(CHK_DIR_EXISTS) $(BUILD_DIR) || $(MKDIR) $(BUILD_DIR)
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $(BUILD_DIR)/$(PROJECT_NAME).o src/code_saturne.cpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $(BUILD_DIR)/callbacks.o src/callbacks.cpp
	$(LINK) $(LFLAGS) -o $(BUILD_DIR)/$(TARGET) $(BUILD_DIR)/$(PROJECT_NAME).o $(BUILD_DIR)/callbacks.o
	rm -f $(BUILD_DIR)/*.o
	if ! test -h $(BUILD_DIR)/lib$(TARGET) ; then \
	  \ln -s -r $(BUILD_DIR)/$(TARGET) $(BUILD_DIR)/lib$(TARGET) ; \
	fi

fmu: build
	$(ZIP) -r FMU_CS.fmu modelDescription.xml binaries

clean:
	rm -f $(BUILD_DIR)/*.o

distclean:

	rm -rf binaries

