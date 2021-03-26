#-----------------------------------------------------------------------------
PROGRAM = UnitedModelsScatteringGeant

SRC_DIR := src
OBJ_DIR := obj
HDR_DIR := include

SRCS := $(wildcard $(SRC_DIR)/*.cpp)
OBJS := $(SRCS:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)
HDRS := $(wildcard $(HDR_DIR)/*.h)

#-----------------------------------------------------------------------------
ObjSuf = o
SrcSuf = cpp
ExeSuf = exe
DllSuf = so
OutPutOpt = -o

ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)
ROOTGLIBS = $(shell root-config --glibs)

# debug option: -g
# no exceptions: -fno-exceptions

CXX = g++
CXXFLAGS = -O2 -Wall -fPIC -g -std=c++17
LD = g++
LDFLAGS = -O2 -Wall -g 
SOFLAGS = -shared
ARCHFLAGS =

#
CXXFLAGS += $(ROOTCFLAGS) -I$(HDR_DIR)
LIBS = $(ROOTLIBS)
GLIBS = $(ROOTGLIBS)

$(PROGRAM).so: $(OBJS) dict_$(PROGRAM).o
	@echo "Linking $(PROGRAM).so ..."
	$(LD) $(ARCHFLAGS) $(LDFLAGS) $(SOFLAGS) $(OBJS) dict_$(PROGRAM).o $(ROOTLIBS) -o$(PROGRAM).so
	@echo "linking done"

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	@echo "Making .o files ..."
	$(CXX) -c $(CXXFLAGS) $< -o $@ #$(SRCS) #*.cpp
	@echo "Done!"

dict_$(PROGRAM).o: dict_$(PROGRAM).cpp
	@echo "Compiling dictionary..."
	$(CXX) -c $(CXXFLAGS) $< -o $@

dict_$(PROGRAM).cpp: $(HDRS)
	@echo "Generating dictionary ..."
	rootcint -f dict_$(PROGRAM).cpp -c -p -I$(HDR_DIR) $(HDRS) Linkdef.h
	@echo "Done!"

# clean:
# 	@echo "Cleaning... "
# 	@rm -f $(OBJS) $(PROGRAM).so dict_*


#------------------------------------------------------------------------
