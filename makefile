CXX := g++
RCINT := rootcint
RC := root-config

BDIR :=bin
LDIR :=lib

TOBJS := PlottingUtils.cxx PlottingIO.cxx
TOBJH := $(TOBJS:.cxx=.hxx)
TOBJH += PlottingTypes.hxx
TOBJSO := $(TOBJS:.cxx=.so)

TARGET := MakeSingleTransversePlots.exe
TARGETSRC := $(TARGET:.exe=.cxx)

ROOTCFLAGS := `$(RC) --cflags`
ROOTLDFLAGS := `$(RC) --libs --glibs`

CXXFLAGS := -fPIC $(ROOTCFLAGS) -g -std=c++11 -Wall
LDFLAGS := $(ROOTLDFLAGS) -Wl,-rpath=.

.PHONY: all clean

all: $(TARGET)
	@echo ""
	@echo "*********************************************************************"
	@echo "Success. Built plotting."
	@echo "*********************************************************************"
	mkdir -p $(BDIR) $(LDIR)
	mv $(TARGET) $(BDIR)/
	mv PlottingUtils.so PlottingIO.so $(LDIR)/

$(TARGET): $(TARGETSRC) $(TOBJSO) $(TOBJH)
	$(CXX) -o $@ $(CXXFLAGS) $(LDFLAGS) $(TARGETSRC) $(TOBJSO)

PlottingUtils.so: PlottingUtils.cxx PlottingUtils.hxx PlottingTypes.hxx
	$(CXX) $(CXXFLAGS) -shared PlottingUtils.cxx -o $@
PlottingIO.so: PlottingIO.cxx PlottingIO.hxx PlottingTypes.hxx PlottingUtils.hxx
	$(CXX) $(CXXFLAGS) -shared PlottingIO.cxx -o $@

clean:
	rm -rf $(TOBJSO)\
        $(TARGET)\
        $(BDIR) \
        $(LDIR)
