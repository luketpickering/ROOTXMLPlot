CXX := g++
RCINT := rootcint
RC := root-config

BDIR :=bin
LDIR :=lib

TOBJS := PlottingIO.cxx PlottingSelections.cxx PlottingUtils.cxx  style_Xianguo.cxx PlottingAutoPlots.cxx
TOBJH := $(TOBJS:.cxx=.hxx)
TOBJH += PlottingTypes.hxx
TOBJ := $(TOBJS:.cxx=.o)

TARGET := MakeSingleTransversePlots.exe
TARGETSRC := $(TARGET:.exe=.cxx)

ROOTCFLAGS := `$(RC) --cflags`
ROOTLDFLAGS := `$(RC) --libs --glibs` -lXMLIO

UTILSLOC=../utils
UTILSDEPSO=libPureGenUtils.so

CXXFLAGS := -fPIC $(ROOTCFLAGS) -g -std=c++11 -Wall -I$(UTILSLOC)
LDFLAGS := $(ROOTLDFLAGS) -Wl,-rpath=.

.PHONY: all clean

all: $(TARGET)
	@echo ""
	@echo "*********************************************************************"
	@echo "Success. Built plotting."
	@echo "*********************************************************************"
	mkdir -p $(BDIR) $(LDIR)
	mv $(TARGET) $(BDIR)/
	mv $(UTILSDEPSO) $(LDIR)/

$(TARGET): $(TARGETSRC) $(TOBJ) $(TOBJH) $(UTILSDEPSO)
	$(CXX) -o $@ $(CXXFLAGS) $(LDFLAGS) $(TARGETSRC) $(TOBJ) $(UTILSDEPSO)

$(UTILSDEPSO):
	ln -s `readlink -f $(UTILSLOC)/lib/$@` $@

PlottingUtils.o: PlottingUtils.cxx PlottingUtils.hxx PlottingTypes.hxx
	$(CXX) $(CXXFLAGS) -c PlottingUtils.cxx -o $@

PlottingIO.o: PlottingIO.cxx PlottingIO.hxx PlottingTypes.hxx PlottingUtils.hxx
	$(CXX) $(CXXFLAGS) -c PlottingIO.cxx -o $@

PlottingSelections.o: PlottingSelections.cxx PlottingSelections.hxx PlottingUtils.hxx
	$(CXX) $(CXXFLAGS) -c PlottingSelections.cxx -o $@

PlottingAutoPlots.o: PlottingAutoPlots.cxx PlottingAutoPlots.hxx PlottingUtils.hxx
	$(CXX) $(CXXFLAGS) -c PlottingAutoPlots.cxx -o $@

style_Xianguo.o: style_Xianguo.cxx style_Xianguo.hxx
	$(CXX) $(CXXFLAGS) -c style_Xianguo.cxx -o $@

clean:
	rm -rf $(TOBJ)\
        $(TARGET)\
        $(BDIR) \
        $(LDIR)
