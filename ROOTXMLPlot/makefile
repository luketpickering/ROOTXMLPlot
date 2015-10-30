CXX := g++
RCINT := rootcint
RC := root-config

BDIR :=bin
LDIR :=lib

TOBJS := PlottingIO.cxx PlottingSelections.cxx PlottingUtils.cxx  style_Xianguo.cxx PlottingAutoPlots.cxx
TOBJH := $(TOBJS:.cxx=.hxx)
TOBJH += PlottingTypes.hxx
TOBJ := $(TOBJS:.cxx=.o)

TARGET := ROOTXMLPlot.exe
TARGETSRC := $(TARGET:.exe=.cxx)

ROOTCFLAGS := `$(RC) --cflags`
ROOTLDFLAGS := `$(RC) --libs --glibs` -lXMLIO

UTILSLOC := ../utils
LIBUTILS := $(UTILSLOC)/lib/libPureGenUtils.a
LIBUTILSLD := -L$(UTILSLOC)/lib -lPureGenUtils

CXXFLAGS := -fPIC $(ROOTCFLAGS) -g -std=c++11 -Wall -I$(UTILSLOC)
LDFLAGS := $(ROOTLDFLAGS)  $(LIBUTILSLD) -Wl,-rpath=.

.PHONY: all clean

all: $(TARGET)
	@echo ""
	@echo "*********************************************************************"
	@echo "Success. Built plotting."
	@echo "*********************************************************************"
	mkdir -p $(BDIR)
	mv $(TARGET) $(BDIR)/

$(TARGET): $(TARGETSRC) $(TOBJ) $(TOBJH)
	$(CXX) -o $@ $(CXXFLAGS) $(TARGETSRC) $(TOBJ) $(LDFLAGS)

$(LIBUTILS):
	cd $(UTILSLOC); $(MAKE)

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
	rm -f $(TOBJ)\
        $(TARGET)\

distclean: clean
	rm -rf $(BDIR) \
	cd $(UTILSLOC); $(MAKE) clean
