ARNOLDLOC  := 
GXXLOC     := 
TARGETLOC  := 

INCLUDES := -I$(ARNOLDLOC)/include
LDFLAGS  := -L$(ARNOLDLOC)/bin -lai
CXXFLAGS := -g -O3 -fPIC -w -fvisibility=hidden -fopenmp
CXX      := $(GXXLOC)g++
SOURCES  := source/arnoldSkyShader.cpp
OBJECTS  := $(SOURCES:.cpp=.o)
TARGET   := $(TARGETLOC)aaPhysicalSky.so

$(TARGET) : $(OBJECTS)
    $(CXX) $(inputs) -shared $< -o $@ $(LDFLAGS)

%.o : %.cpp
    $(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@
    
all: $(TARGET)

clean:
    rm $(OBJECTS) $(TARGET)