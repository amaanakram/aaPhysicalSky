ARNOLDLOC  := 
GXXLOC     := 
TARGETLOC  := 

INCLUDES := -I$(ARNOLDLOC)/include
LDFLAGS  := -L$(ARNOLDLOC)/bin -lai
CXXFLAGS := -g -fPIC -w -fvisibility=hidden -fopenmp
CXX      := $(GXXLOC)g++
SOURCES  := arnoldSkyShader.cpp
OBJECTS  := $(SOURCES:.cpp=.o)
TARGET   := $(TARGETLOC)aaPhysicalSky.so

$(TARGET) : $(OBJECTS)
    $(CXX) $(inputs) -shared $< -o $@ $(LDFLAGS)

%.o : %.cpp
    $(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@
    
all: $(TARGET)

clean:
    rm $(OBJECTS) $(TARGET)