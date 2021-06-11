
# Location of our plugins
#
PLUGINDIR	= plugins

# Location of LibXtract
#
LIBXTRACTDIR	= LibXtract

# Compile flags
#
CFLAGS		:= $(CFLAGS) -DXTRACT_FFT=1 -DUSE_OOURA=1 -DNDEBUG -O3 -ffast-math -Wall -fPIC -I. -I$(LIBXTRACTDIR)
CXXFLAGS	:= $(CXXFLAGS) $(CFLAGS) 

# Libraries required for the plugins.
#
PLUGIN_LIBS	= -Wl,-Bstatic -lvamp-sdk -lfftw3f -Wl,-Bdynamic

# Flags required to tell the compiler to make a dynamically loadable object
#
PLUGIN_LDFLAGS	= $(LDFLAGS) -shared -Wl,-Bsymbolic -Wl,--version-script=vamp-plugin.map

# File extension for a dynamically loadable object
#
PLUGIN_EXT	= .so

## For OS/X with g++:
#PLUGIN_LDFLAGS	= -dynamiclib -exported_symbols_list=vamp-plugin.list
#PLUGIN_EXT	= .dylib


### End of user-serviceable parts

PLUGIN_OBJECTS	= libmain.o $(patsubst %.cpp,%.o,$(wildcard $(PLUGINDIR)/*.cpp))
XTRACT_OBJECTS	= $(patsubst %.c,%.o,$(wildcard $(LIBXTRACTDIR)/src/*.c $(LIBXTRACTDIR)/src/*/*.c))
PLUGIN_HEADERS	= $(patsubst %.cpp,%.h,$(wildcard $(PLUGINDIR)/*.cpp))
PLUGIN_TARGET	= vamp-libxtract$(PLUGIN_EXT)

all:		$(PLUGIN_TARGET)

$(PLUGIN_TARGET):	$(PLUGIN_OBJECTS) $(XTRACT_OBJECTS) $(PLUGIN_HEADERS)
		$(CXX) $(LDFLAGS) $(PLUGIN_LDFLAGS) -o $@ $(PLUGIN_OBJECTS) $(XTRACT_OBJECTS) $(PLUGIN_LIBS)

clean:		
		rm -f $(PLUGIN_OBJECTS) $(XTRACT_OBJECTS)

distclean:	clean
		rm -f $(PLUGIN_TARGET) *~ */*~


