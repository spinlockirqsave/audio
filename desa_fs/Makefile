CXX			= gcc
LDFLAGS		= -pthread -lm
SOURCES		= desa3_fs.c desa2_fs.c
INCLUDES	= -I.
OBJECTS		= $(SOURCES:.c=.o)
TARGET		= desa3_fs

all: $(SOURCES) $(TARGET)

debug:	CXXFLAGS += -DDEBUG -g3 -O0 -Wall -D_GNU_SOURCE -std=gnu99 -pthread
debug:	$(SOURCES) $(TARGET)

release:	CXXFLAGS += -O3 -Wall -D_GNU_SOURCE -std=gnu99 -pthread
release:	$(SOURCES) $(TARGET)

$(TARGET): $(OBJECTS) 
	$(CXX) -o $(TARGET) $(OBJECTS) $(LDFLAGS)

.c.o:
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $<

.PHONY:
clean:
	rm -rf $(OBJECTS) $(TARGET)
