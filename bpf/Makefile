CXX			= gcc
LDFLAGS		= -pthread -lm
SOURCES		= main.c
INCLUDES	= -I.
OBJECTS		= $(SOURCES:.c=.o)
TARGET		= bpf

all: $(SOURCES) $(TARGET)

debug:	CXXFLAGS += -DDEBUG -ggdb -g3 -O0 -Wall -D_GNU_SOURCE -std=gnu99 -pthread 
debug:	$(SOURCES) $(TARGET)

release:	CXXFLAGS += -O3 -Wall -D_GNU_SOURCE -std=gnu99 -pthread 
release:	$(SOURCES) $(TARGET)

install:	$(SOURCES) $(TARGET)
	sudo cp $(TARGET) /usr/local/bin

$(TARGET): $(OBJECTS) 
	$(CXX) -o $(TARGET) $(OBJECTS) $(LDFLAGS)

.c.o:
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $<

.PHONY:
clean:
	rm -rf $(OBJECTS) $(TARGET)
