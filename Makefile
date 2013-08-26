CC=g++
CFLAGS=-O3 -Wno-write-strings
SOURCE=source/photodynam.cpp source/n_body.cpp source/n_body_state.cpp source/n_body_lc.cpp source/elliptic.c source/icirc.c source/scpolyint.c source/mttr.c
EXECUTABLE=photodynam

all: $(EXECUTABLE)

$(EXECUTABLE): $(SOURCE)
	$(CC) -o $(EXECUTABLE) $(SOURCE)

