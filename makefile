CC		= g++
DEBUG		= -g
OPTIMIZE	= -O
CFLAGS		= $(OPTIMIZE) -fpermissive -w
LDFLAGS		= $(OPTIMIZE)
BOOM		= BOOM
OBJ		= obj
LIBS		= -LBOOM -lBOOM -lgsl -lm -lgslcblas

.PHONY : clean
clean:
	@rm -f obj/*.o

#---------------------------------------------------------
$(OBJ)/phd.o:\
		phd.C
	$(CC) $(CFLAGS) -o $(OBJ)/phd.o -c \
		phd.C
#---------------------------------------------------------
phd: \
		$(OBJ)/SamReader.o \
		$(OBJ)/SamRecord.o \
		$(OBJ)/phd.o
	$(CC) $(LDFLAGS) -o phd \
		$(OBJ)/SamReader.o \
		$(OBJ)/SamRecord.o \
		$(OBJ)/phd.o \
		$(LIBS)
#---------------------------------------------------------
$(OBJ)/sort-gff.o:\
		sort-gff.C
	$(CC) $(CFLAGS) -o $(OBJ)/sort-gff.o -c \
		sort-gff.C
#---------------------------------------------------------
sort-gff: \
		$(OBJ)/sort-gff.o
	$(CC) $(LDFLAGS) -o sort-gff \
		$(OBJ)/sort-gff.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/SamReader.o:\
		SamReader.C\
		SamReader.H
	$(CC) $(CFLAGS) -o $(OBJ)/SamReader.o -c \
		SamReader.C
#--------------------------------------------------------
$(OBJ)/SamRecord.o:\
		SamRecord.C\
		SamRecord.H
	$(CC) $(CFLAGS) -o $(OBJ)/SamRecord.o -c \
		SamRecord.C
#---------------------------------------------------------
