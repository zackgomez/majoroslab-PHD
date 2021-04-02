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
		$(OBJ)/ReadVariants.o \
		$(OBJ)/ConnectedComponent.o \
		$(OBJ)/Allele.o \
		$(OBJ)/Phase.o \
		$(OBJ)/SamReader.o \
		$(OBJ)/SamRecord.o \
		$(OBJ)/Variant.o \
		$(OBJ)/VariantGraph.o \
		$(OBJ)/VariantInRead.o \
		$(OBJ)/phd.o
	$(CC) $(LDFLAGS) -o phd \
		$(OBJ)/ReadVariants.o \
		$(OBJ)/ConnectedComponent.o \
		$(OBJ)/Allele.o \
		$(OBJ)/Phase.o \
		$(OBJ)/SamReader.o \
		$(OBJ)/SamRecord.o \
		$(OBJ)/Variant.o \
		$(OBJ)/VariantGraph.o \
		$(OBJ)/VariantInRead.o \
		$(OBJ)/phd.o \
		$(LIBS)
#---------------------------------------------------------
$(OBJ)/subset-vcf-by-samples.o:\
		subset-vcf-by-samples.C
	$(CC) $(CFLAGS) -o $(OBJ)/subset-vcf-by-samples.o -c \
		subset-vcf-by-samples.C
#---------------------------------------------------------
subset-vcf-by-samples: \
		$(OBJ)/subset-vcf-by-samples.o
	$(CC) $(LDFLAGS) -o subset-vcf-by-samples \
		$(OBJ)/subset-vcf-by-samples.o \
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

#--------------------------------------------------------
$(OBJ)/Variant.o:\
		Variant.C\
		Variant.H
	$(CC) $(CFLAGS) -o $(OBJ)/Variant.o -c \
		Variant.C
#---------------------------------------------------------

#--------------------------------------------------------
$(OBJ)/VariantGraph.o:\
		VariantGraph.C\
		VariantGraph.H
	$(CC) $(CFLAGS) -o $(OBJ)/VariantGraph.o -c \
		VariantGraph.C
#---------------------------------------------------------

#--------------------------------------------------------
$(OBJ)/VariantInRead.o:\
		VariantInRead.C\
		VariantInRead.H
	$(CC) $(CFLAGS) -o $(OBJ)/VariantInRead.o -c \
		VariantInRead.C
#---------------------------------------------------------

#--------------------------------------------------------
$(OBJ)/ReadVariants.o:\
		ReadVariants.C\
		ReadVariants.H
	$(CC) $(CFLAGS) -o $(OBJ)/ReadVariants.o -c \
		ReadVariants.C
#---------------------------------------------------------

#--------------------------------------------------------
$(OBJ)/Phase.o:\
		Phase.C\
		Phase.H
	$(CC) $(CFLAGS) -o $(OBJ)/Phase.o -c \
		Phase.C
#--------------------------------------------------------
$(OBJ)/Allele.o:\
		Allele.C\
		Allele.H
	$(CC) $(CFLAGS) -o $(OBJ)/Allele.o -c \
		Allele.C
#---------------------------------------------------------

#--------------------------------------------------------
$(OBJ)/ConnectedComponent.o:\
		ConnectedComponent.C\
		ConnectedComponent.H
	$(CC) $(CFLAGS) -o $(OBJ)/ConnectedComponent.o -c \
		ConnectedComponent.C
#---------------------------------------------------------
#---------------------------------------------------------
$(OBJ)/count-isoform-variants.o:\
		count-isoform-variants.C
	$(CC) $(CFLAGS) -o $(OBJ)/count-isoform-variants.o -c \
		count-isoform-variants.C
#---------------------------------------------------------
count-isoform-variants: \
		$(OBJ)/Variant.o \
		$(OBJ)/count-isoform-variants.o
	$(CC) $(LDFLAGS) -o count-isoform-variants \
		$(OBJ)/Variant.o \
		$(OBJ)/count-isoform-variants.o \
		$(LIBS)

#--------------------------------------------------------
$(OBJ)/count-combinations-in-fragments.o:\
		count-combinations-in-fragments.C
	$(CC) $(CFLAGS) -o $(OBJ)/count-combinations-in-fragments.o -c \
		count-combinations-in-fragments.C
#---------------------------------------------------------
count-combinations-in-fragments: \
		$(OBJ)/count-combinations-in-fragments.o
	$(CC) $(LDFLAGS) -o count-combinations-in-fragments \
		$(OBJ)/count-combinations-in-fragments.o \
		$(LIBS)
#---------------------------------------------
