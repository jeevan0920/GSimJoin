COMPILER = g++
CCFLAGS = -O3
CCINCLUDES = -I.
ALLOBJECTS = mai.o io.o index.o allpairs.o ppjoin.o ppjoinplus.o
PPJOINPLUSOBJECTS = main.o io.o index.o ppjoinplus.o 
PPJOINOBJECTS = main.o io.o index.o ppjoin.o 
ALLPAIRSOBJECTS = main.o io.o index.o allpairs.o

all: $(PPJOINPLUSOBJECTS) $(PPJOINOBJECTS) $(ALLPAIRSOBJECTS)
	${COMPILER} ${CCFLAGS} -o ppjoinplus $(PPJOINPLUSOBJECTS)
	${COMPILER} ${CCFLAGS} -o ppjoin $(PPJOINOBJECTS)
	${COMPILER} ${CCFLAGS} -o allpairs $(ALLPAIRSOBJECTS)

ppjoinplus: $(PPJOINPLUSOBJECTS)
	${COMPILER} ${CCFLAGS} -o ppjoinplus $(PPJOINPLUSOBJECTS)

ppjoin: $(PPJOINOBJECTS)
	${COMPILER} ${CCFLAGS} -o ppjoin $(PPJOINOBJECTS)

allpairs: $(ALLPAIRSOBJECTS)
	${COMPILER} ${CCFLAGS} -o allpairs $(ALLPAIRSOBJECTS)

main.o: main.cpp
	${COMPILER} ${CCFLAGS} -c main.cpp

io.o: io.cpp
	${COMPILER} ${CCFLAGS} -c io.cpp

index.o: index.cpp
	${COMPILER} ${CCFLAGS} -c index.cpp

ppjoinplus.o: ppjoinplus.cpp
	${COMPILER} ${CCFLAGS} ${CCINCLUDES} -c ppjoinplus.cpp

ppjoin.o: ppjoin.cpp
	${COMPILER} ${CCFLAGS} ${CCINCLUDES} -c ppjoin.cpp

allpairs.o: allpairs.cpp
	${COMPILER} ${CCFLAGS} ${CCINCLUDES} -c allpairs.cpp

clean:
	rm -f ppjoinplus ppjoin allpairs *.o
