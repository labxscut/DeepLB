CC=g++
CFLAGS= -w -O3 -static 
#CFLAGS= -Wall -O3 -static 
getoptDIR=./my_getopt-1.5
INCLUDE= -I/home/cmb-panasas2/wel/lib/libs_cpp/boost_1_52_0/include -I/home/cmb-panasas2/wel/lib/libs_c/gsl-1.15/include
OBJECTS1= wig_binning.o data_types.o utils.o my_getopt.o
OBJECTS2= reads_binning.o sam.o data_types.o utils.o my_getopt.o
OBJECTS3= calc_reads_prob.o matrix.o utils.o data_types.o my_getopt.o
OBJECTS4= em.o matrix.o utils.o data_types.o my_getopt.o
OBJECTS5= calc_reads_prob_grid_methy.o matrix.o utils.o data_types.o my_getopt.o
LIBS= -L/home/cmb-panasas2/wel/lib/libs_cpp/boost_1_52_0/lib -L/home/cmb-panasas2/wel/lib/libs_c/gsl-1.15/lib -L/home/cmb-panasas2/wel/lib/libs_c/zlib-1.2.8/lib -lboost_iostreams -lz -lm -lgsl -lgslcblas
PROG1=wig_binning
PROG2=reads_binning
PROG3=calc_reads_prob
PROG4=em
PROG5=calc_reads_prob_grid_methy

all: ${PROG1} ${PROG2} ${PROG3} ${PROG4} ${PROG5}

${PROG1}: ${OBJECTS1}
	${CC} ${CFLAGS} ${OBJECTS1} ${LIBS} -o ${PROG1}

${PROG2}: ${OBJECTS2}
	${CC} ${CFLAGS} ${OBJECTS2} ${LIBS} -o ${PROG2}

${PROG3}: ${OBJECTS3}
	${CC} ${CFLAGS} ${OBJECTS3} ${LIBS} -o ${PROG3}

${PROG4}: ${OBJECTS4}
	${CC} ${CFLAGS} ${OBJECTS4} ${LIBS} -o ${PROG4}

${PROG5}: ${OBJECTS5}
	${CC} ${CFLAGS} ${OBJECTS5} ${LIBS} -o ${PROG5}

%.o: %.cpp %.h
	${CC} ${CFLAGS} ${INCLUDE} -c $< -o $@

my_getopt.o: $(getoptDIR)/my_getopt.c $(getoptDIR)/my_getopt.h $(getoptDIR)/getopt.h
	cd $(getoptDIR)
	$(CC) $(Flags) -I$(getoptDIR) -c $(getoptDIR)/my_getopt.c -o my_getopt.o

# --- remove binary and executable files
clean:
	rm -f ${OBJECTS1} ${OBJECTS2} ${OBJECTS3} ${OBJECTS4} ${OBJECTS5}
	rm -f ${PROG1} ${PROG2} ${PROG3} ${PROG4} ${PROG5}

copy:
	/bin/cp -f ${PROG1} ../bin/
	/bin/cp -f ${PROG2} ../bin/
	/bin/cp -f ${PROG3} ../bin/
	/bin/cp -f ${PROG4} ../bin/
	/bin/cp -f ${PROG5} ../bin/

