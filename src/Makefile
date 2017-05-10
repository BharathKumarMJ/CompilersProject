INC = -I./include
HEADERDIR = ./include

all: mmio.o spmm.o read_matrix.o merge_sort.o  main.o 
	gcc -O3 mmio.o spmm.o read_matrix.o merge_sort.o main.o -o main

main.o: main.c 
	gcc -O3 $(INC) -c main.c 
mmio.o: mmio.c 
	gcc -O3 $(INC) -c mmio.c  
spmm.o: spmm.c 
	gcc -O3 $(INC) -c spmm.c
read_matrix.o: read_matrix.c
	gcc -O3 $(INC) -c read_matrix.c 
merge_sort.o: merge_sort.c 
	gcc -O3 $(INC) -c merge_sort.c 

clean:
	rm -rf *.o
	rm -rf main