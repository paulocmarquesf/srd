RM = del

all:
	 gcc -ansi -O3 -o mcmc -lm lib/util.c lib/linalg.c lib/prob.c mcmc.c
	
clean:
	 $(RM) mcmc.exe
