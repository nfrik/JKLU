# JAVA_HOME=
#JAVA_HOME=/usr/lib/jvm/java-8-openjdk-amd64/

all:   libsuperLUWrapper.so

libsuperLUWrapper.so: edu_ufl_cise_klu_wrapper_SuperLUWrapper.o libsuperlu_mt_PTHREAD.a
	gcc -shared -o $@ edu_ufl_cise_klu_wrapper_SuperLUWrapper.o  libsuperlu_mt_PTHREAD.a \
-Wl,--no-whole-archive  -Wl,-rpath,/usr/lib/x86_64-linux-gnu/ -L/usr/lib/x86_64-linux-gnu/  -lpthread -lm -lblas

edu_ufl_cise_klu_wrapper_SuperLUWrapper.o: edu_ufl_cise_klu_wrapper_SuperLUWrapper.c
	gcc -I"$(JAVA_HOME)/include" -I"$(JAVA_HOME)/include/linux"  -I../../superlu-mt/SRC -fPIC -O2 -c $^ -o $@

.PHONY: clean
clean:
	rm -rf *.o libsuperLUWrapper.so
