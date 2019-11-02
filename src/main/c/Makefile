# JAVA_HOME=
#JAVA_HOME=/usr/lib/jvm/java-8-openjdk-amd64/

# if you need to rebuild JNI , it's needed to export JAVA_HOME variable
#
# I'm doing it as:
# export JAVA_HOME=/usr/lib/jvm/java-8-openjdk-amd64
# to point to Java-8 files
#
# you will also need superlu-mt headers. Makefile is written so that sources are expected in
# JKLU/src/superlu-mt

all:   libsuperLUWrapper.so

libsuperLUWrapper.so: edu_ufl_cise_klu_wrapper_SuperLUWrapper.o libsuperlu_mt_PTHREAD.a
	gcc -shared -o $@ edu_ufl_cise_klu_wrapper_SuperLUWrapper.o  libsuperlu_mt_PTHREAD.a \
-Wl,--no-whole-archive  -Wl,-rpath,/usr/lib/x86_64-linux-gnu/ -L/usr/lib/x86_64-linux-gnu/  -lpthread -lm -lblas

edu_ufl_cise_klu_wrapper_SuperLUWrapper.o: edu_ufl_cise_klu_wrapper_SuperLUWrapper.c
	gcc -shared -I"$(JAVA_HOME)/include" -I"$(JAVA_HOME)/include/linux" -D_LONGINT -I../../superlu-mt/SRC -fPIC -O3 -c $^ -o $@

.PHONY: clean
clean:
	rm -rf *.o libsuperLUWrapper.so
