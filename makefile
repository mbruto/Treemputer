
default: libc/Library_iaf.c
	gcc -c libc/Library_iaf.c -o Library_iaf.so -fPIC `python3-config --cflags`



clean :
	rm *.so


