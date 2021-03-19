
default: libc/Library_iaf.c
	python3 setup.py build_ext --inplace && mv Library_iaf.*.so ./Library_iaf.so && rm -Rf build



clean :
	rm *.so


