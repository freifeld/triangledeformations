
TEST_PY_FILENAE = ~/Dropbox/code/pyHumanBody/debug_log.py
PYX_FILENAME = SO3so3

all:
	python setup.py build_ext --inplace
test:   
	python $(TEST_PY_FILENAE)
annotate:
	cython -a $(PYX_FILENAME).pyx  
	firefox $(PYX_FILENAME).html
clean:
	rm -rf build
	rm *.c *.so
	rm *.c *.so
	 
