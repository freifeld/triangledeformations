
 
PYX_FILENAME = rot_btwn_two_vecs

all:
	python setup.py build_ext --inplace
test: 
	python ~/Dropbox/code/pymonkey/pose/demo_pose_Isaac_into_the_scene.py 

annotate:
	cython -a $(PYX_FILENAME).pyx  
	firefox $(PYX_FILENAME).html
clean:
	rm -rf build
	rm *.c *.so
	 
