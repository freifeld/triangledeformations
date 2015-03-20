README.md
--------------------------

_Code for computing triangle deformations._


OS/packages/versions/OS/versions requirements: tested on 

	 Ubuntu 12.04

	 python 2.7
	 
	 cython 0.19.1
	 
	 numpy 1.8.0
 
Instructions:

You first need to compile some cython files.

	cd <where-you-put-this-repo>/triangledeformations/cy/rodrigues/

	make

	cd c<where-you-put-this-repo>/triangledeformations/rot_btwn_two_vecs/

	make

Then:

	cd <where-you-put-this-repo>/triangledeformations/

	python simple.py

If you get no errors, it means that the code is working. 
Inside simple.py, at the end of the file, you will see how to use the code. 
The important thing is the call:

	Simple.compute_deformations(Xs,Ys,Rxs,Rys,As,Ss,Qs)    

That's probably the only function you will need to use.

Oren Freifeld 





