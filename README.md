README.md: Code for computing triangle deformations
--------------------------

Author: Oren Freifeld.
Email: freifeld.csail.mit.edu

This code repository has become publicly available: 03/20/2015.

Usage of this code repository is subject to the terms specified in the license.txt file.
_____________________________________________________________________________



It is based on [Freifeld and Black, ECCV 2012].
For technical reasons, it is not the same code we used in that paper. However, it produces the same deformations.

At the moment, note this is only a partial version of what we presented in the paper: at the moment, the only thing it does is taking pairs of traingles and compute the "Q" matrices and the corresponding (R,A,S) decompositions (see paper for details). 

Remark (03/20/2015): 
We intend to add other parts of the code (e.g., mesh synthesis, Lie algebra, exp/log, etc.) soon.

Requirements:

	Cython, numpy, scipy, opencv

Tested on 

	 Ubuntu 12.04 

	 python 2.7
	 
	 cython 0.19.1
	 
	 numpy 1.8.0
	 
	 opencv 2.4
	 
 
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






