Versions of python and related pkgs:
	I used the following versions, but it may work for others.
	
	 python 2.7
	 
	 cython 0.19.1
	 
	 numpy 1.8.0
 
Instructions:


In the terminal: 

cd cy/rodrigues/

make

cd ../rot_btwn_two_vecs/

make

Now, in the main directory, try:
python simple.py

If you get no errors, it means that the code is working. 
Inside simple.py, at the end of the file
you will see how to use the code. 
The important thing is the call:
Simple.compute_deformations(Xs,Ys,Rxs,Rys,As,Ss,Qs)    

That's the only function you need to use.

Please let me know if you get stuck.   








