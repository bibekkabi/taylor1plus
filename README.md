# taylor1plus
Zonotope abstract domain

Taylor1+ [1] is an implementation of affine sets based or zontopes abstract domain [2] in [Apron](http://apron.cri.ensmp.fr/library/).
It was originally implemented by Khalil Ghorbal. The link to taylor1+ implementation by Khalil is: https://github.com/kghorbal/taylor1plus. 
Since taylor1+ is an Apron complaint library it contains all the basic operations required by Apron. Additionally, it defined the meet operation of two affine sets as a logical product of standard zonotopes and boxes [3]. This meet operation is a functional interpretation of the intersection of a zonotope with a guard. 

In the current version of taylor1+ few new operations have been included. They are:

1. The **inclusion test (t1p_is_leq)**: to test if one zonotope is enclosed inside another. 
   It defines a cheaper test based on the result of Lemma 3 & 4 of [3]. 
   
2. The **intersection test (t1p_is_intersect)**: to test if one zonotope collides with another. 
   Earlier, this operation didn't exist in Apron. So, it is mostly likely that one may not find the same function implemented for    other abstract domains like boxes, octagons and polyhedra

3. The **vertex representation of a zonotope (t1p_to_generator_array)**: to compute all the vertices of a zonotope.
   We enumerate the vertices of a zonotope as sign vectors of the so-called hyperplane arrangement [3] corresponding to a zonotope
   This function is very useful if someone wants to plot a zonotope. One can a write a SVG function and call this fucntion to access the vertices.  
   
   All the above functions are inside **t1p_constructor.c**
   
4. The **split function for zonotopes (t1p_tilings)**: to split a zonotope into a set of parallelotopes such that the union of the  parallelotopes is the original zonotope. Note: there is no overlap between any two parallelotopes [SWIM2016](https://swim2016.sciencesconf.org/data/pages/Kabi_Goubault_Putot.pdf)[Poster@Marktoberdorf](https://asimod.in.tum.de/2017/posters/Kabi_Bibek.pdf).
  This function is inside **t1p_tilings.c**
  This is the first time a split operation is inclued in Apron. 
  However, this function is only restricted to taylor1+. 
  We wanted to call this function from an OCaml protoptype analyzer [Prototype_analyzerwithApron](https://github.com/bibekkabi/Prototype_analyzerwithApron). The OCaml bindings of this function can be found inside **t1p_idl.c**
  
5. The **meet operation between two zonotopes (t1p_meet)**: a new geometric meet operation on zonotopes. 
   This function is inside **t1p_meetjoin.c**
   
**Note: The steps to install Apron can be found in Pointers_apron.pdf** 
   
 
# References

[1] Ghorbal, K., Goubault, E., & Putot, S. (2009, June). The zonotope abstract domain taylor1+. In International Conference on Computer Aided Verification (pp. 627-633). Springer, Berlin, Heidelberg.

[2] Goubault, E., & Putot, S. (2006, August). Static analysis of numerical algorithms. In International Static Analysis Symposium (pp. 18-34). Springer, Berlin, Heidelberg.

[3] Ghorbal, K., Goubault, E., & Putot, S. (2010, July). A logical product approach to zonotope intersection. In International Conference on Computer Aided Verification (pp. 212-226). Springer, Berlin, Heidelberg.

[4] Goubault, E., & Putot, S. (2015). A zonotopic framework for functional abstractions. Formal Methods in System Design, 47(3), 302-360.
