# Treemputer

## Imputation of missing distances in Supertree matrices

Treemputer is a method allowing the reconstruction of Supertree with branch length. It takes multiple source trees as input and perform the following steps :

  * Create a sparse distance matrix for all species present in the source tree. Distances of the matrix correspond to the mean patristic distance.

  * Impute the missing distances based on known distances.

  * Output a complete matrix if the necessary conditions are respected.

This method is highly inspired from previous work of Alain Guénoche and colleagues and is based on (i) the four-point condition of trees and (ii) distance imputation using quartet distances.

Consider two species I and J with unknown distance in larger distance matrix. The program selects two additional species K and L for which all distances are known. For this example, let's consider that I;K and J;L are closely related (have a small distance) while I;L and J;K are distantly related (have a long distance). We obtain the following quartet :

```
I                   J    I/J/K/L are the species considered in the quartet with :
 \                 /        
  \               /         - d(IJ) unknown
   \a___________b/          - d(IK); d(IL); d(JK); d(JL) and d(KL) known   
   /             \        
  /               \      a and b are the internal nodes of the quartet
 /                 \
K                   L
```

Using this quartet one can calculate the distances d(Ia) ; d(Jb) and d(ab) for which the sum will be equal to d(IJ). These distances are calculated using the formulae :

 * d(Ia) = 1/2 * ( d(IK) + d(IL) - d(KL) )
 * d(Jb) = 1/2 * ( d(JL) + d(JK) - d(KL) )
 * d(ab) = 1/2 * ( d(IL) + d(JK) - d(IK) - d(KL) )
 * d(IJ) = d(Ia) + d(Jb) + d(ab)

Previous attempts have been made to perform imputation of missing distances but our method comes with multiple heuristics that ameliorate unknown distance imputation.

First, 

## Limitations

In general, different kind of source trees can be given to reconstruct a Supertree. Classically, each source tree has been reconstructed with different loci and are further amalgamated. **Our method has not been designed to work on this kind of datasets and users should not use Treemputer in that context**. Other programs have been specifically designed to use these data such as ASTRAL. **As an input, users should provide overlapping source trees that have been inferred on a similar set of markers (either single marker or multiple marker) as distances are important for the imputation process.**

**Another important limitation of the approach is the impossibility to impute missing distances of sister nodes. This means that users are awaited to give source trees covering closely related species in order to take full advantage of the method.**

  
## Build Treemputer

```
git clone https://github.com/mbruto/Treemputer
cd Treemputer
make
```
