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

using this quartet one can calculate the distances d(Ia) ; d(Ib) and d(ab) for which the sum will be equal to d(IJ). These distances are calculated using the formulae :

 * d(Ia) = 1/2 * ( d(IK) + d(IL) - d(KL) )


  
## Build Treemputer

```
git clone https://github.com/mbruto/Treemputer
cd Treemputer
make
```
