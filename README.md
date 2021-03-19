# Treemputer

## Imputation of missing distances in Supertree matrices

Treemputer is a method allowing the reconstruction of Supertree with branch length. It takes multiple source trees as input and perform the following steps :

  * Create a sparse distance matrix for all species present in the source tree. Distances of the matrix correspond to the mean patristic distance.

  * Impute the missing distances based on known distances.

  * Output a complete matrix if the necessary conditions are respected.

This method is highly inspired from previous work of Alain Guénoche and colleagues and is based on (i) the additivity property of trees and (ii) distance imputation using quartet distances.

Consider two species I and J with unknown distance in larger distance matrix. The program selects two additional species K and L for which all distances are known. For this example, let's consider that I;K and J;L are closely related (have a small distance) while I;L and J;K are distantly related (have a long distance). We obtain the following quartet :

I                   J
 \                 /
  \               /
   \a___________b/
   /             \
  /               \
 /                 \
K                   L


  
## Build Treemputer

```
git clone https://github.com/mbruto/Treemputer
cd Treemputer
make
```
