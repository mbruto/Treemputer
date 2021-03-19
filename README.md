# Treemputer

## Imputation of missing distances in Supertree matrices

Treemputer is a mathod allowing the reconstruction of Supertree with branch length. It takes multiple source trees as input and perform the following steps :

    - Create a sparse distance matrix for all species present in the source tree. Distances of the matrix correspond to the mean patristic distance.
    - Impute the missing distances based on known distances.
    - Output a complete matrix if the necessary conditions are respected.

## Build Treemputer

```
git clone https://github.com/mbruto/Treemputer
cd Treemputer
make
```
