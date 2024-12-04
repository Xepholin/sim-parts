# Code de Simulation Microscopique

## Compilation

### Pre-requis
- CMake 3.16+
- C17 conforming compiler

### Build
```sh
cmake -S . -B <BUILD_DIR> [CMAKE_FLAGS ...]
cmake --build <BUILD_DIR> [-j]
```

### Run
```sh
<BUILD_DIR>/spart [PARTICULES_FILE_PATH]
```
## À propos

**Objectif : développer un code de simulation microscopique basé sur un potentiel de Lennard-Jones pour étudier un fluide de particules homogène (i.e. toutes les particules sont identiques).**

- **Langage: au choix**
- **Processeur visé: au choix.**

Le fichier ``particule.xyz`` contient les coordonnées cartésiennes (x,y,z) d’un système de 1 000 particules
identiques, avec le format :

    > 0 1
    > 2 x y z
    > 2 x y z

La première ligne est un commentaire, dans les suivantes ‘ 2 ‘ est le type de particule, pour notre sujet sans importance. Ecrire un code qui permet de lire ce fichier et de stocker les coordonnées x,y,z dans un (ou des) vecteurs ou tableaux. On introduira les paramètre N_particules_total = 1000 et N_particules_local, le dernier permettant de réaliser des calculs sur un sous ensemble des N_particules_local premières particules (N_particules_local < N_particules_total).


Sur la base des coordonnées (x,y,z) des 1000 particules, calculer l’énergie microscopique de ce système et les forces agissant sur chacune des particules dans le cas d’un potentiel de Lennard Jones :

$$
U^{LJ} = 4 \sum_{i=1}^N \sum_{j>i}^N \varepsilon^* \left[ \left(\frac{r^*}{r_{ij}^*}\right)^{12} - 2 \left(\frac{r^*}{r_{ij}^*}\right)^6 \right] = 4 \sum_{i=1}^N \sum_{j>i}^N u_{ij}
$$

On prendra $r^* = 3.0$ et $\varepsilon^* = 0.2$.

Rappel pour un terme de Lennard Jones

$$
U^{LJ} = 4 \sum_{i=1}^N \sum_{j>i}^N \varepsilon^* \left[ \left(\frac{r^*}{r_{ij}^*}\right)^{12} - 2 \left(\frac{r^*}{r_{ij}^*}\right)^6 \right] = \sum_{i=1}^N \sum_{j>i}^N u_{ij}
$$

Les forces des particules se calculent à partir du gradient analytique associé aux fonctions élémentaires $u_{ij}$

$$
\frac{\partial u_{ij}}{\partial x_i} = -4 \varepsilon^* \left[ 12 \left(\frac{r^*}{r_{ij}}\right)^{13} - 2 \times 6 \left(\frac{r^*}{r_{ij}}\right)^7 \right] \times \frac{\partial r_{ij}}{\partial x_i} = -48 \varepsilon_{ij}^* \left[ \left(\frac{r^*}{r_{ij}}\right)^{14} - \left(\frac{r^*}{r_{ij}}\right)^8 \right] \times (x_i - x_j)
$$

Comme pour les coordonnées (x,y,z), les forces (fx,fy,fz) seront stockées dans un (des) vecteur(s) ou tableau(x).