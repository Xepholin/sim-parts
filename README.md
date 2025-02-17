# Simulation de Dynamique Moléculaire

Ce projet est une simulation de dynamique moléculaire implémentée en C++. Il simule le comportement de particules sous l'influence des forces du potentiel de Lennard-Jones. La simulation prend en charge les conditions aux limites périodiques et peut utiliser l'algorithme de la liste de Verlet pour des calculs de force efficaces.

## Fonctionnalités

- **Potentiel de Lennard-Jones** : Simule les interactions entre les particules en utilisant le potentiel de Lennard-Jones.
- **Conditions aux Limites Périodiques** : Assure que les particules interagissent à travers les limites de manière périodique.
- **Liste de Verlet** : Optimise les calculs de force en maintenant une liste des voisins proches pour chaque particule.
- **Intégration de Verlet** : Utilise l'algorithme de Verlet pour l'intégration des équations du mouvement.
- **Contrôle de la Température** : Permet de maintenir la température du système à une valeur souhaitée grâce à un facteur de correction.
- **Sortie au Format PDB** : Génère des fichiers de sortie au format PDB pour visualiser les configurations des particules.

## Structure du Projet

- **main.cpp** : Point d'entrée du programme, gère les arguments de la ligne de commande et lance la simulation.
- **part.cpp** : Implémentation des méthodes de la classe `Simulator`, y compris les calculs de force, l'intégration et la gestion des conditions aux limites.
- **tools.cpp** : Contient des fonctions utilitaires pour les calculs statistiques et d'initialisation.
- **const.hpp** : Définit les constantes physiques et les paramètres de simulation.
- **part.hpp** : Déclaration de la classe `Simulator` et des structures de données utilisées.

## Compilation

Le projet utilise CMake pour la compilation. Pour compiler le projet, exécutez les commandes suivantes :

```bash
mkdir build
cd build
cmake ..
make
```
## Utilisation

Le programme prend plusieurs arguments en ligne de commande pour configurer la simulation :
```bash
./spart n_particles n_sym t_min t_max dt T0 gamma use_verlet input_file [output_file]
```
- n_particles : Nombre de particules dans la simulation.

- n_sym : Nombre de translations périodiques à considérer.

- t_min : Temps initial de la simulation.

- t_max : Temps final de la simulation.

- dt : Pas de temps pour l'intégration.

- T0 : Température cible du système.

- gamma : Facteur de correction pour le contrôle de la température.

- use_verlet : Utiliser la liste de Verlet (1 pour oui, 0 pour non).

- input_file : Fichier d'entrée contenant les positions initiales des particules.

- output_file : Fichier de sortie pour enregistrer les configurations (optionnel, par défaut `../output.pdb`).

## Exemple

Pour lancer une simulation avec 1000 particules, en utilisant la liste de Verlet, et enregistrer les résultats dans `output.pdb` :
```bash
./spart 1000 27 0 100 0.01 300 0.1 1 initial_positions.txt output.pdb
```
## Visualisation des Résultats

Les fichiers de sortie au format PDB peuvent être visualisés à l'aide de logiciels de visualisation moléculaire tels que VMD.