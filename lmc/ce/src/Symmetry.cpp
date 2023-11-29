/**************************************************************************************************
 * Copyright (c) 2023-2023. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date: 7/18/23 4:41 PM                                                                         *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 9/27/23 11:29 AM                                                          *
 **************************************************************************************************/

#include "Symmetry.h"
#include "spglib.h"

void FindPointGroupSymmetry(const Config &config) {

  SpglibDataset *dataset;

  double lattice_c[3][3];
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      lattice_c[i][j] = config.GetBasis()(static_cast<int>(i), static_cast<int>(j));

  size_t num_atom = config.GetNumLattices();
  auto position = new double[num_atom][3];
  for (size_t i = 0; i < num_atom; ++i)
    for (size_t j = 0; j < 3; ++j)
      position[i][j] = config.GetRelativePositionOfLattice(i)(static_cast<int>(j));

  int *types = new int[num_atom];
  for (size_t i = 0; i < num_atom; i++) { types[i] = 1; }
  // SpglibDataset has to be freed after use.
  dataset = spg_get_dataset(lattice_c, position, types, static_cast<int>(num_atom), 1e-5);
  // delete[] position;
  // delete[] types;

  printf("International symbol: %s (%d)\n", dataset->international_symbol,
         dataset->spacegroup_number);
  printf("Hall symbol:   %s\n", dataset->hall_symbol);
  printf("Wyckoff letters:\n");

  printf("\n");
  printf("Equivalent atoms:\n");
  for (size_t i = 0; i < dataset->n_atoms; i++) {
    printf("%zu -> %d\n", i, dataset->equivalent_atoms[i]);
  }
  printf("Space group operations:\n");
  for (size_t i = 0; i < dataset->n_operations; i++) {
    printf("--- %zu ---\n", i + 1);
    for (size_t j = 0; j < 3; j++) {
      printf("%2d %2d %2d\n", dataset->rotations[i][j][0],
             dataset->rotations[i][j][1], dataset->rotations[i][j][2]);
    }
    printf("%f %f %f\n", dataset->translations[i][0],
           dataset->translations[i][1], dataset->translations[i][2]);
  }
  printf("Point group symbol: %s\n", dataset->pointgroup_symbol);

  // Deallocate SpglibDataset, otherwise induce memory leak.
  spg_free_dataset(dataset);
}
