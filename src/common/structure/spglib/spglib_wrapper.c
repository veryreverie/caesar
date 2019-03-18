/*
 * See corresponding interfaces in spglib_wrapper.f90 for comments.
 */

#include <stdbool.h> // bool
#include <stddef.h>  // NULL
#include <stdio.h>   // FILE
#include <string.h>  // strcpy

#include "spglib.h"

void spglib_calculate_symmetries
(
        double       (*  lattice )[3][3],
        double       (*  position)[][3],
  const int          (*  types   )[],
  const int           *  num_atom,
  const double        *  symprec,
        bool          *  success,
        SpglibDataset ** spg_dataset,
        int           *  n_operations
)
{
  *spg_dataset = spg_get_dataset
  (
    *lattice,
    *position,
    *types,
    *num_atom,
    *symprec
  );
  
  if (*spg_dataset == NULL)
  {
    *success = false;
  }
  else
  {
    *success = true;
    *n_operations = (*spg_dataset)->n_operations;
  }
}

void spglib_retrieve_symmetries
(
  const SpglibDataset  ** spg_dataset,
        int               spacegroup_number,
        char           *  international_symbol,
        double            transformation[3][3],
        double            origin_shift[3],
        int               n_operations,
        int           (*  tensors)[][3][3],
        double        (*  translations)[][3],
        int               n_atoms,
        char           *  pointgroup_symbol
)
{
  spacegroup_number = (*spg_dataset)->spacegroup_number;
  strcpy(international_symbol, (*spg_dataset)->international_symbol);
  for (int i=0;i<3;++i)
  {
    for (int j=0;j<3;++j)
    {
      // Transpose transformation from C layout to Fortran layout.
      transformation[j][i] = (*spg_dataset)->transformation_matrix[i][j];
    }
  }
  for (int i=0;i<3;++i)
  {
    origin_shift[i] = (*spg_dataset)->origin_shift[i];
  }
  n_operations = (*spg_dataset)->n_operations;
  for (int i=0;i<n_operations;++i)
  {
    for (int j=0;j<3;++j)
    {
      for (int k=0;k<3;++k)
      {
        // Transpose tensors from C layout to Fortran layout.
        (*tensors)[i][j][k] = (*spg_dataset)->rotations[i][k][j];
      }
      (*translations)[i][j] = (*spg_dataset)->translations[i][j];
    }
  }
  n_atoms = (*spg_dataset)->n_atoms;
  strcpy(pointgroup_symbol, (*spg_dataset)->pointgroup_symbol);
}

void drop_spg_dataset
(
  SpglibDataset ** spg_dataset
)
{
  spg_free_dataset(*spg_dataset);
}

int spglib_standardize_cell
(
        double (*  lattice)[3][3],
        double (*  position)[][3],
        int    (*  types)[],
  const int     *  num_atom,
  const int     *  to_primitive,
  const int     *  no_idealize,
  const double  *  symprec
)
{
  return spg_standardize_cell
  (
    *lattice,
    *position,
    *types,
    *num_atom,
    *to_primitive,
    *no_idealize,
    *symprec
  );
}
