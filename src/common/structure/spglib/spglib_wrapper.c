/*
 * See corresponding interfaces in calculate_symmetry.f90 for comments.
 */

#include <stdbool.h> // bool
#include <stddef.h>  // NULL

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
  const SpglibDataset ** spg_dataset,
        int           spacegroup_number,
        int           hall_number,
        char          international_symbol[11],
        char          hall_symbol[17],
        char          choice[6],
        double        transformation[3][3],
        double        origin_shift[3],
        int           n_operations,
        int           (* tensors)[][3][3],
        double        (* translations)[][3],
        int           n_atoms,
        char          pointgroup_symbol[6]
)
{
  spacegroup_number = (*spg_dataset)->spacegroup_number;
  hall_number = (*spg_dataset)->hall_number;
  international_symbol = (*spg_dataset)->international_symbol;
  hall_symbol = (*spg_dataset)->hall_symbol;
  choice = (*spg_dataset)->choice;
  for (int i=0;i<3;++i)
  {
    for (int j=0;j<3;++j)
    {
      // Transpose transformation from C layout to Fortran layout.
      transformation[j][i] = (*spg_dataset)->transformation_matrix[i][j];
    }
  }
  origin_shift = (*spg_dataset)->origin_shift;
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
  pointgroup_symbol = (*spg_dataset)->pointgroup_symbol;
}

void drop_spg_dataset
(
  SpglibDataset ** spg_dataset
)
{
  spg_free_dataset(*spg_dataset);
}
