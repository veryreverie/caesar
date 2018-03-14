/*
 * See corresponding interfaces in calculate_symmetry.f90 for comments.
 */

#include <stdbool.h> // bool
#include <stddef.h>  // NULL

#include "spglib.h"

void spglib_calculate_symmetries
(
  // Inputs.
        double       (*  lattice )[3][3],
        double       (*  position)[][3],
  const int          (*  types   )[],
  const int           *  num_atom,
  const double        *  symprec,
  // Outputs.
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
  // Inputs.
  const SpglibDataset ** spg_dataset,
  // Outputs.
        int          (* rotations)[][3][3],
        double       (* translations)[][3]
)
{
  for (int i=0;i<(*spg_dataset)->n_operations;++i)
  {
    for (int j=0;j<3;++j)
    {
      for (int k=0;k<3;++k)
      {
        // N.B. Transposes rotations from C layout to Fortran layout.
        (*rotations)[i][j][k] = (*spg_dataset)->rotations[i][k][j];
      }
      (*translations)[i][j] = (*spg_dataset)->translations[i][j];
    }
  }
}

void drop_spg_dataset
(
  // Input.
  SpglibDataset ** spg_dataset
)
{
  spg_free_dataset(*spg_dataset);
}
