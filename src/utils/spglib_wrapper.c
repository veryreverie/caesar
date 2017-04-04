#include "spglib.h"

void spglib_get_no_symmetries
(
  // Inputs.
        double (* lattice )[3][3],
        double (* position)[][3],
  const int    (* types   )[],
  const int     * num_atom,
  const double  * symprec,
  // Outputs.
        int     * n_operations
)
{
  SpglibDataset spg_dataset = *spg_get_dataset
  (
    *lattice,
    *position,
    *types,
    *num_atom,
    *symprec
  );
  
  *n_operations = spg_dataset.n_operations;
}

void spglib_get_symmetries
(
  // Inputs.
        double (* lattice )[3][3],
        double (* position)[][3],
  const int    (* types   )[],
  const int     * num_atom,
  const double  * symprec,
  // Outputs.
        int    (* rotations)[][3][3],
        double (* translations)[][3]
)
{
  SpglibDataset spg_dataset = *spg_get_dataset
  (
    *lattice,
    *position,
    *types,
    *num_atom,
    *symprec
  );
  
  for (int i=0;i<spg_dataset.n_operations;++i)
  {
    for (int j=0;j<3;++j)
    {
      for (int k=0;k<3;++k)
      {
        (*rotations)[i][j][k] = spg_dataset.rotations[i][j][k];
      }
      (*translations)[i][j] = spg_dataset.translations[i][j];
    }
  }
}
