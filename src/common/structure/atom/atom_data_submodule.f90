submodule (caesar_atom_data_module) caesar_atom_data_submodule
  use caesar_atom_module
contains

module procedure new_AtomData
  this%lattice_ = lattice
  this%recip_lattice_ = recip_lattice
  
  this%species_ = basic_atom%species
  this%mass_ = basic_atom%mass
  
  call this%set_cartesian_position(basic_atom%cartesian_position)
  
  this%id_ = id
  this%prim_id_ = prim_id
  this%rvec_id_ = rvec_id
end procedure

module procedure new_BasicAtom_AtomData
  output = BasicAtom(this%species(), this%mass(), this%cartesian_position())
end procedure

module procedure species
  output = this%species_
end procedure

module procedure mass
  output = this%mass_
end procedure

module procedure fractional_position
  output = this%fractional_position_
end procedure

module procedure cartesian_position
  output = this%cartesian_position_
end procedure

module procedure id
  output = this%id_
end procedure

module procedure prim_id
  output = this%prim_id_
end procedure

module procedure rvec_id
  output = this%rvec_id_
end procedure

module procedure set_species
  this%species_ = species
end procedure

module procedure set_mass
  this%mass_ = mass
end procedure

module procedure set_fractional_position
  this%fractional_position_ = fractional_position
  this%cartesian_position_ = transpose(this%lattice_) * fractional_position
end procedure

module procedure set_cartesian_position
  this%cartesian_position_ = cartesian_position
  this%fractional_position_ = this%recip_lattice_ * cartesian_position
end procedure
end submodule
