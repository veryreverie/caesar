submodule (caesar_physical_constants_module) caesar_physical_constants_submodule
  use caesar_structure_module
  
  type(String), allocatable :: ATOMIC_SYMBOLS(:)
contains

module procedure initialise_atomic_symbols
  if (.not. allocated(ATOMIC_SYMBOLS)) then
    ATOMIC_SYMBOLS = [ str('H'),  &
                     & str('He'), &
                     & str('Li'), &
                     & str('Be'), &
                     & str('B'),  &
                     & str('C'),  &
                     & str('N'),  &
                     & str('O'),  &
                     & str('F'),  &
                     & str('Ne'), &
                     & str('Na'), &
                     & str('Mg'), &
                     & str('Al'), &
                     & str('Si'), &
                     & str('P'),  &
                     & str('S'),  &
                     & str('Cl'), &
                     & str('Ar'), &
                     & str('K'),  &
                     & str('Ca'), &
                     & str('Sc'), &
                     & str('Ti'), &
                     & str('V'),  &
                     & str('Cr'), &
                     & str('Mn'), &
                     & str('Fe'), &
                     & str('Co'), &
                     & str('Ni'), &
                     & str('Cu'), &
                     & str('Zn'), &
                     & str('Ga'), &
                     & str('Ge'), &
                     & str('As'), &
                     & str('Se'), &
                     & str('Br'), &
                     & str('Kr'), &
                     & str('Rb'), &
                     & str('Sr'), &
                     & str('Y'),  &
                     & str('Zr'), &
                     & str('Nb'), &
                     & str('Mo'), &
                     & str('Tc'), &
                     & str('Ru'), &
                     & str('Rh'), &
                     & str('Pd'), &
                     & str('Ag'), &
                     & str('Cd'), &
                     & str('In'), &
                     & str('Sn'), &
                     & str('Sb'), &
                     & str('Te'), &
                     & str('I'),  &
                     & str('Xe'), &
                     & str('Cs'), &
                     & str('Ba'), &
                     & str('La'), &
                     & str('Ce'), &
                     & str('Pr'), &
                     & str('Nd'), &
                     & str('Pm'), &
                     & str('Sm'), &
                     & str('Eu'), &
                     & str('Gd'), &
                     & str('Tb'), &
                     & str('Dy'), &
                     & str('Ho'), &
                     & str('Er'), &
                     & str('Tm'), &
                     & str('Yb'), &
                     & str('Lu'), &
                     & str('Hf'), &
                     & str('Ta'), &
                     & str('W'),  &
                     & str('Re'), &
                     & str('Os'), &
                     & str('Ir'), &
                     & str('Pt'), &
                     & str('Au'), &
                     & str('Hg'), &
                     & str('Tl'), &
                     & str('Pb'), &
                     & str('Bi'), &
                     & str('Po'), &
                     & str('At'), &
                     & str('Rn'), &
                     & str('Fr'), &
                     & str('Ra'), &
                     & str('Ac'), &
                     & str('Th'), &
                     & str('Pa'), &
                     & str('U'),  &
                     & str('Np'), &
                     & str('Pu'), &
                     & str('Am'), &
                     & str('Cm'), &
                     & str('Bk'), &
                     & str('Cf'), &
                     & str('Es'), &
                     & str('Fm'), &
                     & str('Md'), &
                     & str('No'), &
                     & str('Lr'), &
                     & str('Rf'), &
                     & str('Db'), &
                     & str('Sg'), &
                     & str('Bh'), &
                     & str('Hs'), &
                     & str('Mt'), &
                     & str('Ds'), &
                     & str('Rg'), &
                     & str('Cn'), &
                     & str('Nh'), &
                     & str('Fl'), &
                     & str('Mc'), &
                     & str('Lv'), &
                     & str('Ts'), &
                     & str('Og')  ]
  endif
end procedure

module procedure atomic_number_to_symbol
  call initialise_atomic_symbols()
  
  if (input<1 .or. input>size(ATOMIC_SYMBOLS)) then
    call print_line(ERROR//': '//input//' is not a valid atomic number.')
    call err()
  endif
  
  output = ATOMIC_SYMBOLS(input)
end procedure

module procedure atomic_symbol_to_number
  call initialise_atomic_symbols()
  
  output = first(ATOMIC_SYMBOLS==input, default=0)
end procedure
end submodule
