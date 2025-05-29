file(REMOVE_RECURSE
  "libscifor.a"
  "libscifor.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang Fortran)
  include(CMakeFiles/scifor.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
