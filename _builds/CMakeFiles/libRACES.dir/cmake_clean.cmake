file(REMOVE_RECURSE
  ".0.1.0"
  "RACES.0.1.0.dylib"
  "RACES.dylib"
  "RACES.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang CXX)
  include(CMakeFiles/libRACES.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
