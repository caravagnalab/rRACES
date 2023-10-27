file(REMOVE_RECURSE
  "libRACES.a"
  "libRACES.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang CXX)
  include(CMakeFiles/libRACES_static.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
