file(REMOVE_RECURSE
  "../../build_x64/Debug/lib/libSWIFT.a"
  "../../build_x64/Debug/lib/libSWIFT.pdb"
  "SWIFT_autogen"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/SWIFT.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()