file(REMOVE_RECURSE
  "../../build_x64/Test/lib/libFFTW3.a"
  "../../build_x64/Test/lib/libFFTW3.pdb"
  "FFTW3_autogen"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/FFTW3.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
