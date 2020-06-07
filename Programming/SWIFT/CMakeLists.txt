cmake_minimum_required(VERSION 3.15)

# Set the project name to the name of the directory replacing any spaces by _
get_filename_component(PROJECT_NAME ${CMAKE_CURRENT_LIST_DIR} NAME)
string(REPLACE " " "_" ProjectId ${PROJECT_NAME})

message(${PROJECT_NAME})
project(${PROJECT_NAME})

set(SOURCES
        HestonCalibrator.cpp
        distributions.h
        option_contracts.cpp
        option_contracts.h
        swift.h
        my_math.cpp
        my_math.h
        known_distribution_contract_combinations.h c_m_calculators.cpp c_m_calculators.h swift.cpp distributions.cpp)

# For any nested file, replace "/" by  "\"
foreach(source ${SOURCES})
    get_filename_component(source_path "${source}" PATH)
    string(REPLACE "/" "\\" source_path_msvc "${source_path}")
    source_group("${source_path_msvc}" FILES "${source}")
endforeach()

add_library(${PROJECT_NAME} STATIC ${SOURCES})
set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(${PROJECT_NAME} PUBLIC ${BOOST_STATIC_LIBRARIES})

add_subdirectory(UnitTests)
add_subdirectory(Z_Exec_CuiExperiment)
add_subdirectory(Z_Exec_CuiExperimentOriginal)