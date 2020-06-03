cmake_minimum_required(VERSION 3.15)

# Set the project name to the name of the directory replacing any spaces by _
get_filename_component(PROJECT_NAME ${CMAKE_CURRENT_LIST_DIR} NAME)
string(REPLACE " " "_" ProjectId ${PROJECT_NAME})

message(${PROJECT_NAME})
project(${PROJECT_NAME})

set(SOURCES
        include_levmar.h
    )

# For any nested file, replace "/" by  "\"
foreach(source ${SOURCES})
    get_filename_component(source_path "${source}" PATH)
    string(REPLACE "/" "\\" source_path_msvc "${source_path}")
    source_group("${source_path_msvc}" FILES "${source}")
endforeach()

add_library(${PROJECT_NAME} STATIC ${SOURCES})
set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE CXX)

include_directories(/home/eudald/Downloads/levmar-2.6)
target_link_libraries(${PROJECT_NAME} /home/eudald/Downloads/levmar-2.6/liblevmar.a)
add_subdirectory(UnitTests)