cmake_minimum_required(VERSION 3.2)

if (CMAKE_BUILD_TYPE MATCHES Test)
    set(PROJECT_NAME X_UT_SWIFT)

    message(${PROJECT_NAME})
    project(${PROJECT_NAME})

    file(GLOB_RECURSE SOURCES *.cpp *.h *.ipp)
    add_executable(${PROJECT_NAME} ${SOURCES})
    target_link_libraries(${PROJECT_NAME} PUBLIC SWIFT FFTW3 ${BOOST_STATIC_LIBRARIES})
endif()