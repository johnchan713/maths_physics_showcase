add_library(ns_cascade INTERFACE)
target_include_directories(ns_cascade INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_compile_features(ns_cascade INTERFACE cxx_std_11)

add_executable(navier_stokes_cascade
    src/main.cpp
)
target_link_libraries(navier_stokes_cascade PRIVATE ns_cascade)

if(BUILD_TESTING)
    add_executable(navier_stokes_cascade_tests
        tests/test_galerkin.cpp
    )
    target_link_libraries(navier_stokes_cascade_tests PRIVATE ns_cascade)
    add_test(NAME navier_stokes_cascade_invariants
             COMMAND navier_stokes_cascade_tests)
endif()
