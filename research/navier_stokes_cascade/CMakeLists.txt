add_library(ns_cascade INTERFACE)
target_include_directories(ns_cascade INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_compile_features(ns_cascade INTERFACE cxx_std_11)

add_executable(navier_stokes_cascade
    src/main.cpp
)
target_link_libraries(navier_stokes_cascade PRIVATE ns_cascade)

add_executable(navier_stokes_convergence
    src/convergence.cpp
)
target_link_libraries(navier_stokes_convergence PRIVATE ns_cascade)

add_executable(navier_stokes_fft
    src/fft_main.cpp
)
target_link_libraries(navier_stokes_fft PRIVATE ns_cascade)

add_executable(navier_stokes_search
    src/search.cpp
)
target_link_libraries(navier_stokes_search PRIVATE ns_cascade)

add_executable(navier_stokes_optimize
    src/optimize.cpp
)
target_link_libraries(navier_stokes_optimize PRIVATE ns_cascade)

add_executable(navier_stokes_tangent_check
    src/tangent_check.cpp
)
target_link_libraries(navier_stokes_tangent_check PRIVATE ns_cascade)

add_executable(navier_stokes_adjoint_check
    src/adjoint_check.cpp
)
target_link_libraries(navier_stokes_adjoint_check PRIVATE ns_cascade)

add_executable(navier_stokes_state_optimize
    src/state_optimize.cpp
)
target_link_libraries(navier_stokes_state_optimize PRIVATE ns_cascade)

option(NS_CASCADE_REQUIRE_FFTW_REFERENCE
       "Fail configuration when the independent FFTW oracle is unavailable" OFF)
find_library(NS_CASCADE_FFTW3_LIBRARY NAMES fftw3 fftw3.so.3 libfftw3.so.3)
if(NS_CASCADE_FFTW3_LIBRARY)
    add_executable(navier_stokes_fftw_compare
        src/fftw_compare.cpp
    )
    target_link_libraries(navier_stokes_fftw_compare PRIVATE
        ns_cascade
        ${NS_CASCADE_FFTW3_LIBRARY}
    )

    add_executable(navier_stokes_fftw_tests
        tests/test_fftw_crosscheck.cpp
    )
    target_link_libraries(navier_stokes_fftw_tests PRIVATE
        ns_cascade
        ${NS_CASCADE_FFTW3_LIBRARY}
    )
    if(BUILD_TESTING)
        add_test(NAME navier_stokes_independent_fftw
                 COMMAND navier_stokes_fftw_tests)
    endif()
elseif(NS_CASCADE_REQUIRE_FFTW_REFERENCE)
    message(FATAL_ERROR "FFTW3 is required for the independent research oracle")
else()
    message(STATUS "FFTW3 not found; independent FFT oracle target is skipped")
endif()

if(BUILD_TESTING)
    add_executable(navier_stokes_candidate_evidence_tests
        tests/test_candidate_evidence.cpp
    )
    target_link_libraries(navier_stokes_candidate_evidence_tests PRIVATE ns_cascade)
    add_test(NAME navier_stokes_candidate_evidence
             COMMAND navier_stokes_candidate_evidence_tests)
    find_package(Python3 COMPONENTS Interpreter)
    if(Python3_Interpreter_FOUND)
        add_test(NAME navier_stokes_robust_search
                 COMMAND ${Python3_EXECUTABLE}
                         ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_robust_search.py)
    endif()

    add_executable(navier_stokes_cascade_tests
        tests/test_galerkin.cpp
    )
    target_link_libraries(navier_stokes_cascade_tests PRIVATE ns_cascade)
    add_test(NAME navier_stokes_cascade_invariants
             COMMAND navier_stokes_cascade_tests)

    add_executable(navier_stokes_fft_tests
        tests/test_pseudospectral.cpp
    )
    target_link_libraries(navier_stokes_fft_tests PRIVATE ns_cascade)
    add_test(NAME navier_stokes_fft_oracle
             COMMAND navier_stokes_fft_tests)

    add_executable(navier_stokes_candidate_score_tests
        tests/test_candidate_score.cpp
    )
    target_link_libraries(navier_stokes_candidate_score_tests PRIVATE ns_cascade)
    add_test(NAME navier_stokes_candidate_score
             COMMAND navier_stokes_candidate_score_tests)

    add_executable(navier_stokes_parameter_optimizer_tests
        tests/test_parameter_optimizer.cpp
    )
    target_link_libraries(navier_stokes_parameter_optimizer_tests PRIVATE
        ns_cascade
    )
    add_test(NAME navier_stokes_parameter_optimizer
             COMMAND navier_stokes_parameter_optimizer_tests)

    add_executable(navier_stokes_state_optimizer_tests
        tests/test_state_optimizer.cpp
    )
    target_link_libraries(navier_stokes_state_optimizer_tests PRIVATE
        ns_cascade
    )
    add_test(NAME navier_stokes_state_optimizer
             COMMAND navier_stokes_state_optimizer_tests)
endif()
