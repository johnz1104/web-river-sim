# Library sources of the lattice Boltzmann solver, shared by the native build
# (fluids/CMakeLists.txt) and the WebAssembly build (lbm/wasm/CMakeLists.txt).
#   LBM_CORE_SOURCES      the river and its solver; everything the website needs
#   LBM_RESEARCH_SOURCES  native-only ML research support (reads model files, so
#                         it is kept out of the website build)

set(LBM_RESEARCH_SOURCES
    ${CMAKE_CURRENT_LIST_DIR}/src/LearnedCollision.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/Coarsen.cpp
)
set(LBM_CORE_SOURCES
    ${CMAKE_CURRENT_LIST_DIR}/src/BankSpline.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/Domain.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/River.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/RiverOptions.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/Solver.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/Stochastic.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/Tracers.cpp
)
set(LBM_INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR}/include)
