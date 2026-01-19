cmake_minimum_required(VERSION 3.20)
project(${SKBUILD_PROJECT_NAME} VERSION ${SKBUILD_PROJECT_VERSION} LANGUAGES C)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/src/scripts/cmake")

# --- Find PyHMMER -------------------------------------------------------------

find_package(PyHMMER CONFIG REQUIRED)
set(HMMER_IMPL "${PyHMMER_HMMER_IMPL}")
set(INFERNAL_IMPL "${PyHMMER_HMMER_IMPL}")

# --- Setup Python folders -----------------------------------------------------

cmake_path(APPEND CMAKE_CURRENT_LIST_DIR "src" OUTPUT_VARIABLE SRC_DIR)
set_property(GLOBAL PROPERTY PYTHON_EXTENSIONS_SOURCE_DIR "${SRC_DIR}")
cmake_path(APPEND CMAKE_CURRENT_LIST_DIR "include" OUTPUT_VARIABLE CYTHON_HEADERS_DIR)

# Install shared library objects to `pyinfernal.libs` (like auditwheel repair)
if(PYINFERNAL_INSTALL_LIBS)
    set(PYINFERNAL_INSTALL_LIBS_DIR "${CMAKE_PROJECT_NAME}.libs")
    message(STATUS "Installing shared libraries to ${PYINFERNAL_INSTALL_LIBS_DIR}")
    install(DIRECTORY ${CYTHON_HEADERS_DIR} DESTINATION ${PYINFERNAL_INSTALL_LIBS_DIR})
endif()

# --- Compile source code ------------------------------------------------------

add_subdirectory(src)
