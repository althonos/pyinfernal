find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)

execute_process(
    COMMAND 
        ${Python_EXECUTABLE} -c "import pyhmmer; print(pyhmmer.__path__[0], end='.libs')"
    OUTPUT_VARIABLE 
        PYHMMER_LIBS_DIR
)

add_library(libeasel INTERFACE)
target_include_directories(libeasel INTERFACE ${PYHMMER_LIBS_DIR})
target_link_directories(libeasel INTERFACE ${PYHMMER_LIBS_DIR})
target_link_libraries(libeasel INTERFACE -llibeasel)
