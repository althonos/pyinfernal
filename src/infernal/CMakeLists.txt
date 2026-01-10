# --- Handle config.h generation *a la* ./configure ----------------------------

include(CheckIncludeFile)
include(CheckFunctionExists)

CHECK_INCLUDE_FILE("alloca.h" HAVE_ALLOCA_H)

set(INFERNAL_DATE "Sep 2023")
set(INFERNAL_COPYRIGHT "Copyright (C) 2023 Howard Hughes Medical Institute.")
set(INFERNAL_LICENSE  "Freely distributed under the BSD open source license.")
set(INFERNAL_VERSION "0.49")
set(INFERNAL_URL "http://bioeasel.org/")

configure_file(config.h.cmake config.h)
file(READ ${CMAKE_CURRENT_BINARY_DIR}/config.h CONFIG_H)

message(DEBUG "-- Show configuration ")
message(DEBUG "HAVE_ALLOCA_H: ${HAVE_ALLOCA_H}")
message(DEBUG "-- Dumping config.h --")
message(DEBUG "${CONFIG_H}")
message(DEBUG "-- End config.h --")

set(INFERNAL_PATCHED_SOURCES
    ${INFERNAL_PATCHED_SOURCES}
    ${CMAKE_CURRENT_BINARY_DIR}/config.h
)

# --- Manually declare required sources only -----------------------------------

set(INFERNAL_SOURCE_DIR
    "${PROJECT_SOURCE_DIR}/vendor/infernal/src"
)

set(INFERNAL_SOURCES
    infernal.h
    cm.c
    cm_alidisplay.c
    cm_alndata.c
    cm_dpalign.c
    cm_dpalign_trunc.c
    cm_dpsearch.c
    cm_dpsearch_trunc.c
    cm_dpsmall.c
    cm_file.c
    cm_modelconfig.c
    cm_modelmaker.c
    cm_mx.c
    cm_parsetree.c
    cm_pipeline.c
    cm_qdband.c
    cm_submodel.c
    cm_tophits.c
    cm_trunc.c
    cm_p7_band.c
    cm_p7_domaindef.c
    cm_p7_modelconfig_trunc.c
    cm_p7_modelmaker.c
    cp9.c
    cp9_dp.c
    cp9_modelmaker.c
    cp9_mx.c
    cp9_trace.c
    alphabet.c
    display.c
    errors.c
    eweight.c
    hmmband.c
    logsum.c
    mpisupport.c
    prior.c
    rnamat.c
    stats.c
    truncyk.c
)

# --- Detect platform-accelerated implementation -------------------------------

if(INFERNAL_IMPL STREQUAL "SSE")
    set(INFERNAL_SOURCES
        ${INFERNAL_SOURCES}
        impl_sse/impl_sse.h
        impl_sse/cm_optimized.c
        impl_sse/sse_cm_dpsearch.c
        impl_sse/sse_cm_dpsmall.c
        impl_sse/sse_cmcons_hitmx.c
        impl_sse/sse_cmcons_mscyk.c
        impl_sse/sse_cmsearch.c
        impl_sse/sse_util.c
    )
else()
    message(FATAL_ERROR "No Infernal implementation defined")
endif()

# --- Copy/patch source files if needed ----------------------------------------

foreach(_file IN ITEMS ${INFERNAL_SOURCES})
    if(EXISTS ${PROJECT_SOURCE_DIR}/patches/${_file}.patch)
        if(NOT EXISTS ${INFERNAL_SOURCE_DIR}/${_file})
            add_custom_command(
                OUTPUT
                    ${_file}
                COMMENT
                    "Creating ${_file}"
                COMMAND
                    cmake -E touch ${CMAKE_CURRENT_BINARY_DIR}/empty.txt
                COMMAND
                    ${Python_EXECUTABLE} ${PROJECT_SOURCE_DIR}/src/scripts/apply_patch.py
                        --input ${CMAKE_CURRENT_BINARY_DIR}/empty.txt
                        --patch ${PROJECT_SOURCE_DIR}/patches/${_file}.patch
                        --output ${CMAKE_CURRENT_BINARY_DIR}/${_file}
                DEPENDS
                    ${PROJECT_SOURCE_DIR}/patches/${_file}.patch
            )
        else()
            add_custom_command(
                OUTPUT
                    ${_file}
                COMMENT
                    "Patching ${_file}"
                COMMAND
                    ${Python_EXECUTABLE} ${PROJECT_SOURCE_DIR}/src/scripts/apply_patch.py
                        --input ${INFERNAL_SOURCE_DIR}/${_file}
                        --patch ${PROJECT_SOURCE_DIR}/patches/${_file}.patch
                        --output ${CMAKE_CURRENT_BINARY_DIR}/${_file}
                DEPENDS
                    ${INFERNAL_SOURCE_DIR}/${_file}
                    ${PROJECT_SOURCE_DIR}/patches/${_file}.patch
            )
        endif()
    else()
        add_custom_command(
            OUTPUT
                ${_file}
            COMMAND
                cmake -E copy ${INFERNAL_SOURCE_DIR}/${_file} ${CMAKE_CURRENT_BINARY_DIR}/${_file}
            DEPENDS
                ${INFERNAL_SOURCE_DIR}/${_file}
        )
    endif()
    set(INFERNAL_PATCHED_SOURCES ${INFERNAL_PATCHED_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/${_file})
endforeach()

# --- Build libinfernal with copied/patched sources ----------------------------

add_library(libinfernal ${INFERNAL_PATCHED_SOURCES})
target_include_directories(libinfernal PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(libinfernal PUBLIC ${PyHMMER_INCLUDE_DIRS})
target_link_libraries(libinfernal PUBLIC ${PyHMMER_LIBRARIES})

# --- Install shared libraries & headers if required ---------------------------

if(DEFINED PYINFERNAL_INSTALL_LIBS_DIR)
    install(TARGETS libinfernal DESTINATION "${SKBUILD_PLATLIB_DIR}/${PYINFERNAL_INSTALL_LIBS_DIR}")
    cmake_path(APPEND SKBUILD_PLATLIB_DIR "${PYINFERNAL_INSTALL_LIBS_DIR}" "include" "libinfernal" OUTPUT_VARIABLE INCLUDE_DEST_FOLDER)

    # install(FILES "${CMAKE_CURRENT_BINARY_DIR}/config.h" "${INFERNAL_SOURCE_DIR}/infernal.h" DESTINATION "${INCLUDE_DEST_FOLDER}")
    # install(FILES "${IMPL_SLUG}/${IMPL_SLUG}.h" DESTINATION "${INCLUDE_DEST_FOLDER}/${IMPL_SLUG}")

    foreach(_file IN ITEMS ${INFERNAL_PATCHED_SOURCES})
         cmake_path(GET _file EXTENSION LAST_ONLY _ext)
         if(_ext STREQUAL ".h")
            cmake_path(RELATIVE_PATH _file BASE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" OUTPUT_VARIABLE _local)
            cmake_path(GET _local PARENT_PATH _dest)
            install(FILES "${_file}" DESTINATION "${INCLUDE_DEST_FOLDER}/${_dest}")
         endif()
    endforeach()

    set_target_properties(libinfernal PROPERTIES INSTALL_RPATH "\$ORIGIN/.")
endif()