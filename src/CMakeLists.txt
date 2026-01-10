file(COPY ../vendor/infernal DESTINATION .)
file(COPY ./infernal DESTINATION .)

add_subdirectory("${CMAKE_CURRENT_BINARY_DIR}/infernal" "${CMAKE_CURRENT_BINARY_DIR}/infernal/build")
add_subdirectory(pyinfernal)
