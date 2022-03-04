set(TRIANGLE_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/triangle.c)

add_definitions(-DTRILIBRARY)
add_definitions(-w)
add_definitions(-DLINUX)

add_library(triangle STATIC ${TRIANGLE_SOURCES})

