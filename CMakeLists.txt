cmake_minimum_required(VERSION 3.1)
project(poly)

add_subdirectory(ext)

# include directories
include_directories(${PROJECT_SOURCE_DIR}/src/
                    ${PROJECT_SOURCE_DIR}/include/
                    ${PROJECT_SOURCE_DIR}/ext/
                    ${PROJECT_SOURCE_DIR}/ext/stb/
                    ${PROJECT_SOURCE_DIR}/ext/eigen/
)

set(ext_stb_srcs
    ext/stb/stb_image_write.h
    ext/stb/stb_image.h
)

add_definitions(-DDATA_DIR=R\"\(${PROJECT_SOURCE_DIR}/data\)\")

set(common_srcs
    include/poly/common.h
)

# Compile common library
#add_library(common_lib STATIC
#    ${common_srcs}
#    ${ext_stb_srcs}
#)

#set_target_properties(common_lib PROPERTIES LINKER_LANGUAGE CXX)

#SOURCE_GROUP("ext\\stb" FILES ${ext_stb_srcs})

set(poly_srcs
    ext/stb/stb_image_write.h
    ext/stb/stb_image.h

    include/poly/common.h

    src/main.cpp
)    # a5

add_executable(poly ${poly_srcs})                                               # a5

#target_link_libraries(poly common_lib)                                        # a5
#SOURCE_GROUP("src" FILES ${poly_srcs})                                        # a5
                                                                            # a5

# Create the output data directory
add_custom_command(TARGET poly POST_BUILD COMMAND ${CMAKE_COMMAND} -E make_directory "${PROJECT_SOURCE_DIR}/data/output") # a5
