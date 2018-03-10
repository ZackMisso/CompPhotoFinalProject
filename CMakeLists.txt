cmake_minimum_required(VERSION 3.1)

project(poly)

if(APPLE)
    set(CMAKE_MACOSX_RPATH ON)
endif()

include(CheckCXXCompilerFlag)
include(CheckCXXSourceRuns)

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  # Quench annoying deprecation warnings when compiling GLFW on OSX
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-deprecated-declarations")
endif()

if (CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|Intel)$")
  CHECK_CXX_COMPILER_FLAG("-std=c++14" HAS_CPP14_FLAG)
  CHECK_CXX_COMPILER_FLAG("-std=c++11" HAS_CPP11_FLAG)

  if (HAS_CPP14_FLAG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
  elseif (HAS_CPP11_FLAG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  else()
    message(FATAL_ERROR "Unsupported compiler -- pybind11 requires C++11 support!")
  endif()
endif()

if (APPLE)
  find_library(cocoa_library Cocoa)
  find_library(opengl_library OpenGL)
  find_library(corevideo_library CoreVideo)
  find_library(iokit_library IOKit)
  list(APPEND NANOGUI_EXTRA_LIBS ${cocoa_library} ${opengl_library} ${corevideo_library} ${iokit_library})
  list(APPEND LIBNANOGUI_EXTRA_SOURCE src/darwin.mm)
 endif()

# Disable building extras we won't need (pure C++ project)
set(NANOGUI_BUILD_EXAMPLE OFF CACHE BOOL " " FORCE)
set(NANOGUI_BUILD_PYTHON  OFF CACHE BOOL " " FORCE)
set(NANOGUI_INSTALL       OFF CACHE BOOL " " FORCE)

add_subdirectory(ext/nanogui)
add_subdirectory(ext)

set(GLFW_INCLUDE_DIR
  ${CMAKE_CURRENT_SOURCE_DIR}/ext/nanogui/ext/glfw/include)
set(GLEW_INCLUDE_DIR
  ${CMAKE_CURRENT_SOURCE_DIR}/ext/nanogui/ext/glew/include)
set(NANOVG_INCLUDE_DIR
  ${CMAKE_CURRENT_SOURCE_DIR}/ext/nanogui/ext/nanovg/src)
set(NANOGUI_INCLUDE_DIR
  ${CMAKE_CURRENT_SOURCE_DIR}/ext/nanogui/include)
set(EIGEN_INCLUDE_DIR
  ${CMAKE_CURRENT_SOURCE_DIR}/ext/nanogui/ext/eigen)

find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

# include directories
include_directories(# Include directory
                    #${CMAKE_CURRENT_SOURCE_DIR}/include/
                    # Random number generator
                    #${PCG32_INCLUDE_DIR}/
                    # GLFW library for OpenGL context creation
                    ${GLFW_INCLUDE_DIR}
                    # GLEW library for accessing OpenGL functions
                    ${GLEW_INCLUDE_DIR}
                    # NanoVG drawing library
                    ${NANOVG_INCLUDE_DIR}
                    # NanoGUI user interface library
                    ${NANOGUI_INCLUDE_DIR}
                    ${NANOGUI_EXTRA_INCS}
                    ${PROJECT_SOURCE_DIR}/src/
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
    include/poly/image.h
    include/poly/utils.h
    include/poly/gui.h
    include/poly/indexmap.h
    include/poly/subscreen.h
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
    include/poly/image.h
    include/poly/utils.h
    include/poly/gui.h
    include/poly/indexmap.h
    include/poly/subscreen.h

    src/main.cpp
    src/image.cpp
    src/gui.cpp
    src/indexmap.cpp
    src/subscreen.cpp
)    # a5

add_executable(poly ${poly_srcs})                                               # a5

target_link_libraries(poly nanogui ${NANOGUI_EXTRA_LIBS} ${PYTHON_LIBRARIES})                                    # a5
#SOURCE_GROUP("src" FILES ${poly_srcs})

set(CompilerFlags
      CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
      CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO CMAKE_C_FLAGS
      CMAKE_C_FLAGS_DEBUG CMAKE_C_FLAGS_RELEASE CMAKE_C_FLAGS_MINSIZEREL
      CMAKE_C_FLAGS_RELWITHDEBINFO COMPILE_DEFINITIONS U_CMAKE_BUILD_TYPE
      CMAKE_MACOSX_RPATH
       PCG32_INCLUDE_DIR
      GLFW_INCLUDE_DIR GLEW_INCLUDE_DIR
      NANOVG_INCLUDE_DIR NANOGUI_INCLUDE_DIR EIGEN_INCLUDE_DIR
      NANOGUI_EXTRA_LIBS
)                                      # a5
                                                                            # a5

# Create the output data directory
add_custom_command(TARGET poly POST_BUILD COMMAND ${CMAKE_COMMAND} -E make_directory "${PROJECT_SOURCE_DIR}/data/output") # a5
