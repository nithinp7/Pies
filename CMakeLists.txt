cmake_minimum_required(VERSION 3.24 FATAL_ERROR)

project(
    Pies
    VERSION 0.1.0
    LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include_directories(Include/Pies Src)

function(glob_files out_var_name regexes)
    set(files "")
    foreach(arg ${ARGV})
        list(APPEND regexes_only "${arg}")
    endforeach()
    list(POP_FRONT regexes_only)
    if (GLOB_USE_CONFIGURE_DEPENDS)
        file(GLOB_RECURSE files CONFIGURE_DEPENDS ${regexes_only})
    else()
        file(GLOB files ${regexes_only})
    endif()
    set(${ARGV0} "${files}" PARENT_SCOPE)
endfunction()

glob_files(SRC_FILES_LIST Src/*.cpp)
add_library(Pies ${SRC_FILES_LIST})

if (NOT TARGET glm)
  add_subdirectory(Extern/glm)
endif()

add_subdirectory(Extern/tetgen)
add_subdirectory(Extern/parallel-hashmap)

target_compile_definitions(
  ${PROJECT_NAME}
  PUBLIC
      GLM_FORCE_XYZW_ONLY # Disable .rgba and .stpq to make it easier to view values from debugger
      GLM_FORCE_EXPLICIT_CTOR # Disallow implicit conversions between dvec3 <-> dvec4, dvec3 <-> fvec3, etc
      GLM_FORCE_SIZE_T_LENGTH # Make vec.length() and vec[idx] use size_t instead of int
)

target_include_directories (Pies
    PUBLIC
      Include
      Extern/glm
      Extern/tetgen
      Extern/parallel-hashmap
      Extern/eigen)

target_link_libraries(${PROJECT_NAME} PUBLIC glm)
target_link_libraries(${PROJECT_NAME} PUBLIC tetgen)
target_link_libraries(${PROJECT_NAME} PUBLIC phmap)


