cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

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

target_link_libraries(${PROJECT_NAME} PUBLIC glm)

