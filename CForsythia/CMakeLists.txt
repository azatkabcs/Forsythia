# CMakeList.txt: проект CMake для CForsythia; включите исходный код и определения,
# укажите здесь логику для конкретного проекта.
#
cmake_minimum_required (VERSION 3.8)
project ("CForsythia")

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMake/Modules)
find_package(GMP REQUIRED)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/out")

INCLUDE_DIRECTORIES(
    "${GMP_INCLUDES}"
)
# Добавьте источник в исполняемый файл этого проекта.
add_executable (CForsythia "utils.cpp" "utils.h" "gfp2.cpp" "gfp2.h" "parameters.h" "parameters.cpp" "montgomery.h" "montgomery.cpp" "CForsythia.cpp" "CForsythia.h" )
# Линковка необходимых библиотек
target_link_libraries(CForsythia ${GMP_LIBRARIES})

# TODO: Добавьте тесты и целевые объекты, если это необходимо.
