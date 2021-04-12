# CMakeList.txt: проект CMake для CForsythia; включите исходный код и определения,
# укажите здесь логику для конкретного проекта.
#
cmake_minimum_required (VERSION 3.8)

project ("CForsythia")

add_definitions("-DEBUG")
# Добавьте источник в исполняемый файл этого проекта.
add_executable (CForsythia "utils.cpp" "utils.h" "CForsythia.cpp" "CForsythia.h" )

# TODO: Добавьте тесты и целевые объекты, если это необходимо.
