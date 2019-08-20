cmake_minimum_required(VERSION 2.6)
project(suffixtree)

add_executable(suffixtree main.cpp suffixtree.cpp)

install(TARGETS suffixtree RUNTIME DESTINATION bin)
