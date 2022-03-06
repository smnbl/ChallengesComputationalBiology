cmake_minimum_required(VERSION 2.6)
project(alignment)

add_executable(globalNW globalNW.cpp)
add_executable(bandedNW bandedNW.cpp)
add_executable(hirschberg hirschberg.cpp)
add_executable(localSW localSW.cpp)
add_executable(overlap overlap.cpp)
add_executable(affine affine.cpp)

add_executable(personal_global personal/global.cpp)
add_executable(personal_dHdV personal/dHdVNW.cpp)

install(TARGETS globalNW RUNTIME DESTINATION bin)
install(TARGETS bandedNW RUNTIME DESTINATION bin)
install(TARGETS hirschberg RUNTIME DESTINATION bin)
install(TARGETS localSW RUNTIME DESTINATION bin)
install(TARGETS overlap RUNTIME DESTINATION bin)
install(TARGETS affine RUNTIME DESTINATION bin)

install(TARGETS personal_global RUNTIME DESTINATION bin)
install(TARGETS personal_dHdV RUNTIME DESTINATION bin)
