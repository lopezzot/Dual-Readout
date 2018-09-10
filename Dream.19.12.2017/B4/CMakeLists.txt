# $Id: CMakeLists.txt 68058 2013-03-13 14:47:43Z gcosmo $

#---Adding example N3 subdirectories explicitly 
#   and a custom target to for building all example N3 options ----------

cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

add_subdirectory(B4a)

add_custom_target(B4 DEPENDS exampleB4a)
