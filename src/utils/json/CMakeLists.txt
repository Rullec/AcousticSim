cmake_minimum_required(VERSION 3.5)
project(DeepMimicCore)

add_definitions("-O3 -Wall -fPIC")
include_directories(..)
aux_source_directory(. DIR_SRCS)
add_library(JSON ${DIR_SRCS})