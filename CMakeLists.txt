CMAKE_MINIMUM_REQUIRED(VERSION 2.8.7 FATAL_ERROR)
PROJECT(TRICAL)

INCLUDE_DIRECTORIES(./include ./src)

set(CMAKE_C_FLAGS "-O3 -Weverything -Werror -Wno-documentation -Wno-padded -Wno-unknown-pragmas -fPIC")

ADD_LIBRARY(TRICAL SHARED src/TRICAL.c src/filter.c)
ADD_LIBRARY(TRICALstatic STATIC src/TRICAL.c src/filter.c)

ENABLE_TESTING()
ADD_SUBDIRECTORY(test EXCLUDE_FROM_ALL)
