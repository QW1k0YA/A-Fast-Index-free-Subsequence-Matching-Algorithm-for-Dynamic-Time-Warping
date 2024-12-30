cmake_minimum_required(VERSION 3.26)
project(scale_dtw_optimize)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
set(INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(LIB_DIR "${CMAKE_CURRENT_SOURCE_DIR}/lib")

MATH(EXPR stack_size "100*1024*1024")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--stack,${stack_size}")
set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")

set(MY_BIN_DIR "${CMAKE_BINARY_DIR}/bin")
message(STATUS "MY_BIN_DIR is set to: ${MY_BIN_DIR}")

link_directories(${LIB_DIR})


add_executable(scale_dtw_optimize LB_q.cpp
        lower_upper_lemire.cpp
        def.h
        funcs.h
        normal_cdf.cpp
        mvmean.cpp
        ele_multi.cpp
        sip.cpp
        fftw3.h
        ED.cpp
        ED.cpp
        IP.cpp

        LB_t.cpp
        main.cpp
        ele_plus.cpp
        mvsum.cpp
        LB_KK.cpp
        dtw.cpp
        mask_T.cpp
        table_of_q.cpp
        LB_PETITJEAN.cpp
        lb_data.cpp
        dot_mul.cpp
        compute_M.cpp

)

target_link_libraries(scale_dtw_optimize libfftw3-3.lib libfftw3f-3.lib  libfftw3l-3.lib)
