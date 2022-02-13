add_library(mps SHARED 
    filesio.cc
)

target_include_directories(mps PUBLIC ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/src)