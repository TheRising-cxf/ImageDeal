# CMake 最低版本号要求
cmake_minimum_required (VERSION 2.8)

# 项目信息
project (ImageDeal)

option (EXPORTS "Use EXPORTS" ON)
# 加入一个配置头文件，用于处理 CMake 对源码的设置
configure_file (
  "${PROJECT_SOURCE_DIR}/config.h.in"
  "${PROJECT_SOURCE_DIR}/include/config.h"
  )
  
include_directories(./include)

aux_source_directory(./src DIR_SRCS)


# 生成链接库
add_library (${PROJECT_NAME} ${DIR_SRCS})