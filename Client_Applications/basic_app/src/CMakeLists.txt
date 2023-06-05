CMAKE_MINIMUM_REQUIRED(VERSION 3.2)

PROJECT (BasicApplication)


## include directories
set(includeDir
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../../../Libraries/Explicit-cpp/include
  ${CMAKE_CURRENT_SOURCE_DIR}/../../../Libraries/myFRIClient/include
  ${CMAKE_CURRENT_SOURCE_DIR}/../../../Libraries/Explicit-cpp/Eigen
)
include_directories(${includeDir})

## make all source files visible in the IDE
set(sourceFiles "")
file(GLOB_RECURSE sourceFiles ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp )

## add executable
add_executable(${PROJECT_NAME} ${sourceFiles} ${headerFiles})


##### LINKING #####

set(EXPLICIT_CPP_LIBRARY "${CMAKE_CURRENT_SOURCE_DIR}/../../../Libraries/Explicit-cpp/lib/libExplicit-cpp.so" CACHE STRING "Explicit-cpp library")
set(FRI_LIBRARY "${CMAKE_CURRENT_SOURCE_DIR}/../../../Libraries/myFRIClient/lib/libFRIClient.a" CACHE STRING "FRI library")

set(linkLibraries
${EXPLICIT_CPP_LIBRARY}
${FRI_LIBRARY}
)


##### SOME PARTS OF THE CODE GO INTO AN ADDITIONAL LIBRARY #####

target_link_libraries(${PROJECT_NAME} ${linkLibraries})


## the output directory
set_target_properties(${PROJECT_NAME}
        PROPERTIES
DEBUG_POSTFIX "${CMAKE_DEBUG_POSTFIX}"
  RUNTIME_OUTPUT_DIRECTORY_DEBUG "${CMAKE_CURRENT_SOURCE_DIR}/../build/Debug"
  RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_CURRENT_SOURCE_DIR}/../build/Release"
)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/ DESTINATION ${CMAKE_INSTALL_PREFIX}/include
  FILES_MATCHING PATTERN "*.h"
)
set(libOutputDir ${CMAKE_INSTALL_PREFIX}/lib)
INSTALL (TARGETS ${PROJECT_NAME}
                LIBRARY DESTINATION ${libOutputDir}
                ARCHIVE DESTINATION ${libOutputDir}
                RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/application
)
