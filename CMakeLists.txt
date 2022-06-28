cmake_minimum_required(VERSION 3.0)


# name of application
set(APP_NAME QHOSYN)

project(${APP_NAME} VERSION 0.8
DESCRIPTION "Quantum Harmonic Oscillator Synthesizer"
LANGUAGES CXX)

# add allolib as a subdirectory to the project
add_subdirectory(external/allolib)
add_subdirectory(external/al_ext)

if (NOT TARGET al_soundfile)
  message(FATAL_ERROR "al_soundfile not built. Likely libsndfile not found")
endif()

set(includes
src/QHOSYN.hpp
src/shadercode.hpp
src/wavefunction.hpp
src/utility.hpp
external/fparser4.5.2/fparser.cc
external/fparser4.5.2/fpoptimizer.cc
)

set(src
src/QHOSYN.cpp
)

set(LIB_NAME QHOSYNLib)
add_library(${LIB_NAME}
  STATIC 
    ${includes}
    ${srcs}
)



target_include_directories(${LIB_NAME} PRIVATE 
  ${CMAKE_SOURCE_DIR}/external
  ${CMAKE_SOURCE_DIR}/src
)

set_target_properties(${LIB_NAME} PROPERTIES
  CXX_STANDARD 17
  CXX_STANDARD_REQUIRED ON
)

# Move libsndfile-1.dll to bin/Resources
if(WIN32)
  get_filename_component(sndfile_dir ${SNDFILE_LIBRARY} DIRECTORY)
  set(${sndfile_dir}/../bin/libsndfile-1.dll)
  message("WINDOWS: Attempting to copy sndfile dll to bin/")
  file(COPY ${sndfile_dir}/../bin/libsndfile-1.dll
        DESTINATION ${CMAKE_SOURCE_DIR}/bin)
endif()

# OSX library needed for native file dialog.
if(APPLE)
  find_library(AppKit Appkit)
endif()

# Linux library needed for native file dialog.
if(UNIX AND NOT APPLE)
  set(LINUX TRUE)
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(GTK3 REQUIRED gtk+-3.0)
  include_directories(${GTK3_INCLUDE_DIRS})
  link_directories(${GTK3_LIBRARY_DIRS})
endif()

if (EXISTS ${CMAKE_CURRENT_LIST_DIR}/al_ext)
  message("Buiding extensions in al_ext")
  add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/al_ext)
  get_target_property(AL_EXT_LIBRARIES al_ext AL_EXT_LIBRARIES)
  target_link_libraries(${APP_NAME} PRIVATE ${AL_EXT_LIBRARIES})
endif()

# bunch of stuff for packaging MacOS Apps
if(APPLE)
  find_library(AppKit Appkit)
  target_link_libraries(${LIB_NAME} PUBLIC ${Appkit})

  set(RESOURCE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/bin/QHOSYN.app/Contents/Resources")
  set(CMAKE_MACOSX_RPATH 1)

  set(MACOSX_BUNDLE_ICON_FILE "QHOSYN.icns")
  set(ICON ${CMAKE_CURRENT_SOURCE_DIR}/deployment/icons/QHOSYN.icns)
  set_source_files_properties(${ICON} PROPERTIES
       MACOSX_PACKAGE_LOCATION "Resources")

  set(CMAKE_INSTALL_RPATH "@executable_path/../Resources/")

  set(RESOURCE_LIST
    "libsndfile"
  )

  foreach(RESOURCE_DIR_NAME ${RESOURCE_LIST})
    SET(RESOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/deployment/externalResources/${RESOURCE_DIR_NAME}")
    if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
      file(GLOB CHILDREN RELATIVE "${RESOURCE_DIR}" ${RESOURCE_DIR}/*)
      foreach(CHILD ${CHILDREN})
        IF(IS_DIRECTORY "${RESOURCE_DIR}/${CHILD}")
          # be wary of .DS_Store files. glob them out by using [A-Za-z]*.*
          FILE(GLOB RESOURCE_DIR_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "${RESOURCE_DIR}/${CHILD}/[A-Za-z]*.*")
          SET_SOURCE_FILES_PROPERTIES(
            ${RESOURCE_DIR_FILES}
            PROPERTIES MACOSX_PACKAGE_LOCATION Resources/${RESOURCE_DIR_NAME}/${CHILD}
          )
          SET(RESOURCE_FILES ${RESOURCE_FILES} "deployment/externalResources/${RESOURCE_DIR_FILES}")
        ELSEIF(NOT ${CHILD} MATCHES "\\.DS_Store|.*\\.am")
          SET_SOURCE_FILES_PROPERTIES(
            "${RESOURCE_DIR}/${CHILD}"
            PROPERTIES MACOSX_PACKAGE_LOCATION Resources/${RESOURCE_DIR_NAME}
          )
          SET(RESOURCE_FILES ${RESOURCE_FILES} "deployment/externalResources/${RESOURCE_DIR_NAME}/${CHILD}")
        ENDIF()
    endforeach()
  else(${CMAKE_SYSTEM_NAME} MATCHES "Darwin") #WORK ON
    INSTALL( DIRECTORY ${RESOURCE_DIR} DESTINATION "${DATA_DIR}" PATTERN "*.am" EXCLUDE)
  endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  endforeach()

  add_executable(${APP_NAME} MACOSX_BUNDLE
    src/QHOSYN.cpp
    ${RESOURCE_FILES}
    ${ICON}
  )
  
endif(APPLE)

if(WIN32)
  set(ICON ${CMAKE_CURRENT_SOURCE_DIR}/deployment/icons/QHOSYN-windows.rc)
  # path to main source file
  add_executable(${APP_NAME} 
    src/QHOSYN.cpp
    ${ICON}
  )
endif(WIN32)

if(UNIX AND NOT APPLE)
  # path to main source file
  add_executable(${APP_NAME} src/QHOSYN.cpp)
endif()

target_include_directories(${APP_NAME} PRIVATE
${CMAKE_SOURCE_DIR}/external
${CMAKE_SOURCE_DIR}/external/nativefiledialog/src/include
)

find_library(NFD_LIBRARY nfd HINTS 
  ${CMAKE_SOURCE_DIR}/external/nativefiledialog/build/lib/Release/x64
)

# link allolib to project
get_target_property(AL_EXT_LIBRARIES al_ext AL_EXT_LIBRARIES)

if(NOT WIN32)
    find_package(GSL REQUIRED)
    target_link_libraries(${LIB_NAME} PRIVATE al GSL::gsl GSL::gslcblas fftw3)
endif()

if(APPLE)
    find_package(Boost REQUIRED)
    target_link_libraries(${LIB_NAME} PRIVATE ${Boost_LIBRARIES})
endif()

if(WIN32)
    if (NOT DEFINED ENV{VCPKG_ROOT})
        message(FATAL_ERROR "VCPKG_ROOT environment variable not defined")
    endif()
    set(VCPKG_ROOT $ENV{VCPKG_ROOT})
    if(CMAKE_SIZEOF_VOID_P MATCHES 8)
        # Set 64-bit stuff
        set(FFTW_ROOT "${VCPKG_ROOT}/packages/fftw3_x64-windows")
    else()
        # Set 32-bit stuff
        set(FFTW_ROOT "${VCPKG_ROOT}/packages/fftw3_x86-windows")
    endif()
    set(FFTW3_DIR "${FFTW_ROOT}/share/fftw3")
    include_directories("${FFTW_ROOT}/include")
    # find_library(FFTW3 "${FFTW_ROOT}/lib")
    find_package(FFTW3 REQUIRED)
    find_package(GSL REQUIRED)
    target_link_libraries(${LIB_NAME} PRIVATE al GSL::gsl GSL::gslcblas FFTW3::fftw3)
endif()

target_link_libraries(${LIB_NAME} PUBLIC 
  al 
  "${AL_EXT_LIBRARIES}"
  ${NFD_LIBRARY}
) 

if(APPLE)
  target_link_libraries(${LIB_NAME} PUBLIC ${Appkit})
endif()

if(LINUX)
  target_link_libraries(${LIB_NAME} PUBLIC ${GTK3_LIBRARIES})
endif()

# Link windows libraries.
if(WIN32)
  target_link_libraries(${LIB_NAME} PUBLIC comctl32)
endif()


if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    # On macOS if jack is present, RtMidi will link to jack but not provide the
    # /usr/local/lib location as link directory. This statement only serves the
    # purpose of working around that bug. Perhaps can be removed on future
    # updates of RtMidi
    link_directories("/usr/local/lib")
    link_directories("/Library/Frameworks/Jackmp.framework/Versions/A/")
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

target_link_libraries(${APP_NAME} PUBLIC ${LIB_NAME} ${LIBRARIES})

set_target_properties(${LIB_NAME} PROPERTIES
  CXX_STANDARD 17
  CXX_STANDARD_REQUIRED ON
)

# binaries are put into the ./bin directory by default
set_target_properties(${APP_NAME} PROPERTIES
  BUILD_WITH_INSTALL_RPATH 1
  INSTALL_NAME_DIR "@rpath"
  CXX_STANDARD 17
  CXX_STANDARD_REQUIRED ON
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/bin
  RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_CURRENT_LIST_DIR}/bin
  RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_CURRENT_LIST_DIR}/bin
)

if(APPLE)
  set_target_properties(${APP_NAME} PROPERTIES
  MACOSX_BUNDLE_BUNDLE_NAME "${APP_NAME}"
  MACOSX_BUNDLE_GUI_IDENTIFIER "${APP_NAME}"
  MACOSX_BUNDLE_INFO_PLIST ${CMAKE_CURRENT_SOURCE_DIR}/deployment/MacOS/MacOSBundleInfo.plist.in
  )
  add_custom_command (TARGET ${APP_NAME}
    POST_BUILD 
    COMMAND "${CMAKE_INSTALL_NAME_TOOL}"
      "-change" "/usr/local/opt/libsndfile/lib/libsndfile.1.dylib" 
      "@rpath/libsndfile/libsndfile.1.dylib"
      "$<TARGET_FILE:${APP_NAME}>" 

    COMMAND "${CMAKE_INSTALL_NAME_TOOL}"
      "-change" "/usr/local/opt/libvorbis/lib/libvorbis.0.dylib" 
      "@rpath/libsndfile/libvorbis.0.dylib"
      "${RESOURCE_PATH}/libsndfile/libsndfile.1.dylib"
    COMMAND "${CMAKE_INSTALL_NAME_TOOL}"
      "-change" "/usr/local/opt/flac/lib/libFLAC.8.dylib" 
      "@rpath/libsndfile/libFLAC.8.dylib"
      "${RESOURCE_PATH}/libsndfile/libsndfile.1.dylib"
    COMMAND "${CMAKE_INSTALL_NAME_TOOL}"
      "-change" "/usr/local/opt/libogg/lib/libogg.0.dylib" 
      "@rpath/libsndfile/libogg.0.dylib"
      "${RESOURCE_PATH}/libsndfile/libsndfile.1.dylib"
    COMMAND "${CMAKE_INSTALL_NAME_TOOL}"
      "-change" "/usr/local/opt/libvorbis/lib/libvorbisenc.2.dylib" 
      "@rpath/libsndfile/libvorbisenc.2.dylib"
      "${RESOURCE_PATH}/libsndfile/libsndfile.1.dylib" 

    COMMAND "${CMAKE_INSTALL_NAME_TOOL}"
    "-change" "/usr/local/opt/libogg/lib/libogg.0.dylib" 
    "@rpath/libsndfile/libogg.0.dylib"
    "${RESOURCE_PATH}/libsndfile/libFLAC.8.dylib"

    COMMAND "${CMAKE_INSTALL_NAME_TOOL}"
      "-change" "/usr/local/opt/libogg/lib/libogg.0.dylib" 
      "@rpath/libsndfile/libogg.0.dylib"
      "${RESOURCE_PATH}/libsndfile/libvorbis.0.dylib"

    COMMAND "${CMAKE_INSTALL_NAME_TOOL}"
      "-change" "/usr/local/opt/libogg/lib/libogg.0.dylib" 
      "@rpath/libsndfile/libogg.0.dylib"
      "${RESOURCE_PATH}/libsndfile/libvorbisenc.2.dylib" 
    COMMAND "${CMAKE_INSTALL_NAME_TOOL}"
      "-change" "/usr/local/Cellar/libvorbis/1.3.6/lib/libvorbis.0.dylib" 
      "@rpath/libsndfile/libvorbis.0.dylib"
      "${RESOURCE_PATH}/libsndfile/libvorbisenc.2.dylib" 
    VERBATIM
  )

endif(APPLE)