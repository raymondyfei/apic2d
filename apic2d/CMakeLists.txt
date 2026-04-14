# APIC2D Executable

append_files (Headers "h" .)
append_files (Sources "cpp" .)
append_files (Sources "c" .)

option(SINGLE_THREADED "simulate with single thread only" OFF)

if(CMAKE_CXX_COMPILER MATCHES "/em\\+\\+")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -s LEGACY_GL_EMULATION=1 -s GL_FFP_ONLY=1")
  set(CMAKE_EXECUTABLE_SUFFIX ".html")
  set(SINGLE_THREADED ON CACHE BOOL "" FORCE)
else()
  # Locate OpenGL
  find_package (OpenGL REQUIRED)
  if (OPENGL_FOUND)
    include_directories (${OPENGL_INCLUDE_DIR})
    set (APIC2D_LIBRARIES ${APIC2D_LIBRARIES} ${OPENGL_LIBRARIES})
  else (OPENGL_FOUND)
    message (SEND_ERROR "Unable to locate OpenGL")
  endif (OPENGL_FOUND)

  if(WIN32)
    set (APIC2D_LIBRARIES_RELEASE ${APIC2D_LIBRARIES} "${THIRDPARTY_WIN64_COMMON_LIBRARIES_PATH}/freeglut.lib")
    set (APIC2D_LIBRARIES_DEBUG ${APIC2D_LIBRARIES} "${THIRDPARTY_WIN64_COMMON_LIBRARIES_PATH}/freeglutd.lib")
  else()
    # Locate GLUT
    find_package (FreeGLUT REQUIRED)
    if (FREEGLUT_FOUND)
      include_directories (${FREEGLUT_INCLUDE_DIR})
      set (APIC2D_LIBRARIES ${APIC2D_LIBRARIES} ${FREEGLUT_LIBRARIES})
    else ()
      message (SEND_ERROR "Unable to locate FreeGLUT")
    endif ()

    if ((NOT SINGLE_THREADED) AND CMAKE_CXX_COMPILER MATCHES "/Xcode")
      set(USE_TBB ON CACHE BOOL "" FORCE)
    endif()
  endif ()

  if(NOT SINGLE_THREADED)
    # Use TBB by default.
    option(USE_TBB "parallelize with TBB" ON)
    if (USE_TBB)
      add_definitions(-DUSE_TBB)
      if(WIN32)
        set (APIC2D_LIBRARIES_RELEASE ${APIC2D_LIBRARIES_RELEASE} "${THIRDPARTY_WIN64_COMMON_LIBRARIES_PATH}/tbb.lib")
        set (APIC2D_LIBRARIES_DEBUG ${APIC2D_LIBRARIES_DEBUG} "${THIRDPARTY_WIN64_COMMON_LIBRARIES_PATH}/tbb_debug.lib")
      else()
        find_package (TBB REQUIRED)
        if (TBB_FOUND)
          include_directories (${TBB_INCLUDE_DIRS})
          set (APIC2D_LIBRARIES ${APIC2D_LIBRARIES} ${TBB_LIBRARIES})
        else (TBB_FOUND)
          message (SEND_ERROR "Unable to locate TBB")
        endif (TBB_FOUND)  
      endif()
    endif()
  endif()

  if(CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "arm64")
    # AArch64 always has NEON; -mfpu is AArch32-only and rejected by Apple clang.
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mcpu=native")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
  endif()  
endif()

if(SINGLE_THREADED)
  add_definitions(-DSINGLE_THREADED)
endif()

add_executable (APIC2D ${Headers} ${Templates} ${Sources})

if (USE_TBB)
  set_property(TARGET APIC2D PROPERTY CXX_STANDARD 14)
else ()
  set_property(TARGET APIC2D PROPERTY CXX_STANDARD 17)
endif ()

if (WIN32)
  add_dependencies (APIC2D DllCopy)
  target_link_libraries (APIC2D optimized ${APIC2D_LIBRARIES_RELEASE} )
  target_link_libraries (APIC2D debug ${APIC2D_LIBRARIES_DEBUG} )
else()
  target_link_libraries (APIC2D ${APIC2D_LIBRARIES} )
endif ()

