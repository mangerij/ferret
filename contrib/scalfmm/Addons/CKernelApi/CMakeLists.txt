# check if compiling into source directories
STRING(COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" insource)
if(insource)
    MESSAGE(FATAL_ERROR "${PROJECT_NAME} requires an out of source build. Goto ./Build and tapes cmake ../")
endif(insource)

project(ADDONS_CKERNELAPI_SCALFMM CXX C)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SCALFMM_CXX_FLAGS}")

# Active language
# -----------------------
ENABLE_LANGUAGE(CXX C)
MESSAGE(STATUS " CXX  ${CMAKE_CXX_COMPILER_ID}" )


# if ask to build addon
if(SCALFMM_ADDON_CKERNELAPI)
        # first build lib scalfmmckernelapi
        set(LIBRARY_OUTPUT_PATH ../lib/${CMAKE_BUILD_TYPE})

        # Searching all cpp file
        file( GLOB_RECURSE source_lib_files Src/*.cpp )

        # Adding cpp files to project
        add_library( scalfmmckernelapi STATIC ${source_lib_files} )

        # Add blas library (even if it is set to off)
        target_link_libraries( scalfmmckernelapi scalfmm)

        # Adding the entire project dir as an include dir
        INCLUDE_DIRECTORIES(
         ${CMAKE_BINARY_DIR}/Src
         ${CMAKE_SOURCE_DIR}/Src
         ${SCALFMM_INCLUDES}
    )

        # Install lib
        install( TARGETS scalfmmckernelapi ARCHIVE DESTINATION lib )

        # Install header
        SET(my_include_dirs "Src")

        file( GLOB hpp_in_dir Src/*.hpp Src/*.h)
        INSTALL( FILES ${hpp_in_dir} DESTINATION include/ScalFmm/CKernelApi )

        file( GLOB_RECURSE source_tests_files Tests/*.c )
        INCLUDE_DIRECTORIES( ${CMAKE_BINARY_DIR}/Src )

        # Then build test files
        foreach(exec ${source_tests_files})
                get_filename_component(
                        execname ${exec}
                        NAME_WE
                )

            set(compile_exec "TRUE")

            foreach(fuse_key ${FUSE_LIST})
                file(STRINGS "${exec}" lines_fuse REGEX "@FUSE_${fuse_key}")
                if(lines_fuse)
                        if( NOT SCALFMM_USE_${fuse_key} )
                                MESSAGE( STATUS "This needs ${fuse_key} = ${exec}" )
                                set(compile_exec "FALSE")
                        endif()
                endif()
            endforeach()

                # Dependency are OK
                if( compile_exec )
                        add_executable(	${execname} ${exec} )
                        # link to scalfmm and scalfmmckernelapi
                        target_link_libraries(
                                ${execname}
                                ${scalfmm_lib}
                                scalfmmckernelapi
                                # ${BLAS_LIBRARIES}
                                # ${LAPACK_LIBRARIES}
                            ${SCALFMM_LIBRARIES}
                                )
                endif()
        endforeach(exec)
endif()
