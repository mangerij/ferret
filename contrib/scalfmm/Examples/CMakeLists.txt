# check if compiling into source directories
STRING(COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" insource)
if(insource)
    MESSAGE(FATAL_ERROR "${PROJECT_NAME} requires an out of source build. Goto scalfmm/Build and tapes cmake ../")
endif(insource)

project(EXAMPLES_SCALFMM CXX)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BUILD_TYPE})
ADD_DEFINITIONS( ${ScaLFMM_CXX_FLAGS})
# Find all code files
file(	
	GLOB_RECURSE
	source_tests_files
	./*.cpp
	)


# Adding the project sources dir as an include dir
INCLUDE_DIRECTORIES(
     ${CMAKE_BINARY_DIR}/Src 
     ${CMAKE_SOURCE_DIR}/Src 
     ${SCALFMM_INCLUDES}
     
)


# Add execs - 1 cpp = 1 exec
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
		add_executable(
				${execname}
				${exec}
				)
		SET_TARGET_PROPERTIES(${execname} PROPERTIES ENABLE_EXPORTS TRUE)
		target_link_libraries(
			${execname}
			${scalfmm_lib}
			${SCALFMM_LIBRARIES}
			)
	
                   # Install bin
                   install( TARGETS ${execname}   DESTINATION bin )
		endif()

endforeach(exec)
