# check if compiling into source directories
STRING(COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" insource)
if(insource)
    MESSAGE(FATAL_ERROR "${PROJECT_NAME} requires an out of source build. Goto scalfmm/Build and tapes cmake ../")
endif(insource)

project(UTESTS_SCALFMM CXX)

FIND_PROGRAM( TESTS_GCOV_PATH gcov )
FIND_PROGRAM( TESTS_LCOV_PATH lcov )
FIND_PROGRAM( TESTS_GENHTML_PATH genhtml )
FIND_PROGRAM( TESTS_GCOVR_PATH gcovr PATHS ${CMAKE_SOURCE_DIR}/Modules)

MESSAGE(STATUS  "GG ${TESTS_GCOV_PATH} ${TESTS_LCOV_PATH} ${TESTS_GENHTML_PATH} ${TESTS_GCOVR_PATH}" )
IF( TESTS_GCOV_PATH AND TESTS_LCOV_PATH AND TESTS_GCOVR_PATH )
  #Script for launching lcov and generating an html report
  ADD_CUSTOM_TARGET(coverage)
  ADD_CUSTOM_COMMAND(TARGET coverage 
           COMMAND mkdir -p CoverageReport
          COMMAND ${TESTS_LCOV_PATH} --directory  . --zerocounters
          COMMAND make test
  # Capturing lcov counters and generating report
		COMMAND ${TESTS_LCOV_PATH} --directory . --capture --output-file ${_outputname}.info
		COMMAND ${TESTS_LCOV_PATH} --remove ${_outputname}.info 'tests/*' '/usr/*' --output-file ${_outputname}.info.cleaned
		COMMAND ${TESTS_GENHTML_PATH} -o ${_outputname} ${_outputname}.info.cleaned
		COMMAND ${CMAKE_COMMAND} -E remove ${_outputname}.info ${_outputname}.info.cleaned
           COMMAND echo "Open ${CMAKE_BINARY_DIR}/CoverageReport/index.html to view the report on the coverage of the code."
          WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
          )
MESSAGE(STATUS " GCOV ${COVERAGE_COMMAND}")
MESSAGE(STATUS "CMAKE_CXX_COMPILER_ID ${CMAKE_CXX_COMPILER_ID}")
  IF(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        SET(TESTS_GCOV_CXX_FLAGS "-g -fprofile-arcs -ftest-coverage")
        SET(TESTS_GCOV_LINKER_FLAGS "-g -fprofile-arcs -ftest-coverage -lgcov")
  ENDIF()
 #
 ELSE()
        SET(GCOV_CXX_FLAGS "")
        SET(GCOV_LINKER_FLAGS "")
	   MESSAGE(STATUS  "Could not find gcov in your path - No Coverage option for Tests")
ENDIF()



set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BUILD_TYPE})



ADD_DEFINITIONS( ${SCALFMM_CXX_FLAGS} ${TESTS_GCOV_CXX_FLAGS})

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

		target_link_libraries(
			${execname}
            ${scalfmm_lib}
            ${SCALFMM_LIBRARIES}
			${TESTS_GCOV_LINKER_FLAGS}
		)

        if(lines_mpi)
            add_test(${execname} mpirun -np 4 ${EXECUTABLE_OUTPUT_PATH}/${execname})
        else()
		    add_test(${execname} ${CMAKE_CURRENT_BINARY_DIR}/${execname})
        endif()
	endif()
endforeach(exec)
