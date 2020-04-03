# check if compiling into source directories
STRING(COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" insource)
if(insource)
    MESSAGE(FATAL_ERROR "${PROJECT_NAME} requires an out of source build. Goto scalfmm/Build and tapes cmake ../")
endif(insource)

# Configuration du projet
project(Lib_scalfmm)
set(LIBRARY_OUTPUT_PATH ../lib/${CMAKE_BUILD_TYPE})

ADD_DEFINITIONS(${ScaLFMM_CXX_FLAGS})

# Searching all cpp file
file(
	GLOB_RECURSE
	source_lib_files
	./*.cpp
)

set(source_lib_files_enabled)
foreach(exec ${source_lib_files}) 
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
        list(APPEND source_lib_files_enabled ${exec})
	endif()
endforeach(exec)


# Add CUDA files once they are compiled from cu to .o
if(SCALFMM_USE_CUDA)
    # Find all the CU files in my project
    file(GLOB_RECURSE source_cu_files ./*.cu)
 
    # Iterate and add builind command for each file
    set( SCALFMM_CUDA_SOURCES "" )
    FOREACH (_file ${source_cu_files})
 
        GET_FILENAME_COMPONENT (_filewe ${_file} NAME_WE)
        SET (_filehpp_output ${CMAKE_CURRENT_BINARY_DIR}/${_filewe}.o)
 
        ADD_CUSTOM_COMMAND(OUTPUT  ${_filehpp_output}
                  DEPENDS ${_file}
                  COMMAND echo ARGS Compiling ${_filewe}
                  COMMAND nvcc ARGS -I${CUDA_NEEDED_INCLUDE_DIRS} -c ${_file} -o ${_filehpp_output} ${CUSTOM_CUDA_FLAGS})
 
        LIST (APPEND SCALFMM_CUDA_SOURCES ${_filehpp_output})
    ENDFOREACH ()
 
    MESSAGE( STATUS "SCALFMM_CUDA_SOURCES          = ${SCALFMM_CUDA_SOURCES}" )
endif()

# Adding cpp files to project
add_library(
	scalfmm
	STATIC
	${source_lib_files_enabled}
	${SCALFMM_CUDA_SOURCES}
)

# Add blas library (even if it is set to off)
target_link_libraries(
	scalfmm
	${scalfmm_lib}
	${SCALFMM_LIBRARIES}
)

# Adding the entire project dir as an include dir
INCLUDE_DIRECTORIES(
     ${CMAKE_SOURCE_DIR}/Src 
     ${CMAKE_BINARY_DIR}/Src 
     ${SCALFMM_INCLUDES}
     
)


# Install lib
install( TARGETS scalfmm ARCHIVE DESTINATION lib )

# Install header
SET(my_include_dirs "."
  "Arranger"        "Core"           "Utils"
  "Kernels/Chebyshev"       "Components"       "Extensions"  
  "Containers"       "Files"       "Kernels/Spherical"
  "Kernels/P2P" "Kernels/Taylor" "Kernels/ROtation"
  "Kernels/Uniform" "Kernels/Interpolation")

if(SCALFMM_ONLY_DEVEL)
    SET(my_include_dirs ${my_include_dirs} "GroupTree/Chebyshev"  "GroupTree/Core"  "GroupTree/Cuda"  "GroupTree/OpenCl"
        "GroupTree/Rotation"  "GroupTree/StarPUUtils"  "GroupTree/Taylor"  "GroupTree/TestKernel"  "GroupTree/Uniform")
endif()

FOREACH(my_dir ${my_include_dirs})	
	file(GLOB 	 	 	 	
	 	hpp_in_dir
		${my_dir}/*.hpp ${my_dir}/*.h
	)
	INSTALL( FILES ${hpp_in_dir} DESTINATION include/${my_dir} )
ENDFOREACH()

INSTALL( FILES "${CMAKE_BINARY_DIR}/Src/ScalFmmConfig.h" DESTINATION include/Utils/${my_dir} )

