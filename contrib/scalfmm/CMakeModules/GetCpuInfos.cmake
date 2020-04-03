###########################################################################################
# Berenger Bramas Inria
# This goes with the getCpuInfos.cpp
# This will create one CMAKE value per output option from the cpp file.
# For example the output of the CPP file can be:
# SSE3=TRUE;AVX=FALSE
# Then it will create:
# CPUOPTION_SSE3 = TRUE
# CPUOPTION_AVX = FALSE
#
# The binary should return 0 on success.
###########################################################################################
macro(GetCpuInfos)
# The original CPP file
set(GetCpuInfosFile "${PROJECT_SOURCE_DIR}/CMakeModules/getCpuInfos.cpp")

# Fatal error if the file does not exist
if(NOT EXISTS ${GetCpuInfosFile})
	message(FATAL_ERROR "The GetCpuInfosFile does not exist (${GetCpuInfosFile})")
endif()

# Compile and execute the file
try_run(RUN_RESULT_VAR COMPILE_RESULT_VAR
          ${CMAKE_BINARY_DIR} ${GetCpuInfosFile}  # [CMAKE_FLAGS <Flags>] [COMPILE_DEFINITIONS <flags>]
          COMPILE_OUTPUT_VARIABLE comp
          RUN_OUTPUT_VARIABLE run)

# If it has successfuly compiled an run
if(COMPILE_RESULT_VAR AND (RUN_RESULT_VAR EQUAL 0) )
	set( CPU_OPTIONS ${run} )
	# For each value
	foreach(optionNode ${run})
		# Get name and value
		string(REPLACE "=" ";" optionNameAndValue ${optionNode})
		list(LENGTH optionNameAndValue optionLength)
		# If we get both
		if(optionLength EQUAL 2)
			list(GET optionNameAndValue 0 optionName)
			list(GET optionNameAndValue 1 optionValue)
			# create cmake variable
			set(CPUOPTION_${optionName} ${optionValue})
		else()
			message(WARNING "GetCpuInfosFile wrong format for ${optionNode}.")
		endif()
	endforeach()
	# output the sentence from the binrary
	message(STATUS "CPUOPTION : ${CPU_OPTIONS}")
else()
	message(WARNING "GetCpuInfosFile did not return correctly.")
endif()

endmacro(GetCpuInfos)
