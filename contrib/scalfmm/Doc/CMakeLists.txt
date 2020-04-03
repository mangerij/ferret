# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(
        doc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen" VERBATIM
    )
    # INSTALL(FILES   ${SCALFMM_BINARY_DIR}/Doc/scalfmm.tag DESTINATION doc/  )
    # INSTALL(DIRECTORY   ${SCALFMM_BINARY_DIR}/Doc/html  DESTINATION doc/  )
else()
    message( WARNING "You ask to enable the doc generation but Doxygen cannot be found." )    
endif(DOXYGEN_FOUND)
