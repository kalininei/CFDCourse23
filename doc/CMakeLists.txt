# Documentation build
configure_file(${CMAKE_SOURCE_DIR}/doc/Doxyfile ${CMAKE_BINARY_DIR}/doc/Doxyfile)
configure_file(${CMAKE_SOURCE_DIR}/doc/MainPage.md ${CMAKE_BINARY_DIR}/doc/MainPage.md)

add_custom_target(doc
	COMMAND doxygen Doxyfile
	WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/doc
	COMMENT "Generating API documentation with Doxygen"
	VERBATIM)

