FILE (GLOB src ${CMAKE_CURRENT_SOURCE_DIR}/*.cc)
FILE (GLOB cfg ${CMAKE_CURRENT_SOURCE_DIR}/*.x*)

SET (vd_reconstruction_sources ${src} PARENT_SCOPE)
SET (vd_reconstruction_data ${cfg} PARENT_SCOPE)
