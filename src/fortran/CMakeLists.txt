add_executable(${ExampleName}_Fortran ${ExampleName}.F90)

target_link_libraries(${ExampleName}_Fortran PRIVATE OpenCMISS::libOpenCMISS_Fortran)

install(TARGETS ${ExampleName}_Fortran DESTINATION bin)

configure_file(visualise.cmgui visualise.cmgui)
