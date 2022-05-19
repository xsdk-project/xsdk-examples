# AMReX+SUNDIALS requires amrex@22.04: and sundials@6.2.0:
if(ENABLE_SUNDIALS AND ("${AMReX_RELEASE_NUMBER}" VERSION_GREATER_EQUAL "22.04") AND ("${SUNDIALS_PACKAGE_VERSION}" VERSION_GREATER_EQUAL 6.2.0))
    add_subdirectory(sundials)
else()
    message(STATUS "SKIPPED AMReX+SUNDIALS example because AMReX version is too old or SUNDIALS version is too old (need 22.04: and 6.2.0:)")
endif()
