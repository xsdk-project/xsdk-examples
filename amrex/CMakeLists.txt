# AMReX+SUNDIALS requires amrex@22.04: and sundials@6.2.0:
if(ENABLE_AMREX AND ENABLE_SUNDIALS)
    if("${AMReX_RELEASE_NUMBER}" VERSION_GREATER_EQUAL "220400") AND ("${SUNDIALS_VERSION}" VERSION_GREATER_EQUAL 6.2.0))
        add_subdirectory(sundials)
    else()
        message(STATUS "SKIPPED AMReX+SUNDIALS example because AMReX version is ${AMReX_RELEASE_NUMBER} and SUNDIALS version is ${SUNDIALS_VERSION} (need 22.04: and 6.2.0:)")
    endif()
endif()
