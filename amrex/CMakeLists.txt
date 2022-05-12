
# AMReX+SUNDIALS requires amrex@22.04:
if(ENABLE_SUNDIALS AND ("${AMReX_RELEASE_NUMBER}" VERSION_GREATER_EQUAL "22.04"))
    add_subdirectory(sundials)
else()
    message(STATUS "SKIPPED AMReX+SUNDIALS example because AMReX version is too old")
endif()

