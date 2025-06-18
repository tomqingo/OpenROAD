include("openroad")

set(TEST_NAMES
    macro_only
)

find_library(LIB_PATOH patoh /uac/gds/qluo22/disk/projects/OpenROAD-flow-scripts/tools/OpenROAD/src/mpl2/include/mpl2/patoh)
message(${LIB_PATOH})

foreach(TEST_NAME IN LISTS TEST_NAMES)
    or_integration_test("mpl2" ${TEST_NAME}  ${CMAKE_CURRENT_SOURCE_DIR}/regression)
endforeach()

add_subdirectory(cpp)



