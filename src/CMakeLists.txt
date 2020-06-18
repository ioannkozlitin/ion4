add_subdirectory(mixdata)
add_subdirectory(raizermixsolver)

set(SOURCES_LIST ${SOURCES_LIST}
    ${CMAKE_CURRENT_LIST_DIR}/main.cpp
    PARENT_SCOPE)

set(HEADERS_LIST ${HEADERS_LIST}
    PARENT_SCOPE)