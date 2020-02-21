set(SOURCES_LIST ${SOURCES_LIST}
    ${CMAKE_CURRENT_LIST_DIR}/mixdata.cpp
    ${CMAKE_CURRENT_LIST_DIR}/b.cpp
    ${CMAKE_CURRENT_LIST_DIR}/dfi.cpp
    ${CMAKE_CURRENT_LIST_DIR}/sahaleft.cpp
    ${CMAKE_CURRENT_LIST_DIR}/findroot.cpp
    ${CMAKE_CURRENT_LIST_DIR}/sahamixsolver.cpp
     PARENT_SCOPE)

set(HEADERS_LIST ${HEADERS_LIST}
     ${CMAKE_CURRENT_LIST_DIR}/mixdata.h
     ${CMAKE_CURRENT_LIST_DIR}/b.h
     ${CMAKE_CURRENT_LIST_DIR}/dfi.h
     ${CMAKE_CURRENT_LIST_DIR}/sahaleft.h
     ${CMAKE_CURRENT_LIST_DIR}/findroot.h
     ${CMAKE_CURRENT_LIST_DIR}/sahamixsolver.h
     PARENT_SCOPE)

