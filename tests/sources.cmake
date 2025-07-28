target_sources(controlpp_tests PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/TimeVariantControl_Test.cpp
        ${CMAKE_CURRENT_LIST_DIR}/KalmanFilter_Tests.cpp
        ${CMAKE_CURRENT_LIST_DIR}/Estimators_Tests.cpp
        ${CMAKE_CURRENT_LIST_DIR}/Transformation_Tests.cpp
        ${CMAKE_CURRENT_LIST_DIR}/StateSpace_Tests.cpp
        ${CMAKE_CURRENT_LIST_DIR}/Polynom_Tests.cpp
        ${CMAKE_CURRENT_LIST_DIR}/main.cpp
)