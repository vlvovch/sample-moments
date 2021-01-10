# Now simply link against gtest or gtest_main as needed. Eg
add_executable(test_PartitionsOfSet test_PartitionsOfSet.cpp)
target_link_libraries(test_PartitionsOfSet sample_statistics gtest_main)
set_property(TARGET test_PartitionsOfSet PROPERTY FOLDER tests)
add_test(NAME PartitionsOfSet COMMAND test_PartitionsOfSet)
