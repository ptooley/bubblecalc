set (BUBBLECALCSOURCES
  main
  electron
  bubble
)

add_executable(bubblecalc ${BUBBLECALCSOURCES})
set_target_properties(bubblecalc PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(bubblecalc ${EXTRA_LIBS})
