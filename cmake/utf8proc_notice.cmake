# The bundled header carries authors not all listed in netCDF's COPYRIGHT or
# utf8proc's release license. Retain its complete first comment verbatim.
file(READ "${SOURCE}" header)
string(FIND "${header}" "*/" end)
if(NOT header MATCHES "^/\\*" OR end LESS 0 OR
   NOT header MATCHES "Permission is hereby granted")
  message(FATAL_ERROR "The bundled utf8proc notice could not be identified")
endif()
math(EXPR length "${end} + 2")
string(SUBSTRING "${header}" 0 ${length} notice)
file(WRITE "${OUTPUT}" "${notice}\n")
