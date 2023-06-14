cmake_minimum_required(VERSION 3.16.3...3.19.7 FATAL_ERROR)

project(BiteID)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "https://github.com/agporto/BiteID")
set(EXTENSION_CATEGORY "Examples")
set(EXTENSION_CONTRIBUTORS "Arthur Porto")
set(EXTENSION_DESCRIPTION "Testing ID software")
set(EXTENSION_ICONURL "https://www.example.com/Slicer/Extensions/BiteID.png")
set(EXTENSION_SCREENSHOTURLS "https://www.example.com/Slicer/Extensions/BiteID/Screenshots/1.png")
set(EXTENSION_DEPENDS "NA") # Specified as a list or "NA" if no dependencies

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(BiteID)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})