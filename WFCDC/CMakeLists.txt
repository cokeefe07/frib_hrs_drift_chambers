cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
project(drift_chamber_sims CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_definitions(-D_GLIBCXX_USE_CXX11_ABI=1)

# ---------------------------------------------------------------------
# Default to an optimised build.
# If CMAKE_BUILD_TYPE is empty, CMake passes NO optimisation flags at all --
# not -O2, not even -O1. For a 1.5M-drift sweep that is typically a 3-10x
# slowdown, which can easily exceed everything threading buys back. Override
# with -DCMAKE_BUILD_TYPE=Debug when you actually want symbols.
# ---------------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

find_package(GSL REQUIRED)
find_package(Garfield REQUIRED)

# Deliberately NOT REQUIRED. This is the line that lets the project still
# configure on a machine without OpenMP (e.g. stock AppleClang), falling back
# to a single-threaded build instead of failing outright.
find_package(OpenMP)

# ---------------------------------------------------------------------
# One helper so every target is wired up identically.
#
# It also enforces the rule that broke the previous file: target_link_libraries
# must be called EITHER all-plain OR all-keyword for a given target, never
# both. Calling it once plain and once with PRIVATE is a hard CMake error:
#   "The plain signature for target_link_libraries has already been used
#    with the target ... must be either all-keyword or all-plain."
# Here it is called exactly once per target, always with PRIVATE.
# ---------------------------------------------------------------------
function(add_sim_target name src)
    if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${src}")
        message(STATUS "${name}: ${src} not found - skipping this target")
        return()
    endif()

    add_executable(${name} ${src})

    if(OpenMP_CXX_FOUND)
        target_link_libraries(${name} PRIVATE Garfield::Garfield OpenMP::OpenMP_CXX)
        message(STATUS "${name}: OpenMP found - building multithreaded")
    else()
        target_link_libraries(${name} PRIVATE Garfield::Garfield)
        message(STATUS "${name}: OpenMP NOT found - building single-threaded")
    endif()
endfunction()

add_sim_target(field_cage field_cage.cc)
add_sim_target(strip_cage strip_cage.cc)


# =====================================================================
#  BUILDING
# =====================================================================
#
#   cmake -S . -B build
#   cmake --build build -j
#
# Watch the configure output for these two lines:
#   -- Build type: Release
#   -- field_cage: OpenMP found - building multithreaded
#
# If the build type line says anything else, you are running unoptimised.
#
# =====================================================================
#  RUNNING
# =====================================================================
#
#   export OMP_NUM_THREADS=16
#   ./build/field_cage
#
# The program prints how many threads it actually got, so check that line
# matches what you asked for.
#
# On a SHARED cluster node, always set OMP_NUM_THREADS explicitly. With it
# unset, OpenMP grabs every core it can see -- including cores allocated to
# other people's jobs. If you are under a scheduler, match your allocation:
#
#   SLURM:  export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
#
# On macOS with Homebrew:  brew install libomp
# and configure with:
#   cmake -S . -B build -DOpenMP_ROOT=$(brew --prefix libomp)
#
# =====================================================================
#  VERIFYING THE PARALLEL BUILD
# =====================================================================
#
#   OMP_NUM_THREADS=1  ./build/field_cage      # -> Outputs/<stamp_A>
#   OMP_NUM_THREADS=16 ./build/field_cage      # -> Outputs/<stamp_B>
#
# then diff the per-concentration .tsv files. They will NOT match bit for bit
# (AvalancheMC is Monte Carlo), but each pass_pct should agree to within the
# binomial error, sqrt(p*(1-p)/N) ~ 0.9% at N=2000 and p~0.8. Points differing
# by several percent mean something is genuinely racing -- in which case fall
# back to process-level parallelism (one CF4 concentration per process).

