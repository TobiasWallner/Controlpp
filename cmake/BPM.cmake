cmake_minimum_required(VERSION 3.21)

function(bpm_load_env_var VAR_NAME DEFAULT)
    # VAR_NAME: the variable name to check (e.g., "BPM_CLEAN_INSTALL")
    # DEFAULT: the default value if not found
    
    if(DEFINED ${VAR_NAME})
        # use local CMake variable (highest priority)
        set(_value "${${VAR_NAME}}")
    elseif(DEFINED ENV{${VAR_NAME}})
        # fall back to environment variable
        set(_value "$ENV{${VAR_NAME}}")
    else()
        # fall back to default
        set(_value "${DEFAULT}")
    endif()
    
    set(${VAR_NAME} "${_value}" PARENT_SCOPE)
endfunction()

# @parses a constraint + version string into a closed-open range of allowed version
# 
# input: >=1.2.3 --> output: LIST 1.2.3;inf.inf.inf -- meaning: from version 1.2.3 to upper bound inf.inf.inf
# input: ^1.2.3 --> output: LIST 1.2.3;2.0.0 -- meaning from version 1.2.3 to upper bound 2.0.0
function(bpm_parse_version_string in_pkg_name INPUT out_version_range)
    string(REPLACE ";" "\\;" SAFE_INPUT "${INPUT}")

    set(VERSION_QUALIFIER "")
    set(VERSION_MAJOR "")
    set(VERSION_MINOR "")
    set(VERSION_PATCH "")

    string(REGEX MATCH "^[ \t]*(>=|\\^|~|=)?[ \t]*v?([0-9]+\\.[0-9]+\\.[0-9]+)(.*(<)[ \t]*v?([0-9]+\\.[0-9]+\\.[0-9]+)[^<]*)?[ \t]*$" _match_result "${SAFE_INPUT}")

    # message(STATUS "_match_result: ${_match_result}")
    # message(STATUS "CMAKE_MATCH_0: ${CMAKE_MATCH_0}")
    # message(STATUS "CMAKE_MATCH_1: ${CMAKE_MATCH_1}")
    # message(STATUS "CMAKE_MATCH_2: ${CMAKE_MATCH_2}")
    # message(STATUS "CMAKE_MATCH_3: ${CMAKE_MATCH_3}")
    # message(STATUS "CMAKE_MATCH_4: ${CMAKE_MATCH_4}")
    # message(STATUS "CMAKE_MATCH_5: ${CMAKE_MATCH_5}")

    if(_match_result)
        # has a semver version qualifier
        if(CMAKE_MATCH_1)
            set(version_qualifier "${CMAKE_MATCH_1}")
        else()
            set(version_qualifier "=")
        endif()
        set(lower_bound_version "${CMAKE_MATCH_2}")
        set(upper_bound_version "${CMAKE_MATCH_5}")

        if(lower_bound_version)
            if(upper_bound_version)
                if(NOT "${lower_bound_version}" VERSION_LESS "${upper_bound_version}")
                    message(FATAL_ERROR "BPM [${PROJECT_NAME}:${in_pkg_name}]: Invalid version constraint: ${version_qualifier}${lower_bound_version} <${upper_bound_version}. The lower bound (${lower_bound_version}) has to be strictly less than the uper bound ${upper_bound_version}")
                endif()
            endif()
        endif()

        string(REGEX MATCH "^([0-9]+)\\.([0-9]+)\\.([0-9]+)" _match_result "${lower_bound_version}")
        if(_match_result)
            set(major_lower "${CMAKE_MATCH_1}")
            set(minor_lower "${CMAKE_MATCH_2}")
            set(patch_lower "${CMAKE_MATCH_3}")

            if(version_qualifier STREQUAL ">=")
                set(major_upper "inf")
                set(minor_upper "inf")
                set(patch_upper "inf")
            elseif(version_qualifier STREQUAL "^")
                if(major_lower EQUAL 0)
                    set(major_upper "0")
                    math(EXPR minor_upper "${minor_lower} + 1")
                    set(patch_upper "0")
                else()
                    math(EXPR major_upper "${major_lower} + 1")
                    set(minor_upper "0")
                    set(patch_upper "0")
                endif()
            elseif(version_qualifier STREQUAL "~")
                set(major_upper "${major_lower}")
                math(EXPR minor_upper "${minor_lower} + 1")
                set(patch_upper "0")
            else()
                set(major_upper ${major_lower})
                set(minor_upper ${minor_lower})
                math(EXPR patch_upper "${patch_lower} + 1") # upper bound is non inclusive
            endif()
            
            set(out_version_range_0 "${major_lower}.${minor_lower}.${patch_lower}")
            set(out_version_range_1 "${major_upper}.${minor_upper}.${patch_upper}")

            if(BPM_VERBOSE)
                message(STATUS "BPM [${PROJECT_NAME}:${in_pkg_name}]: range from lowerbound constraint '${version_qualifier}${lower_bound_version}' --> ${out_version_range_0}-${out_version_range_1}")
            endif()

            if(upper_bound_version)
                if("${out_version_range_1}" STREQUAL "inf.inf.inf")
                    set(out_version_range_1 "${upper_bound_version}")
                elseif("${upper_bound_version}" VERSION_LESS out_version_range_1)
                    set(out_version_range_1 "${upper_bound_version}")
                endif()

                if(BPM_VERBOSE)
                    message(STATUS "BPM [${PROJECT_NAME}:${in_pkg_name}]: range from lower- & upper-bound constraints '${version_qualifier}${lower_bound_version} <${upper_bound_version}' --> ${out_version_range_0}-${out_version_range_1}")
                endif()
            endif()
            
            # Valid semver
            set(${out_version_range} "${out_version_range_0}" "${out_version_range_1}" PARENT_SCOPE)
            return()
        else()
            message(FATAL_ERROR "BPM [${PROJECT_NAME}:${in_pkg_name}]: lower_bound_version (${lower_bound_version}) not semver")
        endif()
    else()
        # no semver version qualifier --> assume git-tag or git-commit-hash
        string(REGEX MATCH "^[ \t]*((>=|\\^|~|=)?[ \t]*(.*[^ \t]))[ \t]*$" _match_result "${SAFE_INPUT}")
        
        # message(STATUS "_match_result: ${_match_result}")
        # message(STATUS "CMAKE_MATCH_0: ${CMAKE_MATCH_0}")
        # message(STATUS "CMAKE_MATCH_1: ${CMAKE_MATCH_1}")
        # message(STATUS "CMAKE_MATCH_2: ${CMAKE_MATCH_2}")
        # message(STATUS "CMAKE_MATCH_3: ${CMAKE_MATCH_3}")
        
        if(_match_result)
            set(git_tag "${CMAKE_MATCH_1}")

            set(version_qualifier "${CMAKE_MATCH_2}")
            if(version_qualifier)
                message(WARNING "BPM [${PROJECT_NAME}:${in_pkg_name}]: Version constraints (like '${version_qualifier}') are not supported for named-git-tags or git-commit-hashes '${git_tag}'.")
            endif()

            string(REGEX MATCH "[ \t\r\n]" has_ws "${git_tag}")
            if(has_ws)
                message(FATAL_ERROR "BPM [${PROJECT_NAME}:${in_pkg_name}]: Git-tag or commit-hash '${git_tag}' cannot contain whitespaces.")
            endif()

            set(${out_version_range} "${git_tag}" PARENT_SCOPE)
            return()
        else()
            message(FATAL_ERROR "BPM [${PROJECT_NAME}:${in_pkg_name}]: Could not parse version string ${INPUT}")
        endif()
    endif()

    # --------------------------------------------------------
    # Split qualifier from remainder
    # --------------------------------------------------------
    string(REGEX MATCH "^(>=|\\^|~|=)?(.+)$" _ "${SAFE_INPUT}")
    if(CMAKE_MATCH_1)
        set(VERSION_QUALIFIER "${CMAKE_MATCH_1}")
    else()
        set(VERSION_QUALIFIER "=")
    endif()
    set(VALUE     "${CMAKE_MATCH_2}")

    # --------------------------------------------------------
    # Try semantic version
    # --------------------------------------------------------
    string(REGEX MATCH "^([0-9]+)\\.([0-9]+)\\.([0-9]+)" _ "${SAFE_INPUT}")

    if(CMAKE_MATCH_0)
        set(major_lower "${CMAKE_MATCH_1}")
        set(minor_lower "${CMAKE_MATCH_2}")
        set(patch_lower "${CMAKE_MATCH_3}")

        if(VERSION_QUALIFIER STREQUAL ">=")
            set(major_upper "inf")
            set(minor_upper "inf")
            set(patch_upper "inf")
        elseif(VERSION_QUALIFIER STREQUAL "^")
            if(major_lower EQUAL 0)
                set(major_upper "0")
                math(EXPR minor_upper "${minor_lower} + 1")
                set(patch_upper "0")
            else()
                math(EXPR major_upper "${major_lower} + 1")
                set(minor_upper "0")
                set(patch_upper "0")
            endif()
        elseif(VERSION_QUALIFIER STREQUAL "~")
            set(major_upper "${major_lower}")
            math(EXPR minor_upper "${minor_lower} + 1")
            set(patch_upper "0")
        else()
            set(major_upper ${major_lower})
            set(minor_upper ${minor_lower})
            math(EXPR patch_upper "${patch_lower} + 1") # upper bound is non inclusive
        endif()

        # Valid semver
        set(${out_version_range} 
            "${major_lower}.${minor_lower}.${patch_lower}" 
            "${major_upper}.${minor_upper}.${patch_upper}"
            PARENT_SCOPE)
        return()
    else()
        # check for user input errors
        if(NOT (VERSION_QUALIFIER STREQUAL ">=" OR VERSION_QUALIFIER STREQUAL "^" OR VERSION_QUALIFIER STREQUAL "~" OR VERSION_QUALIFIER STREQUAL "="))
            message(FATAL_ERROR "BPM [${PROJECT_NAME}:${NAME}]: Invalid qualifier '${VERSION_QUALIFIER}'. Allowed are: '>=', '^', '~', '=' and ''")
        endif()

        # ------------------ ----------------------------------
        # Not semver → treat as tag or commit hash
        # ----------------------------------------------------

        # Disallow ^ and ~ for non-semver
        if(VERSION_QUALIFIER STREQUAL ">=" OR VERSION_QUALIFIER STREQUAL "^" OR VERSION_QUALIFIER STREQUAL "~")
            message(FATAL_ERROR "BPM [${PROJECT_NAME}:${NAME}]: Invalid constraint '${VERSION_QUALIFIER}' for non-semver reference '${VALUE}'")
        endif()

        set(${out_version_range} "${VALUE}" PARENT_SCOPE)
        return()
    endif()
endfunction()


function(bpm_parse_short_dependency INPUT out_git_repo out_name out_tag)

    # ------------------------------------------------------------
    # Split into FULL_PATH and optional VERSION_PART using '#'
    # ------------------------------------------------------------
    string(REGEX MATCH "^([^#]+)(#(.+))?$" _ "${INPUT}")

    set(FULL_PATH "${CMAKE_MATCH_1}")
    set(VERSION_PART "${CMAKE_MATCH_3}")

    # ------------------------------------------------------------
    # Extract repository name (after last '/' or '\')
    # ------------------------------------------------------------
    #string(REGEX MATCH "([^/\\\\]+)\\.git$" _ "${FULL_PATH}")
    string(REGEX MATCH "^.*[/\\\\]([^/\\\\]+)[/\\\\]?$" _ "${FULL_PATH}")
    set(NAME "${CMAKE_MATCH_1}")
    if(BPM_VERBOSE)
        message(STATUS "${INPUT} --> Name: ${NAME}")
    endif()

    string(REGEX REPLACE "\\.git[ \t]*$" "" NAME "${NAME}")
    if(BPM_VERBOSE)
        message(STATUS "${INPUT} --> Name (removed '.git'): ${NAME}")
    endif()

    if(NOT NAME)
        message(FATAL_ERROR "BPM [${PROJECT_NAME}]: Could not extract the repository name. Expected: 'path/name' or `NAME ... GIT_REPOSITORY ... GIT_TAG ...` but got: ${FULL_PATH}")
    endif()

    if(NOT VERSION_PART)
        message(FATAL_ERROR "BPM [${PROJECT_NAME}:${NAME}]: No version, git-tag or git-commit-hash provided. Expected: 'path/name#version/tag/hash'")
    endif()

    # ------------------------------------------------------------
    # Export results
    # ------------------------------------------------------------
          
    set(${out_git_repo} "${FULL_PATH}" PARENT_SCOPE)
    set(${out_name} "${NAME}" PARENT_SCOPE)
    set(${out_tag} "${VERSION_PART}" PARENT_SCOPE)

endfunction()

function(bpm_parse_arguments INPUT out_name out_repo out_tag out_options out_packages out_version out_private)
    # clear variables to prevent accidental reuse in the loop
    unset(PKG_NAME)
    unset(PKG_GIT_REPOSITORY)
    unset(PKG_GIT_TAG)
    unset(PKG_PACKAGES)
    unset(PKG_OPTIONS)
    unset(PKG_VERSION)
    unset(PKG_VERSION_QUALIFIER)
    unset(PKG_VERSION_RANGE)
    unset(PKG_INGNORE)

    # Parse arguments in long form
    set(options PRIVATE)
    set(oneValueArgs NAME GIT_REPOSITORY GIT_TAG)
    set(multiValueArgs PACKAGES OPTIONS DEPENDENCIES)
    cmake_parse_arguments(PKG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${INPUT})

    foreach(opt IN LISTS PKG_OPTIONS)
        if(NOT "${opt}" MATCHES "^[^=]+=.+$")
            message(FATAL_ERROR "BPM [${PROJECT_NAME}:${PKG_NAME}]: Option parsing error. Option '${opt}' does not parse to '<name>=<value>'.")
        endif()
        if("${opt}" MATCHES "^-D")
            message(WARNING "BPM [${PROJECT_NAME}:${PKG_NAME}]: Option '${opt}' startis with '-D'. For install packages '-D' will be automatically added.")
        endif()
    endforeach()

    # Parse arguments in short form
    if(NOT PKG_NAME AND NOT PKG_GIT_REPOSITORY AND NOT PKG_GIT_TAG)
        list(GET INPUT 0 FIRST_ARG)
        bpm_parse_short_dependency(${FIRST_ARG} PKG_GIT_REPOSITORY PKG_NAME PKG_GIT_TAG)
    endif()
        
    # Validate arguments
    if(NOT PKG_NAME)
        message(FATAL_ERROR "BPM [${PROJECT_NAME}]: NAME is required")
    endif()

    if(BPM_${PKG_NAME}_ADDED)
        message(STATUS "${BPM_NAME} was already added")
        return()
    endif()

    if(NOT PKG_PACKAGES)
        set(PKG_PACKAGES "${PKG_NAME}")
    endif()
    
    if(NOT PKG_GIT_REPOSITORY)
        message(FATAL_ERROR "BPM [${PROJECT_NAME}:${PKG_NAME}]: GIT_REPOSITORY is required")
    endif()

    if(NOT PKG_GIT_TAG)
        message(FATAL_ERROR "BPM [${PROJECT_NAME}:${PKG_NAME}]: git-tag or version constraint is required")
    else()
        bpm_parse_version_string("${PKG_NAME}" "${PKG_GIT_TAG}" PKG_VERSION)
    endif()
      
    if(PKG_PRIVATE)
        set(${out_private} TRUE PARENT_SCOPE)
    else()
        set(${out_private} FALSE PARENT_SCOPE)
    endif()

    set(${out_name} ${PKG_NAME} PARENT_SCOPE)
    set(${out_repo} ${PKG_GIT_REPOSITORY} PARENT_SCOPE)
    set(${out_tag} ${PKG_GIT_TAG} PARENT_SCOPE)
    set(${out_version_qualifier} ${PKG_VERSION_QUALIFIER} PARENT_SCOPE)
    set(${out_version} ${PKG_VERSION} PARENT_SCOPE)
    set(${out_options} ${PKG_OPTIONS} PARENT_SCOPE)
    set(${out_packages} ${PKG_PACKAGES} PARENT_SCOPE)

endfunction()

function(bpm_version_range_intersection in_version_range_a in_version_range_b out_version_range)
    list(GET in_version_range_a 0 a_lower)
    list(GET in_version_range_a 1 a_upper)

    list(GET in_version_range_b 0 b_lower)
    list(GET in_version_range_b 1 b_upper)


    if(a_lower VERSION_LESS b_lower)
        set(lower_bound ${b_lower})
    else()
        set(lower_bound ${a_lower})
    endif()

    if((a_upper STREQUAL "inf.inf.inf") AND (b_upper STREQUAL "inf.inf.inf"))
        set(upper_bound "inf.inf.inf")
    elseif((NOT a_upper STREQUAL "inf.inf.inf") AND (b_upper STREQUAL "inf.inf.inf"))
        set(upper_bound "${a_upper}")
    elseif((a_upper STREQUAL "inf.inf.inf") AND (NOT b_upper  STREQUAL "inf.inf.inf"))
        set(upper_bound "${b_upper}")
    elseif(a_upper VERSION_LESS b_upper)
        set(upper_bound "${a_upper}")
    else()
        set(upper_bound "${b_upper}")
    endif()

    if((NOT upper_bound STREQUAL "inf.inf.inf") AND (upper_bound VERSION_LESS lower_bound))
        # version conflict: return without result
        return()
    endif()

    set(${out_version_range} "${lower_bound}" "${upper_bound}" PARENT_SCOPE)
endfunction()

function(path_normalise INPUT OUTPUT)
    set(path "${INPUT}")

    # normalize slashes
    file(TO_CMAKE_PATH "${path}" path)

    # remove trailing whitespace
    string(STRIP "${path}" path)

    # remove trailing .git
    string(REGEX REPLACE "\\.git$" "" path "${path}")

    # remove trailing slash
    string(REGEX REPLACE "/+$" "" path "${path}")

    # Decide whether this is a local filesystem path.
    set(is_filesystem_path FALSE)
    if(EXISTS "${path}")
        set(is_filesystem_path TRUE)
    elseif(IS_ABSOLUTE "${path}")
        set(is_filesystem_path TRUE)
    elseif(path MATCHES "^[.][.]?(/.*)?$")
        set(is_filesystem_path TRUE)
    elseif(path MATCHES "^[A-Za-z]:/.*$")
        set(is_filesystem_path TRUE)
    elseif(path MATCHES "^/.*$")
        set(is_filesystem_path TRUE)
    endif()

    if(is_filesystem_path)
        file(REAL_PATH "${path}" path BASE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
        
        string(REGEX REPLACE "/+$" "" path "${path}")
        
        # Windows paths should compare case-insensitively
        if(WIN32)
            string(TOLOWER "${path}" path)
        endif()

        set(${OUTPUT} "${path}" PARENT_SCOPE)
        return()
    else()
        # ssh form: git@host:user/repo
        if(path MATCHES "^[^@]+@([^:]+):(.+)$")
            string(REGEX REPLACE "^[^@]+@([^:]+):(.+)$" "web://\\1/\\2" path "${path}")
        endif()

        # ssh://git@host/user/repo
        if(path MATCHES "^ssh://")
            string(REGEX REPLACE "^ssh://([^@]+@)?" "web://" path "${path}")
        endif()

        # http(s)://host/user/repo
        if(path MATCHES "^http[s]?://")
            string(REGEX REPLACE "^http[s]?://" "web://" path "${path}")
        endif()

        # normalize slashes in case of Windows paths
        file(TO_CMAKE_PATH "${path}" path)

        set(${OUTPUT} "${path}" PARENT_SCOPE)
        return()
    endif()
endfunction()

function(path_equal A B RESULT)
    path_normalise("${A}" A_NORM)
    path_normalise("${B}" B_NORM)

    if(A_NORM STREQUAL B_NORM)
        set(${RESULT} TRUE PARENT_SCOPE)
    else()
        set(${RESULT} FALSE PARENT_SCOPE)
    endif()
endfunction()

function(bpm_is_path INPUT RESULT)
    if(EXISTS "${INPUT}")
        set(${RESULT} TRUE PARENT_SCOPE)
    elseif(IS_ABSOLUTE "${INPUT}")
        set(${RESULT} TRUE PARENT_SCOPE)
    elseif(INPUT MATCHES "^[.][.]?(/.*)?$")
        set(${RESULT} TRUE PARENT_SCOPE)
    elseif(INPUT MATCHES "^[A-Za-z]:/.*$")
        set(${RESULT} TRUE PARENT_SCOPE)
    elseif(INPUT MATCHES "^/.*$")
        set(${RESULT} TRUE PARENT_SCOPE)
    elseif(INPUT MATCHES "^http[s]?://")
        set(${RESULT} TRUE PARENT_SCOPE)
    elseif(INPUT MATCHES "^ssh://")
        set(${RESULT} TRUE PARENT_SCOPE)
    else()
        set(${RESULT} FALSE PARENT_SCOPE)
    endif()
endfunction()

function(bpm_equal_option_value A B OUT)
    if("${A}" STREQUAL "${B}") # check if the have the same value
        set(${OUT} TRUE PARENT_SCOPE)
        return()
    elseif((("${A}" STREQUAL "TRUE") OR ("${A}" STREQUAL "ON")) AND (("${B}" STREQUAL "TRUE") OR ("${B}" STREQUAL "ON")))
        set(${OUT} TRUE PARENT_SCOPE)
        return()
    elseif((("${A}" STREQUAL "FALSE") OR ("${A}" STREQUAL "OFF")) AND (("${B}" STREQUAL "FALSE") OR ("${B}" STREQUAL "OFF")))
        set(${OUT} TRUE PARENT_SCOPE)
        return()
    else()
        # check if it might be a path, compare them with path-aware comparison
        bpm_is_path("${A}" is_A_path)
        bpm_is_path("${B}" is_B_path)
        if(is_A_path AND is_B_path)    
            # options might be paths, compare them with path-aware comparison
            path_equal("${A}" "${B}" eq)
            if(eq)
                set(${OUT} TRUE PARENT_SCOPE)
                return()
            else()
                set(${OUT} FALSE PARENT_SCOPE)
                return()
            endif()
        else()
            set(${OUT} FALSE PARENT_SCOPE)
            return()
        endif()
    endif()
endfunction()


function(bpm_combine_options PKG_OPTIONS SEL_OPTIONS PKG_REQUIRED_FROM SEL_REQUIRED_FROM options_out has_added_out)

    set(combined_options "${SEL_OPTIONS}")
    set(added_options FALSE)
    foreach(pkg_opt IN LISTS PKG_OPTIONS)
        string(REGEX MATCH "^([^=]+)=(.*)$" _match "${pkg_opt}")
        if(_match)
            set(pkg_opt_name ${CMAKE_MATCH_1})
            set(pkg_opt_value ${CMAKE_MATCH_2})
        else()
            set(msg "Option error: '${pkg_opt}' required from '${PKG_REQUIRED_FROM}' does not follow '<name>=<value>'")
            string(APPEND logging "\n  ${msg}")
            set(filename "${CMAKE_BINARY_DIR}/bpm-dependency-solver-log.txt")
            file(WRITE "${filename}" "${logging}")  
            message(FATAL_ERROR "BPM [${PROJECT_NAME}:${PKG_NAME}]: ${msg}. See '${filename}' for details")
        endif()

        set(contains FALSE)
        foreach(sel_opt IN LISTS SEL_OPTIONS)
            string(REGEX MATCH "^([^=]+)=(.*)$" _match "${sel_opt}")
            if(_match)
                set(sel_opt_name ${CMAKE_MATCH_1})
                set(sel_opt_value ${CMAKE_MATCH_2})
            else()
                set(msg "Options error: '${sel_opt}' required from '${SEL_REQUIRED_FROM}' does not follow '<name>=<value>'")
                string(APPEND logging "\n  ${msg}")
                set(filename "${CMAKE_BINARY_DIR}/bpm-dependency-solver-log.txt")
                file(WRITE "${filename}" "${logging}")  
                message(FATAL_ERROR "BPM [${PROJECT_NAME}:${PKG_NAME}]: ${msg}. See '${filename}' for details")
            endif()

            # with the same name
            if("${pkg_opt_name}" STREQUAL "${sel_opt_name}")
                bpm_equal_option_value("${sel_opt_value}" "${pkg_opt_value}" eq)
                if(NOT eq)
                    set(msg "Options conflict: '${sel_opt}' required from '${SEL_REQUIRED_FROM}' != '${pkg_opt}' required from '${PKG_REQUIRED_FROM}'")
                    string(APPEND logging "\n  ${msg}")
                    set(filename "${CMAKE_BINARY_DIR}/bpm-dependency-solver-log.txt")
                    file(WRITE "${filename}" "${logging}")
                    message(FATAL_ERROR "BPM [${PROJECT_NAME}:${PKG_NAME}]: ${msg}. See '${filename}' for details")
                endif()
                set(contains TRUE)
            endif()

        endforeach()
        if(NOT contains)
            set(msg "Added option '${pkg_opt}' to package '${PKG_NAME}' required from '${PKG_REQUIRED_FROM}'")
            string(APPEND logging "\n  ${msg}")
            if(BPM_VERBOSE)
                message(STATUS "  ${msg}")
            endif()
            
            list(APPEND combined_options "${pkg_opt}")
            set(added_options TRUE)
        endif()
    endforeach()

    set(${options_out} "${combined_options}" PARENT_SCOPE)
    set(${has_added_out} ${added_options} PARENT_SCOPE)

endfunction()

function(not_equal_bool A B OUT)
    if(A)
        if(B)
            set(${OUT} FALSE PARENT_SCOPE)
        else()
            set(${OUT} TRUE PARENT_SCOPE)
        endif()
    else()
        if(B)
            set(${OUT} TRUE PARENT_SCOPE)
        else()
            set(${OUT} FALSE PARENT_SCOPE)
        endif()
    endif()
endfunction()

function(bpm_add_package_to_registry PKG_NAME PKG_GIT_REPOSITORY PKG_GIT_TAG PKG_OPTIONS PKG_PACKAGES PKG_VERSION_RANGE PKG_PRIVATE TYPE)
    if(NOT (("${TYPE}" STREQUAL "INSTALL") OR ("${TYPE}" STREQUAL "ADD_SUBDIR")))
        message(FATAL_ERROR "BPM [${PROJECT_NAME}:${PKG_NAME}]: Internal error: TYPE (${TYPE}) should be 'INSTALL' or 'ADD_SUBDIR'")
    endif()

    get_property(BPM_REGISTRY_ GLOBAL PROPERTY BPM_REGISTRY)
    if(BPM_REGISTRY_)
        # uniquely add the package name to the list
        get_property("BPM_${PKG_NAME}_ADDED_" GLOBAL PROPERTY "BPM_${PKG_NAME}_ADDED")
        if(NOT BPM_${PKG_NAME}_ADDED_)
            list(APPEND BPM_REGISTRY_ "${PKG_NAME}")
            set_property(GLOBAL PROPERTY BPM_REGISTRY "${BPM_REGISTRY_}")
        endif()
    else()
        set_property(GLOBAL PROPERTY BPM_REGISTRY "${PKG_NAME}")
    endif()

    get_property(BPM_${PKG_NAME}_ADDED_ GLOBAL PROPERTY BPM_${PKG_NAME}_ADDED)
    if(NOT BPM_${PKG_NAME}_ADDED_)
        set_property(GLOBAL PROPERTY "BPM_REGISTRY_${PKG_NAME}_VERSION_RANGE" "${PKG_VERSION_RANGE}")
        set_property(GLOBAL PROPERTY "BPM_REGISTRY_${PKG_NAME}_GIT_TAG" "${PKG_GIT_TAG}") # TODO rename to constraint
        set_property(GLOBAL PROPERTY "BPM_REGISTRY_${PKG_NAME}_GIT_REPOSITORY" "${PKG_GIT_REPOSITORY}")
        set_property(GLOBAL PROPERTY "BPM_REGISTRY_${PKG_NAME}_PACKAGES" "${PKG_PACKAGES}")
        set_property(GLOBAL PROPERTY "BPM_REGISTRY_${PKG_NAME}_PRIVATE" "${PKG_PRIVATE}")
        set_property(GLOBAL PROPERTY "BPM_REGISTRY_${PKG_NAME}_TYPE" "${TYPE}")

        # sort options
        list(SORT PKG_OPTIONS)
        set_property(GLOBAL PROPERTY "BPM_REGISTRY_${PKG_NAME}_OPTIONS" "${PKG_OPTIONS}")
    else()
        message(WARNING "BPM [${PROJECT_NAME}:${PKG_NAME}]: Package added twice in the same project")

        get_property(REGISTERED_GIT_REPOSITORY GLOBAL PROPERTY "BPM_REGISTRY_${PKG_NAME}_GIT_REPOSITORY")
        if(NOT "${REGISTERED_GIT_REPOSITORY}" STREQUAL "${PKG_GIT_REPOSITORY}")
            message(FATAL_ERROR "BPM [${PROJECT_NAME}:${PKG_NAME}}: Repository Conflict: new: ${PKG_GIT_REPOSITORY}, previously defined: ${REGISTERED_GIT_REPOSITORY}")
        endif()

        get_property(REGISTERED_VERSION_RANGE GLOBAL PROPERTY "BPM_REGISTRY_${PKG_NAME}_VERSION_RANGE")

        bpm_version_range_intersection("${REGISTERED_VERSION_RANGE}" "${PKG_VERSION_RANGE}" intersec_version_range)
        if(NOT intersec_version_range)
            message(FATAL_ERROR 
                "BPM ${PKG_NAME}: Version Conflict, new: ${PKG_VERSION_RANGE}, previously defined: ${REGISTERED_VERSION_RANGE}")
        endif()

        set_property(GLOBAL PROPERTY BPM_REGISTRY_${PKG_NAME}_VERSION_RANGE "${intersec_version_range}")

        # combine options, remove duplicates, check if there are conflicting ones
        list(SORT PKG_OPTIONS)
        get_property(REGISTERED_OPTIONS GLOBAL PROPERTY "BPM_REGISTRY_${PKG_NAME}_OPTIONS")
        bpm_combine_options("${PKG_OPTIONS}" "${REGISTERED_OPTIONS}" "${PROJECT_NAME}" "${PROJECT_NAME}" combined_options _)
        set_property(GLOBAL PROPERTY "BPM_REGISTRY_${PKG_NAME}_OPTIONS" "${combined_options}")

        # resolve packages
        get_property(REGISTERED_PACKAGES GLOBAL PROPERTY "BPM_REGISTRY_${PKG_NAME}_PACKAGES")
        set(JOINED_PACKAGES ${PKG_PACKAGES} ${REGISTERED_PACKAGES})
        list(REMOVE_DUPLICATES JOINED_PACKAGES)
        set_property(GLOBAL PROPERTY "BPM_REGISTRY_${PKG_NAME}_PACKAGES" "${JOINED_PACKAGES}")

        # check if there is a type conflict
        get_property(REGISTERED_TYPE GLOBAL PROPERTY "BPM_REGISTRY_${PKG_NAME}_TYPE")
        if(NOT "${REGISTERED_TYPE}" STREQUAL "${TYPE}")
            message(FATAL_ERROR "BPM [${PROJECT_NAME}:${PKG_NAME}]: Package added, once as 'INSTALL' (find_package) and once as 'SOURCE' (add_subdirectory). Decide on one to prevent target-name conflicts.")
        endif()

        # check if there is an ignore conflict
        get_property(REGISTERED_PRIVATE GLOBAL PROPERTY "BPM_REGISTRY_${PKG_NAME}_PRIVATE")

        not_equal_bool("${REGISTERED_PRIVATE}" "${PKG_PRIVATE}" private_conflict)
        if(private_conflict)
            message(FATAL_ERROR "BPM [${PROJECT_NAME}:${PKG_NAME}]: Package added, once with 'PRIVATE' and once without.")
        endif()

    endif()

    set_property(GLOBAL PROPERTY "BPM_${PKG_NAME}_ADDED" TRUE)
endfunction()


function(BPMAddInstallPackage)

    if(NOT PROJECT_IS_TOP_LEVEL)
        return()
    endif()

    bpm_parse_arguments("${ARGN}"
        PKG_NAME PKG_GIT_REPOSITORY PKG_GIT_TAG 
        PKG_OPTIONS PKG_PACKAGES PKG_VERSION_RANGE PKG_PRIVATE)

    bpm_add_package_to_registry(
        "${PKG_NAME}" "${PKG_GIT_REPOSITORY}" "${PKG_GIT_TAG}" "${PKG_OPTIONS}" 
        "${PKG_PACKAGES}" "${PKG_VERSION_RANGE}" "${PKG_PRIVATE}" "INSTALL")

endfunction()

function(BPMAddSourcePackage)

    if(NOT PROJECT_IS_TOP_LEVEL)
        return()
    endif()

    bpm_parse_arguments("${ARGN}"
        PKG_NAME PKG_GIT_REPOSITORY PKG_GIT_TAG 
        PKG_OPTIONS PKG_PACKAGES PKG_VERSION_RANGE PKG_PRIVATE)

    bpm_add_package_to_registry(
        "${PKG_NAME}" "${PKG_GIT_REPOSITORY}" "${PKG_GIT_TAG}" "${PKG_OPTIONS}" 
        "${PKG_PACKAGES}" "${PKG_VERSION_RANGE}" "${PKG_PRIVATE}" "ADD_SUBDIR")

endfunction()

function(bpm_get_cache_dir RESULT_VAR)
    set(cache_dir "")

    if(DEFINED BPM_CACHE AND NOT "${BPM_CACHE}" STREQUAL "")
        set(cache_dir "${BPM_CACHE}")
        message(STATUS "BPM [${PROJECT_NAME}]: resolve BPM_CACHE - from CMAKE_ARG: ${BPM_CACHE}")
    elseif(DEFINED ENV{BPM_CACHE} AND NOT "$ENV{BPM_CACHE}" STREQUAL "")
        set(cache_dir "$ENV{BPM_CACHE}")
        message(STATUS "BPM [${PROJECT_NAME}]: resolve BPM_CACHE - from environment variable: ${cache_dir}")
    else()
        set(cache_dir "${CMAKE_BINARY_DIR}/_deps")
        message(STATUS "BPM [${PROJECT_NAME}]: resolve BPM_CACHE - no cache provided: use local: ${cache_dir}")
    endif()

    # turn into absolute path
    cmake_path(ABSOLUTE_PATH cache_dir BASE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}" NORMALIZE)

    set(${RESULT_VAR} "${cache_dir}" PARENT_SCOPE)
endfunction()

function(bpm_tag_cache_covers_range_heuristic IN_VERSIONS RANGE OUT)
    list(LENGTH RANGE LIST_SIZE)

    # check if it is a range or a singular version (like for named git tags or hashes)
    if(LIST_SIZE EQUAL 1)
        # it is a named tag (not a version tag) or a commit hash
        list(GET RANGE 0 tag_or_hash)
        message(FATAL_ERROR "Internal error. Range coverage for git tags and commit hashes currently unsupported. tag or hash: '${tag_or_hash}'. versions: '${IN_VERSIONS}'")
        return()
    endif()

    list(GET RANGE 0 version_lower)
    list(GET RANGE 1 version_upper)

    if("${version_upper}" STREQUAL "inf.inf.inf")
        set(${OUT} FALSE PARENT_SCOPE)
        return()
    endif()

    # check if we are searching for an exact match
    set(exact_match FALSE)
    string(REGEX MATCH "([0-9]+)\\.([0-9]+)\\.([0-9]+)" _ "${version_lower}")
    if(CMAKE_MATCH_0)
        set(major_lower ${CMAKE_MATCH_1})
        set(minor_lower ${CMAKE_MATCH_2})
        set(patch_lower ${CMAKE_MATCH_3})
    else()
        message(FATAL_ERROR "BPM [${PROJECT_NAME}]: Internal error: version string is not correct: ${version_lower}")
    endif()

    math(EXPR patch_lower_plus_one "${patch_lower} + 1")
    set(version_lower_plus_one "${major_lower}.${minor_lower}.${patch_lower_plus_one}")

    if(version_lower_plus_one VERSION_EQUAL version_upper)
        # find exact match
        foreach(tag IN LISTS IN_VERSIONS)
            string(REGEX MATCH "([0-9]+\\.[0-9]+\\.[0-9]+)" _ "${tag}")
            if(CMAKE_MATCH_1)
                set(vtag ${CMAKE_MATCH_1})
                if("${version_lower}" VERSION_EQUAL "${vtag}")
                    # found exact match
                    set(${OUT} TRUE PARENT_SCOPE)
                    return()
                endif()
            endif()
        endforeach()

        # did not find exact match
        set(${OUT} FALSE PARENT_SCOPE)
        return()
    else()
        # find range
        set(has_larger FALSE)
        set(contains_at_least_one_larger_than_range FALSE)
        set(contains_at_least_one_in_range FALSE)
        

        if("${version_upper}" STREQUAL "inf.inf.inf")
            set(has_larger FALSE)
        else()
            foreach(tag IN LISTS IN_VERSIONS)
                string(REGEX MATCH "([0-9]+\\.[0-9]+\\.[0-9]+)" _ "${tag}")
                if(CMAKE_MATCH_1)
                    set(vtag ${CMAKE_MATCH_1})
                    if("${version_upper}" VERSION_LESS "${vtag}")
                        set(contains_at_least_one_larger_than_range TRUE)
                    elseif(("${version_lower}" VERSION_LESS_EQUAL "${vtag}") AND ("${vtag}" VERSION_LESS "${version_upper}"))
                        set(contains_at_least_one_in_range TRUE)
                    endif()

                    if(contains_at_least_one_larger_than_range AND contains_at_least_one_in_range)
                        set(${OUT} TRUE PARENT_SCOPE)
                        return()
                    endif()
                endif()
            endforeach()
        endif()
        set(${OUT} FALSE PARENT_SCOPE)
    endif()
endfunction()

#
# @param IN_TAGS a sorted list of tags (newest first)
#
function(bpm_filter_version_tags IN_TAGS IN_RANGE OUT_FILTERED_TAGS)

    list(GET IN_RANGE 0 version_lower)
    list(GET IN_RANGE 1 version_upper)

    # filter tags that match the version range
    set(filtered_version_tags "")

    foreach(tag IN LISTS IN_TAGS)
        string(REGEX MATCH "([0-9]+\\.[0-9]+\\.[0-9]+)" _ "${tag}")
        if(CMAKE_MATCH_1)
            set(vtag ${CMAKE_MATCH_1})
            if("${version_lower}" VERSION_LESS_EQUAL "${vtag}")
                if("${version_upper}" STREQUAL "inf.inf.inf")
                    LIST(APPEND filtered_version_tags "${tag}")
                elseif("${vtag}" VERSION_LESS "${version_upper}")
                    LIST(APPEND filtered_version_tags "${tag}")
                endif()
            endif()
        endif()
    endforeach()

    set(${OUT_FILTERED_TAGS} "${filtered_version_tags}" PARENT_SCOPE)
endfunction()

function(bpm_is_version_in_range in_pkg_name in_mirror in_mirror_lock_file in_tag in_range out)
    list(LENGTH in_range range_size)
    if(range_size EQUAL 1)
        # range is only 1 element.
        # this means we are looking for an exact match
        # also semver compliant git tags will always have a range of 2.
        # so in_range is not version
        # compare git hashes for equality

        file(LOCK "${in_mirror_lock_file}")
            execute_process(COMMAND git --git-dir "${in_mirror}" rev-parse "${in_tag}^{commit}" RESULT_VARIABLE res OUTPUT_VARIABLE PKG_GIT_COMMIT_A OUTPUT_STRIP_TRAILING_WHITESPACE)
            if(NOT res EQUAL 0)
                message(FATAL_ERROR "BPM [${PROJECT_NAME}:${in_pkg_name}]: Cannot convert tag: ${in_tag} to commit-hash")
            endif()
            
            execute_process(COMMAND git --git-dir "${in_mirror}" rev-parse "${in_range}^{commit}" RESULT_VARIABLE res OUTPUT_VARIABLE PKG_GIT_COMMIT_B OUTPUT_STRIP_TRAILING_WHITESPACE)
            if(NOT res EQUAL 0)
                # try again with leading v
                execute_process(COMMAND git --git-dir "${in_mirror}" rev-parse "v${in_range}^{commit}" RESULT_VARIABLE res OUTPUT_VARIABLE PKG_GIT_COMMIT_B OUTPUT_STRIP_TRAILING_WHITESPACE)
                if(NOT res EQUAL 0)
                    message(FATAL_ERROR "BPM [${PROJECT_NAME}:${in_pkg_name}]: Cannot convert tag: '${in_range}' to commit-hash")
                endif()
            endif()
        file(LOCK "${in_mirror_lock_file}" RELEASE)

        if("${PKG_GIT_COMMIT_A}" STREQUAL "${PKG_GIT_COMMIT_B}")
            set(${out} TRUE PARENT_SCOPE)   
            return()
        else()
            set(${out} FALSE PARENT_SCOPE)
            return()
        endif()
    endif()

    if(NOT range_size EQUAL 2)
        message(FATAL_ERROR "'in_range' (${in_range}) has wrong size. should be 1 or 2")
    endif()

    list(GET in_range 0 range_lower)
    list(GET in_range 1 range_upper)

    string(REGEX MATCH "([0-9]+\\.[0-9]+\\.[0-9]+)" _ "${in_tag}")
    if(CMAKE_MATCH_1)
        set(in_version ${CMAKE_MATCH_1})
        
        if("${range_lower}" VERSION_LESS_EQUAL "${in_version}")
            if("${range_upper}" STREQUAL "inf.inf.inf")
                set(${out} TRUE PARENT_SCOPE)
                return()
            elseif("${in_version}" VERSION_LESS "${range_upper}")
                set(${out} TRUE PARENT_SCOPE)
                return()
            else()
                set(${out} FALSE PARENT_SCOPE)
                return()
            endif()
        else()
            set(${out} FALSE PARENT_SCOPE)
            return()
        endif()
    else()
        file(LOCK "${in_mirror_lock_file}")
            execute_process(COMMAND git --git-dir "${in_mirror}" rev-parse "${in_tag}^{commit}" RESULT_VARIABLE res OUTPUT_VARIABLE tag_commit OUTPUT_STRIP_TRAILING_WHITESPACE)
            if(NOT res EQUAL 0)
                message(FATAL_ERROR "BPM [${PROJECT_NAME}:${in_pkg_name}]: Cannot convert tag: '${in_tag}' to commit-hash")
            endif()
            
            execute_process(COMMAND git --git-dir "${in_mirror}" rev-parse "${range_lower}^{commit}" RESULT_VARIABLE res OUTPUT_VARIABLE range_low_commit OUTPUT_STRIP_TRAILING_WHITESPACE)
            if(NOT res EQUAL 0)
                # try again with leading v
                execute_process(COMMAND git --git-dir "${in_mirror}" rev-parse "'v${range_lower}'^{commit}" RESULT_VARIABLE res OUTPUT_VARIABLE range_low_commit OUTPUT_STRIP_TRAILING_WHITESPACE)
                if(NOT res EQUAL 0)
                    message(FATAL_ERROR "BPM [${PROJECT_NAME}:${in_pkg_name}]: Cannot convert tag: '${range_lower}' to commit-hash")
                endif()
            endif()   
        file(LOCK "${in_mirror_lock_file}" RELEASE)

        execute_process(COMMAND git --git-dir "${in_mirror}" merge-base --is-ancestor "${range_low_commit}" "${tag_commit}" RESULT_VARIABLE res)
        
        if(res EQUAL 0)
            # yes commit A is older than commit B

            if("${range_upper}" STREQUAL "inf.inf.inf")
                set(${out} TRUE PARENT_SCOPE)
                return()
            else()
                # convert upper range to commit
                execute_process(COMMAND git --git-dir "${in_mirror}" rev-parse "${range_upper}^{commit}" RESULT_VARIABLE res OUTPUT_VARIABLE range_high_commit OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
                if(NOT res EQUAL 0)
                    # try again with leading v
                    execute_process(COMMAND git --git-dir "${in_mirror}" rev-parse "v${range_upper}^{commit}" RESULT_VARIABLE res OUTPUT_VARIABLE range_low_commit OUTPUT_STRIP_TRAILING_WHITESPACE)
                    if(NOT res EQUAL 0)
                        # check if the tag is semver compliant
                        string(REGEX MATCH "([0-9]+\\.[0-9]+\\.[0-9]+)" _ "${range_upper}")
                        if(CMAKE_MATCH_0)
                            # tag is semver compliant - assume upper bound is a made up tag that does not exist yet --> treat like inf
                            set(${out} TRUE PARENT_SCOPE)
                            return()
                        else()
                            message(FATAL_ERROR "BPM [${PROJECT_NAME}:${in_pkg_name}]: Cannot convert tag: '${range_upper}' to commit-hash")
                        endif()
                    endif()
                endif()

                if("${tag_commit}" STREQUAL "${range_high_commit}")
                    # range high is an open range border --> so if equal --> false
                    set(${out} FALSE PARENT_SCOPE)
                    return()
                else()
                    execute_process(COMMAND git --git-dir "${in_mirror}" merge-base --is-ancestor "${tag_commit}" "${range_high_commit}" RESULT_VARIABLE res)
                    if(res EQUAL 0)
                        # success --> tag is in the range
                        set(${out} TRUE PARENT_SCOPE)
                    else()
                        # failure --> tag is not in the range
                        set(${out} FALSE PARENT_SCOPE)
                    endif()
                endif()
        
            endif()
        else()
            # not commit A is older than commit B
            set(${out} FALSE PARENT_SCOPE)
        endif()   
    endif()
    
endfunction()


 
function(bpm_load_tag_list mirror_dir mirror_lock_file out_tags)
    file(LOCK "${mirror_lock_file}")
        # --sort=-committerdate as a pre sort
        execute_process(COMMAND git --git-dir "${mirror_dir}" tag --sort=version:refname RESULT_VARIABLE res OUTPUT_VARIABLE tags ERROR_QUIET)
    file(LOCK "${mirror_lock_file}" RELEASE)

    if(NOT res EQUAL 0)
        message(FATAL_ERROR "BPM [${PROJECT_NAME}:${PKG_NAME}]: Failed to get tags from mirror: ${mirror_dir}." )
    endif()

    # turn the console output into a CMake list
    string(REPLACE "\r\n" "\n" tags "${tags}")
    string(REPLACE "\n" ";" tags "${tags}")
    set(${out_tags} ${tags} PARENT_SCOPE)
endfunction()

function(bpm_solve_dependencies BPM_CACHE_DIR in_packages out_selected_list)

    set(logging)
    set(logging_conflicts_only)

    set(decision_counter "0")
    set(_todo_list)
    foreach(ipkg IN LISTS in_packages)
        string(STRIP "${ipkg}" ipkg)
        if(NOT ipkg STREQUAL "")
            list(APPEND _todo_list "${ipkg} REQUIRED_FROM ${PROJECT_NAME}")
        endif()
    endforeach()
    set(decision_${decision_counter}_todo_list "${_todo_list}")
    set(decision_${decision_counter}_selected_list "")
    set(decision_${decision_counter}_pkg_name "")
    set(decision_${decision_counter}_version "")
    set(decision_${decision_counter}_range "")
    set(decision_${decision_counter}_tag_wheel "")
    set(decision_${decision_counter}_git_tag "")
    set(decision_${decision_counter}_git_repo "")
    set(decision_${decision_counter}_type "")
    set(decision_${decision_counter}_required_from "")
    set(decision_${decision_counter}_options "")
    set(decision_${decision_counter}_packages "")

    # print current decision
    
    set(short_to_do "")
    foreach(todo IN LISTS decision_${decision_counter}_todo_list)
        # clear to avoid remaining state
        set(TODO_NAME)

        # parse the todo entry
        set(options)
        set(oneValueArgs NAME)
        set(multiValueArgs)
        separate_arguments(todo_tokens UNIX_COMMAND "${todo}")
        cmake_parse_arguments(TODO "${options}" "${oneValueArgs}" "${multiValueArgs}" ${todo_tokens})
        list(APPEND short_to_do "${TODO_NAME}")
        string(REPLACE ";" "," short_to_do "${short_to_do}")
    endforeach()
    string(REPLACE ";" "-" msg_version_range "${decision_${decision_counter}_range}")
    set(msg "decision: ${decision_counter} | todo: [${short_to_do}], pkg: ${decision_${decision_counter}_pkg_name}, version: ${decision_${decision_counter}_version}, range: ${msg_version_range}, wheel: ${decision_${decision_counter}_tag_wheel}, type: ${decision_${decision_counter}_type}, required from: ${decision_${decision_counter}_required_from}")
    if(BPM_VERBOSE)
        message(STATUS "BPM [${PROJECT_NAME}]: ${msg}")
    endif()
    set(logging "${msg}")
    
    

    set(solved_one FALSE)
    set(version_conflict FALSE)

    while(decision_${decision_counter}_todo_list)
    
        set(todo_list ${decision_${decision_counter}_todo_list})
        list(POP_FRONT todo_list pkg)

        # clear vars to avoid remaining state
        set(PKG_NAME)
        set(PKG_VERSION_RANGE)
        set(PKG_GIT_TAG)
        set(PKG_GIT_REPOSITORY)
        set(PKG_TYPE)
        set(PKG_REQUIRED_FROM)
        set(PKG_OPTIONS)
        set(PKG_PACKAGES)

        set(options)
        set(oneValueArgs NAME VERSION_RANGE GIT_TAG GIT_REPOSITORY TYPE REQUIRED_FROM)
        set(multiValueArgs OPTIONS PACKAGES DEPENDENCIES)
        separate_arguments(pkg_tokens UNIX_COMMAND "${pkg}")
        cmake_parse_arguments(PKG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${pkg_tokens})
        string(REPLACE "-" ";" PKG_VERSION_RANGE ${PKG_VERSION_RANGE})

        set(mirror_lock_file "${BPM_CACHE_DIR}/${PKG_NAME}/mirror.lock")
        set(mirror_dir "${BPM_CACHE_DIR}/${PKG_NAME}/mirror")

        if(version_conflict) # if there was a version conflict
            if(decision_${decision_counter}_tag_wheel)
                # get next version from tag wheel
                list(POP_BACK decision_${decision_counter}_tag_wheel top_version)

                set(decision_${decision_counter}_version ${top_version})
                
                set(options_entry)
                if(decision_${decision_counter}_options)
                    string(REPLACE ";" " " options_ "${decision_${decision_counter}_options}")
                    set(options_entry "OPTIONS ${options_}")
                endif()
                set(packages_entry)
                if(decision_${decision_counter}_packages)
                    string(REPLACE ";" " " packages_ "${decision_${decision_counter}_packages}")
                    set(packages_entry "PACKAGES ${packages_}")
                endif()
                set(entry "NAME ${decision_${decision_counter}_pkg_name} VERSION ${top_version} GIT_TAG ${decision_${decision_counter}_git_tag} GIT_REPO ${decision_${decision_counter}_git_repo} TYPE ${decision_${decision_counter}_type} REQUIRED_FROM ${PKG_REQUIRED_FROM} ${options_entry} ${packages_entry}")
                math(EXPR prev_decision_counter "${decision_counter} - 1")
                set(decision_${decision_counter}_selected_list "${decision_${prev_decision_counter}_selected_list}")
                LIST(APPEND decision_${decision_counter}_selected_list "${entry}")

                # print current decision
                set(short_to_do "")
                foreach(todo IN LISTS decision_${decision_counter}_todo_list)
                    # clear to avoid remaining state
                    set(TODO_NAME)

                    set(options)
                    set(oneValueArgs NAME)
                    set(multiValueArgs)
                    separate_arguments(todo_tokens UNIX_COMMAND "${todo}")
                    cmake_parse_arguments(TODO "${options}" "${oneValueArgs}" "${multiValueArgs}" ${todo_tokens})
                    list(APPEND short_to_do "${TODO_NAME}")
                    string(REPLACE ";" "," short_to_do "${short_to_do}")
                endforeach()
                string(REPLACE ";" "-" msg_version_range "${decision_${decision_counter}_range}")
                set(msg "update decision: ${decision_counter} | todo: [${short_to_do}], pkg: ${decision_${decision_counter}_pkg_name}, version: ${decision_${decision_counter}_version}, range: ${msg_version_range}, wheel: ${decision_${decision_counter}_tag_wheel}, TYPE: ${decision_${decision_counter}_type}, OPTIONS: ${decision_${decision_counter}_options}, PACKAGES ${decision_${decision_counter}_packages}, required from: ${decision_${decision_counter}_required_from}")
                
                if(BPM_VERBOSE)
                    message(STATUS "BPM [${PROJECT_NAME}]: ${msg}")
                endif()
                string(APPEND logging "\n${msg}")

                # retry logic with lower version
                set(version_conflict FALSE)
                continue()
                
            else()
                # no more versions to try
                # delete this entry 
                set(msg "    version conflict: pop decision: ${decision_counter}")
                if(BPM_VERBOSE)
                    message(STATUS "BPM [${PROJECT_NAME}]: ${msg}")
                endif()
                string(APPEND logging "\n${msg}")

                unset(decision_${decision_counter}_todo_list)
                unset(decision_${decision_counter}_selected_list)
                unset(decision_${decision_counter}_pkg_name)
                unset(decision_${decision_counter}_version)
                unset(decision_${decision_counter}_range)
                unset(decision_${decision_counter}_tag_wheel)
                unset(decision_${decision_counter}_git_tag)
                unset(decision_${decision_counter}_git_repo)
                unset(decision_${decision_counter}_type)
                unset(decision_${decision_counter}_required_from)
                unset(decision_${decision_counter}_options)
                unset(decision_${decision_counter}_packages)

                # pop
                math(EXPR decision_counter "${decision_counter} - 1")
                if(decision_counter LESS 0)
                    set(filename "${CMAKE_BINARY_DIR}/bpm-dependency-solver-log.txt")
                    file(WRITE "${filename}" "${logging}")
                    list(REMOVE_DUPLICATES logging_conflicts_only)
                    string(REPLACE ";" "\n" logging_conflicts_only "${logging_conflicts_only}")
                    message(FATAL_ERROR "Dependency resolution failed:\n----\n${logging_conflicts_only}\n----\nSee '${filename}' for details")
                endif()
                continue()
            endif()
        endif()

        set(solved_one FALSE)
        set(version_conflict FALSE)

        # check if the package is already in the selected list
        set(selected_count 0)
        foreach(selected IN LISTS decision_${decision_counter}_selected_list)
            #clear to avoid remaining state
            set(SEL_NAME) 
            set(SEL_VERSION)
            set(SEL_GIT_TAG)
            set(SEL_GIT_REPO)
            set(SEL_TYPE)
            set(SEL_OPTIONS)
            set(SEL_PACKAGES)
            set(SEL_REQUIRED_FROM)
            set(SEL_DEPENDENCIES)

            set(options)
            set(oneValueArgs NAME VERSION GIT_TAG GIT_REPO TYPE )
            set(multiValueArgs OPTIONS REQUIRED_FROM PACKAGES DEPENDENCIES)
            separate_arguments(selected_tokens UNIX_COMMAND "${selected}")
            cmake_parse_arguments(SEL "${options}" "${oneValueArgs}" "${multiValueArgs}" ${selected_tokens})

            if(("${SEL_NAME}" STREQUAL "${PKG_NAME}") OR ("${PKG_GIT_REPOSITORY}" STREQUAL "${SEL_GIT_REPO}"))
                # repo conflict?
                path_equal("${PKG_GIT_REPOSITORY}" "${SEL_GIT_REPO}" are_repos_equal)

                if(NOT are_repos_equal)
                    set(filename "${CMAKE_BINARY_DIR}/bpm-dependency-solver-log.txt")
                    file(WRITE "${filename}" "${logging}")
                    message(FATAL_ERROR  "BPM [${PROJECT_NAME}:${PKG_NAME}]: Repository conflict: '${PKG_GIT_REPOSITORY}' required from '${PKG_REQUIRED_FROM}' <-- vs --> '${SEL_GIT_REPO}' required from '${SEL_REQUIRED_FROM}'. See '${filename}'")
                endif()

                # package name conflict ? 
                if(NOT "${PKG_NAME}" STREQUAL "${SEL_NAME}")
                    set(filename "${CMAKE_BINARY_DIR}/bpm-dependency-solver-log.txt")
                    file(WRITE "${filename}" "${logging}")
                    message(FATAL_ERROR  "BPM [${PROJECT_NAME}:${PKG_GIT_REPOSITORY}]: Package name conflict:  '${PKG_NAME}' required from '${PKG_REQUIRED_FROM}' <-- vs -->  '${SEL_NAME}' required from '${SEL_REQUIRED_FROM}'. See '${filename}'")
                endif()

                # type conflict?
                if(NOT "${PKG_TYPE}" STREQUAL "${SEL_TYPE}")
                    set(filename "${CMAKE_BINARY_DIR}/bpm-dependency-solver-log.txt")
                    file(WRITE "${filename}" "${logging}")
                    message(FATAL_ERROR "BPM [${PROJECT_NAME}:${PKG_NAME}]: Package type conflict: '${PKG_TYPE}' required from '${PKG_REQUIRED_FROM}' <-- vs --> '${SEL_TYPE}' required from '${SEL_REQUIRED_FROM}'. See '${filename}'")
                endif()

                # resolve option conflicts
                # combine options, remove duplicates, check if there are conflicting ones
                bpm_combine_options("${PKG_OPTIONS}" "${SEL_OPTIONS}" "${PKG_REQUIRED_FROM}" "${SEL_REQUIRED_FROM}" combined_options added_options)

                # combine packages
                set(combined_packages ${PKG_PACKAGES})
                list(APPEND combined_packages ${SEL_PACKAGES})
                list(REMOVE_DUPLICATES combined_packages)

                # check if the selected version is in the constraint
                bpm_is_version_in_range("${PKG_NAME}" "${mirror_dir}" "${mirror_lock_file}" "${SEL_VERSION}" "${PKG_VERSION_RANGE}" is_in_range)
                
                if(is_in_range)
                    # no version conflict - already part of the list
                    # record decision

                    math(EXPR next_decision_counter "${decision_counter} + 1")
                    set(decision_${next_decision_counter}_todo_list "${todo_list}")
                    
                    # update required from:
                    set(updated_required_from "${SEL_REQUIRED_FROM}")
                    list(APPEND updated_required_from "${PKG_REQUIRED_FROM}")
                    list(REMOVE_DUPLICATES updated_required_from)
                    string(REPLACE ";" " " updated_required_from "${updated_required_from}")

                    # update selected list: options or packages might have changed
                        set(updated_selected_list "${decision_${decision_counter}_selected_list}")
                        list(REMOVE_AT updated_selected_list ${selected_count})

                        # turn option list to string

                        string(REPLACE ";" " " combined_options_ "${combined_options}")
                        set(pkg_options_entry "OPTIONS ${combined_options_}")
                        set(pkg_packages_entry)
                        if(combined_packages)
                            # turn list to multi arg list
                            string(REPLACE ";" " " combined_packages_ "${combined_packages}")
                            set(pkg_packages_entry "PACKAGES ${combined_packages_}")
                        endif()

                        set(pkg_dependency_entry)
                        if(SEL_DEPENDENCIES)
                            string(REPLACE ";" " " sel_dependencies_ "${SEL_DEPENDENCIES}")
                            set(pkg_dependency_entry "DEPENDENCIES ${sel_dependencies_}")
                        endif()
                        set(entry "NAME ${PKG_NAME} VERSION ${SEL_VERSION} GIT_TAG ${PKG_GIT_TAG} GIT_REPO ${PKG_GIT_REPOSITORY} TYPE ${PKG_TYPE} REQUIRED_FROM ${updated_required_from} ${pkg_options_entry} ${pkg_packages_entry} ${pkg_dependency_entry}")
                        list(APPEND updated_selected_list "${entry}")
                        set(decision_${next_decision_counter}_selected_list "${updated_selected_list}")
                    

                    set(decision_${next_decision_counter}_pkg_name "${PKG_NAME}")
                    set(decision_${next_decision_counter}_version "${SEL_VERSION}")
                    set(decision_${next_decision_counter}_tag_wheel "${decision_${decision_counter}_tag_wheel}")
                    set(decision_${next_decision_counter}_range "${PKG_VERSION_RANGE}")
                    set(decision_${next_decision_counter}_git_tag "${PKG_GIT_TAG}")
                    set(decision_${next_decision_counter}_git_repo "${PKG_GIT_REPOSITORY}")
                    set(decision_${next_decision_counter}_type "${PKG_TYPE}")
                    set(decision_${next_decision_counter}_required_from "${PKG_REQUIRED_FROM}")
                    set(decision_${next_decision_counter}_options "${PKG_OPTIONS}")
                    set(decision_${next_decision_counter}_packages "${PKG_PACKAGES}")

                    # increment
                    set(decision_counter "${next_decision_counter}")

                    # print current decision                    
                    set(short_to_do "")
                    foreach(todo IN LISTS decision_${decision_counter}_todo_list)
                        set(options)
                        set(oneValueArgs NAME)
                        set(multiValueArgs)
                        separate_arguments(todo_tokens UNIX_COMMAND "${todo}")
                        cmake_parse_arguments(TODO "${options}" "${oneValueArgs}" "${multiValueArgs}" ${todo_tokens})
                        list(APPEND short_to_do "${TODO_NAME}")
                        string(REPLACE ";" "," short_to_do "${short_to_do}")
                    endforeach()
                    string(REPLACE ";" "-" msg_version_range "${decision_${decision_counter}_range}")
                    set(msg "decision: ${decision_counter} | todo: [${short_to_do}], pkg: ${decision_${decision_counter}_pkg_name}, version: ${decision_${decision_counter}_version}, range: ${msg_version_range}, wheel: ${decision_${decision_counter}_tag_wheel}, TYPE: ${decision_${decision_counter}_type}, OPTIONS: ${decision_${decision_counter}_options}, PACKAGES ${decision_${decision_counter}_packages}, required from: ${decision_${decision_counter}_required_from}")
                    if(BPM_VERBOSE)
                        message(STATUS "BPM [${PROJECT_NAME}]: ${msg}")
                    endif()
                    string(APPEND logging "\n${msg}")

                    set(solved_one TRUE)
                    break() # to leave for loop and then continue in the while loop
                else()
                    # found a version conflict
                    string(REPLACE ";" "-" msg_version_range "${PKG_VERSION_RANGE}")

                    # find the version of the required package
                    set(pkg_req_version)
                    set(sel_req_version)
                    foreach(_sel IN LISTS decision_${decision_counter}_selected_list)
                        set(REQ_NAME)
                        set(REQ_VERSION)

                        set(options)
                        set(oneValueArgs NAME VERSION)
                        set(multiValueArgs)
                        separate_arguments(_sel_tokens UNIX_COMMAND "${_sel}")
                        cmake_parse_arguments(REQ "${options}" "${oneValueArgs}" "${multiValueArgs}" ${_sel_tokens})
                        
                        if("${PKG_REQUIRED_FROM}" STREQUAL "${REQ_NAME}")
                            set(pkg_req_version "${REQ_VERSION}")
                        endif()
                        if("${SEL_REQUIRED_FROM}" STREQUAL "${REQ_NAME}")
                            set(sel_req_version "${REQ_VERSION}")
                        endif()

                        if(pkg_req_version AND sel_req_version)
                            break()
                        endif()
                    endforeach()
                    # print error message
                    if(sel_req_version)
                        set(msg "version conflict: '${PKG_NAME}': ${msg_version_range} required from '${PKG_REQUIRED_FROM}#${pkg_req_version}', existing version: ${SEL_VERSION} required from '${SEL_REQUIRED_FROM}#${sel_req_version}'")
                    else()
                        set(msg "version conflict: '${PKG_NAME}': ${msg_version_range} required from '${PKG_REQUIRED_FROM}#${pkg_req_version}', existing version: ${SEL_VERSION} required from '${SEL_REQUIRED_FROM}'")
                    endif()

                    if(BPM_VERBOSE)
                        message(STATUS "BPM [${PROJECT_NAME}]:   ${msg}")
                    endif()

                    string(APPEND logging "\n  ${msg}")
                    
                    list(APPEND logging_conflicts_only "${msg}")

                    set(version_conflict TRUE)
                    break() # to leave for loop and then continue in the while loop
                endif()
            endif()

            math(EXPR selected_count "${selected_count} + 1")
        endforeach()

        # continue whith the next package if one decision got solved in the for loop
        if(solved_one)
            continue()
        endif()

        if(version_conflict)
            continue()
        endif()

        string(REPLACE ";" "_" temp_version_range_str "${PKG_VERSION_RANGE}")
        string(REPLACE "-" "_" temp_version_range_str "${temp_version_range_str}")
        string(REPLACE "." "_" temp_version_range_str "${temp_version_range_str}")
        if(NOT DEFINED cached_tag_wheel_${PKG_NAME}_${temp_version_range_str})
            # all the following can be skipped if the tag wheel was already cached

            # package is not yet in the selected list
            # build the version wheel for this package
            if(NOT mirror_${PKG_NAME}_up_to_date)
                if(NOT EXISTS "${mirror_dir}/HEAD")
                    set(relative_mirror_dir "\${BPM_CACHE}/${PKG_NAME}/mirror")
                    if(BPM_VERBOSE)
                        message(STATUS "BPM [${PROJECT_NAME}:${PKG_NAME}]: Cloning git repository: '${PKG_GIT_REPOSITORY}' into '${relative_mirror_dir}'")
                    else()
                        message(STATUS "BPM [${PROJECT_NAME}:${PKG_NAME}]: Cloning '${PKG_NAME}' into '${relative_mirror_dir}'")
                    endif()
                    if(BPM_NO_DOWNLOAD)
                        message(FATAL_ERROR "BPM [${PROJECT_NAME}:${PKG_NAME}]: Mirror does not exist: '${relative_mirror_dir}'. Cannot download from '${PKG_GIT_REPOSITORY}', because `NO_DOWNLOAD` was provided.")
                    endif()
                    file(LOCK "${mirror_lock_file}")
                        if(NOT EXISTS "${mirror_dir}/HEAD")
                            if(BPM_VERBOSE)
                                execute_process(COMMAND git clone --mirror "${PKG_GIT_REPOSITORY}" "${mirror_dir}" --recursive -c advice.detachedHead=false RESULT_VARIABLE res)
                            else()
                                execute_process(COMMAND git clone --mirror "${PKG_GIT_REPOSITORY}" "${mirror_dir}" --recursive -c advice.detachedHead=false RESULT_VARIABLE res OUTPUT_QUIET ERROR_QUIET)
                            endif()
                        endif()
                    file(LOCK "${mirror_lock_file}" RELEASE)

                    if(res EQUAL 0)
                        if(BPM_VERBOSE)
                            message(STATUS "BPM [${PROJECT_NAME}:${PKG_NAME}]: Cloning git repository - success")
                        endif()
                        set(mirror_${PKG_NAME}_up_to_date TRUE)
                    else() 
                        message(FATAL_ERROR "BPM [${PROJECT_NAME}:${PKG_NAME}]: Cloning git repository - failed")
                    endif()
                endif()
            endif()

            

            # see if tags have already been acquired
            if(DEFINED BPM_REGISTRY_${PKG_NAME}_GIT_TAGS)
                set(tags "${BPM_REGISTRY_${PKG_NAME}_GIT_TAGS}")
            else()
                # if tags have not been acquired - load them
                bpm_load_tag_list("${mirror_dir}" "${mirror_lock_file}" tags)
                #if(NOT tags)
                #    set(filename "${CMAKE_BINARY_DIR}/bpm-dependency-solver-log.txt")
                #    file(WRITE "${filename}" "${logging}")
                #    message(FATAL_ERROR "BPM [${PKG_NAME}]: Has no tags to checkout. Mirror: ${mirror_dir}. See '${filename}' for details.")
                #endif()
                # cache tags
                set("BPM_REGISTRY_${PKG_NAME}_GIT_TAGS" "${tags}")
            endif()

            if(NOT mirror_${PKG_NAME}_up_to_date)
                # check if the version range is fully contained with the mirrors tags
                list(LENGTH PKG_VERSION_RANGE range_size)
                if(range_size EQUAL 1)
                    # it is a named tag (not a version tag) or a commit hash
                    set(tag_or_hash "${PKG_VERSION_RANGE}")
                    file(LOCK "${mirror_lock_file}")
                        execute_process(COMMAND git --git-dir "${mirror_dir}" cat-file -e "${tag_or_hash}^{commit}" RESULT_VARIABLE res )
                    file(LOCK "${mirror_lock_file}" RELEASE)
                                        
                    if(res EQUAL 0)
                        set(contains TRUE)
                    else()
                        set(contains FALSE)
                    endif()
                else()
                    # is a regular version tag that spans a version range
                    bpm_tag_cache_covers_range_heuristic("${tags}" "${PKG_VERSION_RANGE}" contains)
                endif()


                # fetch if upper version bound is inf or tag is not contained
                if(NOT contains)
                    if(BPM_NO_DOWNLOAD)
                        if(NOT fetch_skipped_due_to_no_downloads)
                            message(STATUS "BPM [${PROJECT_NAME}]: mirrors might be out-of-date - fetch - skipped due to `NO_DOWNLOAD`")
                        endif()
                        set(fetch_skipped_due_to_no_downloads TRUE) # to avoid multiple warnings
                    elseif(BPM_NO_DOWNLOAD_UPDATES)
                        if(NOT fetch_updates_skipped_due_to_no_download_updates)
                            message(STATUS "BPM [${PROJECT_NAME}]: mirrors might be out-of-date - fetch - skipped due to `BPM_NO_DOWNLOAD_UPDATES`")
                        endif()
                        set(fetch_updates_skipped_due_to_no_download_updates TRUE) # to avoid multiple warnings
                    else()
                        message(STATUS "BPM [${PROJECT_NAME}:${PKG_NAME}]: mirror might be out-of-date - fetch for updates")

                        file(LOCK "${mirror_lock_file}")
                            if(BPM_VERBOSE)
                                execute_process(COMMAND git --git-dir "${mirror_dir}" fetch --tags --prune RESULT_VARIABLE res)
                            else()
                                execute_process(COMMAND git --git-dir "${mirror_dir}" fetch --tags --prune RESULT_VARIABLE res OUTPUT_QUIET ERROR_QUIET)
                            endif()
                        file(LOCK "${mirror_lock_file}" RELEASE)

                        if(res EQUAL 0)
                            if(BPM_VERBOSE)
                                message(STATUS "BPM [${PROJECT_NAME}:${PKG_NAME}]: mirror might be out-of-date - fetch for updates - success")
                            endif()
                        else() 
                            message(WARNING "BPM [${PROJECT_NAME}:${PKG_NAME}]: mirror might be out-of-date - fetch for updates - failed")
                        endif()

                        # re-update tags after fetching
                        bpm_load_tag_list("${mirror_dir}" "${mirror_lock_file}" tags)
                        #if(NOT tags)
                        #    set(filename "${CMAKE_BINARY_DIR}/bpm-dependency-solver-log.txt")
                        #    file(WRITE "${filename}" "${logging}")
                        #    message(FATAL_ERROR "BPM [${PROJECT_NAME}:${PKG_NAME}]: Has no tags to checkout. Mirror: ${mirror_dir}. See '${filename}' for details.")
                        #endif()
    
                        # cache tags
                        set("BPM_REGISTRY_${PKG_NAME}_GIT_TAGS" "${tags}")

                        # TODO: check if at least one element is contained --> error if not
                        set(mirror_${PKG_NAME}_up_to_date TRUE)
                    endif()
                endif()
            endif()


            list(LENGTH PKG_VERSION_RANGE range_size)
            if(range_size EQUAL 1)
                # is a single/exact tag or commit hash
                set(tag_wheel ${PKG_VERSION_RANGE})
            else()
                # check if tags have been loaded
                if(NOT tags)
                    set(filename "${CMAKE_BINARY_DIR}/bpm-dependency-solver-log.txt")
                    file(WRITE "${filename}" "${logging}")
                    message(FATAL_ERROR "BPM [${PROJECT_NAME}:${PKG_NAME}]: Failed to load tags for package '${PKG_NAME}' from mirror: ${mirror_dir}. See '${filename}' for details.")
                endif()

                # is an actual range
                bpm_filter_version_tags("${tags}" "${PKG_VERSION_RANGE}" tag_wheel)
                if(NOT tag_wheel)
                    list(GET IN_RANGE 0 version_lower)
                    list(GET IN_RANGE 1 version_upper)
                    set(filename "${CMAKE_BINARY_DIR}/bpm-dependency-solver-log.txt")
                    file(WRITE "${filename}" "${logging}")
                    message(FATAL_ERROR "BPM [${PROJECT_NAME}:${PKG_NAME}]: No tags in range ${version_lower}-${version_upper}. Mirror: ${mirror_dir}. See '${filename}' for details.")
                endif()
            endif()

            set(cached_tag_wheel_${PKG_NAME}_${temp_version_range_str} "${tag_wheel}")
        else()
            set(tag_wheel "${cached_tag_wheel_${PKG_NAME}_${temp_version_range_str}}")
        endif()

        list(POP_BACK tag_wheel top_version)

        # load and cache metadata

        if(NOT DEFINED metadata_${PKG_NAME}_${top_version})
            file(LOCK "${mirror_lock_file}")
                execute_process(COMMAND git --git-dir "${mirror_dir}" cat-file blob "${top_version}:.bpm-registry" RESULT_VARIABLE res OUTPUT_VARIABLE metadata_tmp ERROR_QUIET)
            file(LOCK "${mirror_lock_file}" RELEASE)

            string(REPLACE "\r\n" "\n" metadata_tmp "${metadata_tmp}") # replace new lines windows to unix style
            string(REPLACE "\n" ";" metadata_${PKG_NAME}_${top_version} "${metadata_tmp}") # replace new lines with ; for list seperators
        endif()

        foreach(line IN LISTS metadata_${PKG_NAME}_${top_version})
            if(line)
                string(STRIP "${line}" line)
                if(NOT line STREQUAL "")
                    list(APPEND todo_list "${line} REQUIRED_FROM ${PKG_NAME}")
                endif()
            endif()
        endforeach()

        # no registry found --> no more entries added to the todo list --> make decision and continue
        math(EXPR next_decision_counter "${decision_counter} + 1")
        
        set(pkg_options_entry)
        if(PKG_OPTIONS)
            # turn list to multi arg list
            string(REPLACE ";" " " PKG_OPTIONS_ "${PKG_OPTIONS}")
            set(pkg_options_entry "OPTIONS ${PKG_OPTIONS_}")
        endif()

        set(pkg_packages_entry)
        if(PKG_PACKAGES)
            # turn list to multi arg list
            string(REPLACE ";" " " PKG_PACKAGES_ "${PKG_PACKAGES}")
            set(pkg_packages_entry "PACKAGES ${PKG_PACKAGES_}")
        endif()

        set(pkg_dependencies_entry)
        if(metadata_${PKG_NAME}_${top_version})
            set(pkg_dependencies_entry "DEPENDENCIES")
            foreach(line IN LISTS metadata_${PKG_NAME}_${top_version})
                if(line)
                    string(STRIP "${line}" line)
                    if(NOT line STREQUAL "")
                        # turn list to multi arg list
                        string(REPLACE ";" " " line_ "${line}")

                        # read name of dependency
                        set(options "")
                        set(oneValueArgs NAME)
                        set(multiValueArgs "")
                        separate_arguments(line_tokens UNIX_COMMAND "${line_}")
                        cmake_parse_arguments(LINE "${options}" "${oneValueArgs}" "${multiValueArgs}" ${line_tokens})
                        string(APPEND pkg_dependencies_entry " ${LINE_NAME}")
                    endif()
                endif()
            endforeach()
        endif()
        set(entry "NAME ${PKG_NAME} VERSION ${top_version} GIT_TAG ${PKG_GIT_TAG} GIT_REPO ${PKG_GIT_REPOSITORY} TYPE ${PKG_TYPE} REQUIRED_FROM ${PKG_REQUIRED_FROM} ${pkg_options_entry} ${pkg_packages_entry} ${pkg_dependencies_entry}")

        set(decision_${next_decision_counter}_todo_list "${todo_list}")
        set(decision_${next_decision_counter}_selected_list "${decision_${decision_counter}_selected_list}")
        LIST(APPEND decision_${next_decision_counter}_selected_list "${entry}")
        set(decision_${next_decision_counter}_pkg_name "${PKG_NAME}")
        set(decision_${next_decision_counter}_version "${top_version}")
        set(decision_${next_decision_counter}_tag_wheel "${tag_wheel}")
        set(decision_${next_decision_counter}_range "${PKG_VERSION_RANGE}")
        set(decision_${next_decision_counter}_git_tag "${PKG_GIT_TAG}")
        set(decision_${next_decision_counter}_git_repo "${PKG_GIT_REPOSITORY}")
        set(decision_${next_decision_counter}_type "${PKG_TYPE}")
        set(decision_${next_decision_counter}_required_from "${PKG_REQUIRED_FROM}")
        set(decision_${next_decision_counter}_options "${PKG_OPTIONS}")
        set(decision_${next_decision_counter}_packages "${PKG_PACKAGES}")

        #increment
        set(decision_counter "${next_decision_counter}")

        # print current decision
        set(short_to_do "")
        foreach(todo IN LISTS decision_${decision_counter}_todo_list)
            set(options)
            set(oneValueArgs NAME)
            set(multiValueArgs)
            separate_arguments(todo_tokens UNIX_COMMAND "${todo}")
            cmake_parse_arguments(TODO "${options}" "${oneValueArgs}" "${multiValueArgs}" ${todo_tokens})
            list(APPEND short_to_do "${TODO_NAME}")
            string(REPLACE ";" "," short_to_do "${short_to_do}")
        endforeach()
        string(REPLACE ";" "-" msg_version_range "${decision_${decision_counter}_range}")
        set(msg "new decision: ${decision_counter} | todo: [${short_to_do}], pkg: ${decision_${decision_counter}_pkg_name}, version: ${decision_${decision_counter}_version}, range: ${msg_version_range}, wheel: ${decision_${decision_counter}_tag_wheel}, TYPE: ${decision_${decision_counter}_type}, OPTIONS: ${decision_${decision_counter}_options}, required from: ${decision_${decision_counter}_required_from}")
        if(BPM_VERBOSE)
            message(STATUS "BPM [${PROJECT_NAME}]: ${msg}")
        endif()
        string(APPEND logging "\n${msg}")
        
    endwhile()

    # for each list complete the dependencies recursively
    set("${out_selected_list}" "${decision_${decision_counter}_selected_list}" PARENT_SCOPE)

endfunction()

function(bpm_create_manifest OUT_MANIFEST)
    set(_manifest "")
    foreach(var IN LISTS ARGN)
        if(DEFINED ${var})
            set(_value "${${var}}")
        else()
            set(_value "")
        endif()
        string(APPEND _manifest "${var}=${_value}\n")
    endforeach()
    set(${OUT_MANIFEST} "${_manifest}" PARENT_SCOPE)
endfunction()

function(bpm_try_find_packages lib_name packages lib_install_dir OUT_FOUND_ALL)
    set(all_packages_found TRUE)

    list(LENGTH packages packages_length)
    if(packages_length EQUAL 0)
        message(FATAL_ERROR "BPM [${PROJECT_NAME}:${lib_name}] No packages provided.")
        set(${OUT_FOUND_ALL} FALSE PARENT_SCOPE)
        return()
    endif()

    set(found_packages)
    set(missing_packages)

    foreach(package IN LISTS packages)
        # unset for deterministic find without sideeffects
        unset(${package}_DIR CACHE)
        if(BPM_VERBOSE)
            message(STATUS "BPM [${PROJECT_NAME}:${lib_name}]: Find package: ${package} - looking in: ${lib_install_dir}")
        endif()
        if(BPM_VERBOSE)
            find_package(${package} CONFIG NO_DEFAULT_PATH PATHS "${lib_install_dir}")
        else()
            find_package(${package} CONFIG NO_DEFAULT_PATH PATHS "${lib_install_dir}" QUIET)
        endif()
        unset(${package}_DIR)

        if(${package}_FOUND)

            # export variables to parent scope for later use (e.g. for generating manifest)
            foreach(var IN ITEMS FOUND DIR VERSION CONFIG)
                if(DEFINED ${package}_${var})
                    set(${package}_${var} "${${package}_${var}}" PARENT_SCOPE)
                endif()
            endforeach()

            if(NOT found_packages)
                string(APPEND found_packages "${package}")
            else()
                string(APPEND found_packages ", ${package}")
            endif()
            if(BPM_VERBOSE)
                message(STATUS "BPM [${PROJECT_NAME}:${lib_name}]: Find package: ${package} - found")
            endif()
        else()
            if(NOT missing_packages)
                string(APPEND missing_packages "${package}")
            else()
                string(APPEND missing_packages ", ${package}")
            endif()
            if(BPM_VERBOSE)
                message(STATUS "BPM [${PROJECT_NAME}:${lib_name}]: Find package : ${package} - missing")
            endif()
            set(all_packages_found FALSE)
        endif()
    endforeach()

    if(NOT BPM_VERBOSE)
        if(found_packages)
            message(STATUS "BPM [${PROJECT_NAME}:${lib_name}]: Found packages: ${found_packages}")
        endif()
        if(missing_packages)
            message(STATUS "BPM [${PROJECT_NAME}:${lib_name}]: Missing packages: ${missing_packages}")
        endif()
    endif()

    if(NOT all_packages_found)
        # check which packages are actually in the install dir to give a more detailed error message
        bpm_find_installed_packages("${PKG_NAME}" "${lib_install_dir}" installed_packages)
        if(installed_packages)
            message(STATUS "BPM [${PROJECT_NAME}:${PKG_NAME}]: Installed packages: ${installed_packages}")
        endif()
    endif()

    set(${OUT_FOUND_ALL} ${all_packages_found} PARENT_SCOPE)
endfunction()

#
# @brief Finds all options that contain `test` or `example` (case insensitive) in a file
#
function(bpm_find_test_example_options cmake_file result_var)
    file(READ "${cmake_file}" content)

    set(test_regex "[Tt][Ee][Ss][Tt]")
    set(tests_regex "[Tt][Ee][Ss][Tt][Ss]")
    set(testing_regex "[Tt][Ee][Ss][Tt][Ii][Nn][Gg]")
    set(example_regex "[Ee][Xx][Aa][Mm][Pp][Ll][Ee]")
    set(examples_regex "[Ee][Xx][Aa][Mm][Pp][Ll][Ee][Ss]")
    set(test_or_example_regex "(${test_regex}|${tests_regex}|${testing_regex}|${example_regex})")

    # Find all option(...) or if(...) blocks
    string(REGEX MATCHALL "(option|if)[ \t\r\n]*\\([^)]+\\)" blocks "${content}")

    set(found_options "")

    foreach(block ${blocks})
        # Extract all identifiers from the block (alphanumeric + underscores, starting with letter or underscore)
        string(REGEX MATCHALL "[A-Za-z_][A-Za-z0-9_]*" identifiers "${block}")
        
        foreach(id ${identifiers})
            # Check if the identifier matches the test/example pattern
            if("${id}" MATCHES "(^|[-_ \t\\(\"])${test_or_example_regex}([-_ \t\\)\"]|$)")
                list(APPEND found_options "${id}")
            endif()
        endforeach()
    endforeach()

    list(REMOVE_DUPLICATES found_options)
    set(${result_var} "${found_options}" PARENT_SCOPE)
endfunction()

#
# @brief Finds all options recursively in the provided folder that contain `test` or `example` (case insensitive) in a file
#
function(bpm_find_test_example_options_r search_folder result_var)
    file(GLOB_RECURSE cmake_files "${search_folder}/CMakeLists.txt" "${search_folder}/*.cmake")

    set(all_flags "")

    foreach(file_ ${cmake_files})
        bpm_find_test_example_options(${file_} flags_)
        list(APPEND all_flags ${flags_})
    endforeach()

    list(REMOVE_DUPLICATES all_flags)
    set(${result_var} "${all_flags}" PARENT_SCOPE)
endfunction()

function(bpm_clone_from_mirror lib_name library_mirror_dir lib_src_dir git_tag lib_mirror_lock_file lib_src_lock_file)

    if(BPM_VERBOSE)
        message(STATUS "BPM [${PROJECT_NAME}:${lib_name}]: Cloning mirror '${library_mirror_dir}' into source dir '${lib_src_dir}' at '${git_tag}'")
    endif()

    if(NOT EXISTS ${lib_src_dir}/.git)
        file(LOCK "${lib_src_lock_file}") # lockorder: install -> build -> source -> mirror
            if(NOT EXISTS ${lib_src_dir}/.git)
                file(LOCK "${lib_mirror_lock_file}")
                    # CLONE
                    if(BPM_VERBOSE)
                        execute_process(COMMAND git clone --reference "${library_mirror_dir}" --no-checkout "${library_mirror_dir}" "${lib_src_dir}" -c advice.detachedHead=false RESULT_VARIABLE res)
                    else()
                        execute_process(COMMAND git clone --reference "${library_mirror_dir}" --no-checkout "${library_mirror_dir}" "${lib_src_dir}" -c advice.detachedHead=false RESULT_VARIABLE res OUTPUT_QUIET ERROR_QUIET)
                    endif()
                    if(NOT res EQUAL 0)
                        message(FATAL_ERROR "BPM [${PROJECT_NAME}:${lib_name}]: Failed to clone mirror '${library_mirror_dir}' into source dir '${lib_src_dir}'.")
                    endif()

                    if(BPM_VERBOSE)
                        execute_process(COMMAND git -C "${lib_src_dir}" checkout "${git_tag}" RESULT_VARIABLE res )
                    else()
                        execute_process(COMMAND git -C "${lib_src_dir}" checkout "${git_tag}" RESULT_VARIABLE res OUTPUT_QUIET ERROR_QUIET)
                    endif()
                    if(NOT res EQUAL 0)
                        message(FATAL_ERROR "BPM [${PROJECT_NAME}:${lib_name}]: Failed to checkout '${git_tag}' in source dir '${lib_src_dir}'.")
                    endif()
                
                file(LOCK "${lib_mirror_lock_file}" RELEASE)

                if(res EQUAL 0)
                    if(BPM_VERBOSE)
                        message(STATUS "BPM [${PROJECT_NAME}:${lib_name}]: Cloning mirror into source dir - success")
                    endif()
                else()
                    message(FATAL_ERROR "BPM [${PROJECT_NAME}:${lib_name}]: Cloning mirror '${library_mirror_dir}' into source dir '${lib_src_dir}' - failed")
                endif()
                
                # SUBMODULES
                if(BPM_VERBOSE)
                    message(STATUS "BPM [${PROJECT_NAME}:${lib_name}]: Updating git-submodules")
                endif()

                if(BPM_VERBOSE)
                    execute_process(COMMAND git -C "${lib_src_dir}" submodule update --init --recursive RESULT_VARIABLE res)
                else()
                    execute_process(COMMAND git -C "${lib_src_dir}" submodule update --init --recursive RESULT_VARIABLE res OUTPUT_QUIET)
                endif()
            endif()
        file(LOCK "${lib_src_lock_file}" RELEASE)

        if(res EQUAL 0)
            if(BPM_VERBOSE)
                message(STATUS "Updating git-submodules - success")
            endif()
        else()
            message(FATAL_ERROR "Updating git-submodules - failed")
        endif()
    else()
        if(BPM_VERBOSE)
            message(STATUS "BPM [${PROJECT_NAME}:${lib_name}]: Cloning mirror into source dir - skipped")
        endif()
    endif()

endfunction()

function(bpm_configure_library BPM_CACHE_DIR lib_name lib_src_dir lib_build_dir options dependency_solution lib_src_lock_file lib_build_lock_file )

    message(STATUS "BPM [${PROJECT_NAME}:${lib_name}]: Configuring")

    set(cmake_build_args "")
    if(options)
        foreach(opt IN LISTS options)
            list(APPEND cmake_build_args "-D${opt}")
        endforeach()
    endif()

    set(test_example_options)
    # parse the libraries cmake lists for flags that enable tests and disable them
    bpm_find_test_example_options_r("${lib_src_dir}" test_example_options)
    list(LENGTH test_example_options test_example_options_length)

    set(cmake_disable_test_example_flags)
    foreach(flag ${test_example_options})
        # check if option is already explicitly enabled
        set(found FALSE)
        foreach(opt IN LISTS options)
            if("${opt}" MATCHES "^[ \t]*${flag}[ \t]*=[ \t]*(ON|TRUE|1)$[ \t]*")
                set(found TRUE)
                break()        
            endif()
        endforeach()

        if(NOT found)
            list(APPEND cmake_disable_test_example_flags "-D${flag}=OFF")
        endif()
    endforeach()
    message(STATUS "BPM [${PROJECT_NAME}:${lib_name}]: Found test/example option. Disabiling: '${cmake_disable_test_example_flags}'")

    set(toolchain_args)
    if(CMAKE_TOOLCHAIN_FILE)
        # turn toolchain file path into an absolute path
        set(abs_toolchain_file "${CMAKE_TOOLCHAIN_FILE}")
        cmake_path(ABSOLUTE_PATH abs_toolchain_file BASE_DIRECTORY "${CMAKE_SOURCE_DIR}" NORMALIZE)
        set(toolchain_args "-DCMAKE_TOOLCHAIN_FILE=${abs_toolchain_file}")
    else()
        set(toolchain_args
            "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
            "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
        )
    endif()

    set(config_arg)
    if((NOT CMAKE_CONFIGURATION_TYPES) AND (NOT CMAKE_TOOLCHAIN_FILE))
        set(config_arg "-DCMAKE_BUILD_TYPE=Release")    
    endif()

    set(verbose_arg)
    if(BPM_VERBOSE)
        set(verbose_arg "-DBPM_VERBOSE=ON")
    endif()

    # lockorder: install -> build -> source -> mirror
    file(LOCK "${lib_build_lock_file}")
        file(LOCK "${lib_src_lock_file}")  

            # Check if this library uses BPM (has a .bpm-registry). If yes: pass the version solutions as a variable
            set(dependencies_arg)
            set(bpm_cache_arg)
            if(EXISTS "${lib_src_dir}/.bpm-registry")
                if(BPM_DEPENDENCY_SOLUTION)
                    set(dependencies_arg "-DBPM_DEPENDENCY_SOLUTION=${BPM_DEPENDENCY_SOLUTION}")
                else()
                    set(dependencies_arg "-DBPM_DEPENDENCY_SOLUTION=${CMAKE_BINARY_DIR}/bpm-dependency-solution.txt")
                endif()
                set(bpm_cache_arg "-DBPM_CACHE=${BPM_CACHE_DIR}")
            endif()
            
            set(arg_position_independent_code)
            if(CMAKE_POSITION_INDEPENDENT_CODE)
                set(arg_position_independent_code "-DCMAKE_POSITION_INDEPENDENT_CODE=${CMAKE_POSITION_INDEPENDENT_CODE}")
            endif()

            set(arg_build_shared_libs)
            if(BUILD_SHARED_LIBS)
                set(arg_build_shared_libs "-DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}")
            endif()

            set(quiet)
            if(NOT BPM_VERBOSE)
                set(quiet "OUTPUT_QUIET")
            endif()

            # TODO: Optimisation: skip instead of re-configuring
            if(BPM_VERBOSE)
                message(STATUS "BPM [${PROJECT_NAME}:${lib_name}]: execute command: ${CMAKE_COMMAND} -S \"${lib_src_dir}\" -B \"${lib_build_dir}\" -G \"${CMAKE_GENERATOR}\" ${config_arg} ${bpm_cache_arg} -DCMAKE_INSTALL_PREFIX=\"${lib_install_dir}\" ${arg_position_independent_code} ${arg_build_shared_libs} ${cmake_build_args} ${toolchain_args} ${cmake_disable_test_example_flags} ${dependencies_arg} ${verbose_arg}")
            endif()

            execute_process(
                COMMAND ${CMAKE_COMMAND}
                -S "${lib_src_dir}"
                -B "${lib_build_dir}"
                -G "${CMAKE_GENERATOR}"
                
                "-DCMAKE_MAKE_PROGRAM=${CMAKE_MAKE_PROGRAM}"
                "-DCMAKE_GENERATOR=${CMAKE_GENERATOR}"

                ${bpm_cache_arg}
                "-DCMAKE_INSTALL_PREFIX=${lib_install_dir}"
                ${arg_position_independent_code}
                ${arg_build_shared_libs}
                
                ${cmake_build_args}
                ${toolchain_args}
                ${cmake_disable_test_example_flags}
                ${dependencies_arg}
                
                ${verbose_arg}

                RESULT_VARIABLE res
                ${quiet}
            )

        file(LOCK "${lib_build_lock_file}" RELEASE)
    file(LOCK "${lib_src_lock_file}" RELEASE)

    if(res EQUAL 0)
        if(BPM_VERBOSE)
            message(STATUS "BPM [${PROJECT_NAME}:${lib_name}]: Configuring - done")
        endif()
    else() 
        message(FATAL_ERROR "BPM [${PROJECT_NAME}:${lib_name}]: Configuring - failed")
    endif()

endfunction()

function(bpm_build_library lib_name library_build_dir lib_build_lock_file)

    message(STATUS "BPM [${PROJECT_NAME}:${lib_name}]: Building")

    set(parallel 8)
    if(CMAKE_BUILD_PARALLEL_LEVEL)
        set(parallel ${CMAKE_BUILD_PARALLEL_LEVEL})
    endif()

    set(config_arg)
    if(CMAKE_CONFIGURATION_TYPES)
        set(config_arg --config Release)
    endif()

    file(LOCK "${lib_build_lock_file}")
        message(STATUS "================================ BUILDING ${lib_name} ================================")
        if(BPM_VERBOSE)
            message(STATUS "BPM [${PROJECT_NAME}:${lib_name}]: execute command: ${CMAKE_COMMAND} --build \"${library_build_dir}\" ${config_arg} --parallel ${parallel}")
        endif()
        execute_process(COMMAND ${CMAKE_COMMAND} --build "${library_build_dir}" ${config_arg} --parallel ${parallel} RESULT_VARIABLE res)
        message(STATUS "================================ BUILDING ${lib_name} - FINISHED ================================")
    file(LOCK "${lib_build_lock_file}" RELEASE)

    if(res EQUAL 0)
        if(BPM_VERBOSE)
            message(STATUS "BPM [${PROJECT_NAME}:${lib_name}]: Building - done")
        endif()
    else() 
        message(FATAL_ERROR "BPM [${PROJECT_NAME}:${lib_name}]: Building - failed")
    endif()

endfunction()

function(bpm_install_library lib_name lib_build_dir lib_install_dir lib_build_lock_file lib_install_lock_file)
    if(BPM_VERBOSE)
        message(STATUS "BPM [${PROJECT_NAME}:${lib_name}]: Installing from ${lib_build_dir} into ${lib_install_dir}")
    else()
        message(STATUS "BPM [${PROJECT_NAME}:${lib_name}]: Installing")
    endif()

    set(config_arg)
    if(CMAKE_CONFIGURATION_TYPES)
        set(config_arg --config Release)
    endif()

    # lockorder: install -> build -> source -> mirror
    file(LOCK "${lib_install_lock_file}")
        file(LOCK "${lib_build_lock_file}")
            if(BPM_VERBOSE)
                message(STATUS "BPM [${PROJECT_NAME}:${lib_name}]: execute command: ${CMAKE_COMMAND} --install ${lib_build_dir} --prefix ${lib_install_dir} ${config_arg}")
                execute_process(COMMAND ${CMAKE_COMMAND} --install ${lib_build_dir} --prefix ${lib_install_dir} ${config_arg} RESULT_VARIABLE res)
            else()
                execute_process(COMMAND ${CMAKE_COMMAND} --install ${lib_build_dir} --prefix ${lib_install_dir} ${config_arg} RESULT_VARIABLE res OUTPUT_QUIET ERROR_QUIET)
            endif()
            
        file(LOCK "${lib_build_lock_file}" RELEASE)
    file(LOCK "${lib_install_lock_file}" RELEASE)

    if(res EQUAL 0)
        if(BPM_VERBOSE)
            message(STATUS "BPM [${PROJECT_NAME}:${lib_name}]: Installing - done")
        endif()
    else() 
        # clean install on error
        file(REMOVE_RECURSE "${lib_install_dir}")
        message(FATAL_ERROR "BPM [${PROJECT_NAME}:${lib_name}]: Installing - failed")
    endif()

endfunction()

function(bpm_find_installed_packages PGK_NAME lib_install_dir OUT_PACKAGES)
    if(NOT EXISTS ${lib_install_dir})
        set(${OUT_PACKAGES} "" PARENT_SCOPE)
        return()
    endif()

    file(GLOB_RECURSE config_files "${lib_install_dir}/*Config.cmake" "${lib_install_dir}/*config.cmake")

    set(installed_packages)
    if(config_files)
        foreach(config_file IN LISTS config_files)
            get_filename_component(config_filename "${config_file}" NAME)
            set(package_name "${config_filename}")
            string(REGEX REPLACE "Config\\.cmake$" "" package_name "${package_name}")
            string(REGEX REPLACE "-config\\.cmake$" "" package_name "${package_name}")
            list(APPEND installed_packages "${package_name}")
        endforeach()
    endif()
    set(${OUT_PACKAGES} "${installed_packages}" PARENT_SCOPE)
endfunction()

function(bpm_show_installed_packages PKG_NAME lib_install_dir)
    if(NOT EXISTS ${lib_install_dir})
        message(FATAL_ERROR "BPM [${PROJECT_NAME}:${PKG_NAME}]: Install dir '${lib_install_dir}' does not exist.")
        return()
    endif()
    file(GLOB_RECURSE config_files "${lib_install_dir}/*Config.cmake" "${lib_install_dir}/*config.cmake")
    if(config_files)
        if(BPM_VERBOSE)
            message(STATUS "BPM [${PROJECT_NAME}:${PKG_NAME}]: Installed packages:")
            foreach(config_file IN LISTS config_files)
                get_filename_component(config_dir "${config_file}" DIRECTORY)
                get_filename_component(package_name "${config_dir}" NAME)
                message(STATUS "  - ${package_name}:\t from ${config_file}")
            endforeach()
        else()
            set(installed_package_names)
            foreach(config_file IN LISTS config_files)
                get_filename_component(config_dir "${config_file}" DIRECTORY)
                get_filename_component(package_name "${config_dir}" NAME)
                if(NOT installed_package_names)
                    string(APPEND installed_package_names "${package_name}")
                else()
                    string(APPEND installed_package_names ", ${package_name}")
                endif()
            endforeach()
            message(STATUS "BPM [${PROJECT_NAME}:${PKG_NAME}]: Installed packages: ${installed_package_names}")
        endif()
    else()
        if(BPM_VERBOSE)
            message(STATUS "BPM [${PROJECT_NAME}:${PKG_NAME}]: No packages were found in: ${lib_install_dir}")
        else()
            message(STATUS "BPM [${PROJECT_NAME}:${PKG_NAME}]: No packages were found")
        endif()
    endif()
endfunction()

function(bpm_load_dependencies BPM_CACHE_DIR registry_content master_solution out_solution)

    set(todo_list "${registry_content}")
    set(solution "")

    foreach(todo IN LISTS todo_list)
        # clear name to prevent accidental reuse if parsing fails
        set(TODO_NAME)

        # parse the todo entry
        set(options "")
        set(oneValueArgs NAME)
        set(multiValueArgs)
        separate_arguments(todo_tokens UNIX_COMMAND "${todo}")
        cmake_parse_arguments(TODO "${options}" "${oneValueArgs}" "${multiValueArgs}" ${todo_tokens})

        set(in_todo_${TODO_NAME} TRUE)
    endforeach()

    while(todo_list)
        list(POP_FRONT todo_list todo)

        # clear name to prevent accidental reuse if parsing fails
        set(PKG_NAME)

        set(options "")
        set(oneValueArgs NAME)
        set(multiValueArgs)
        separate_arguments(todo_tokens UNIX_COMMAND "${todo}")
        cmake_parse_arguments(PKG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${todo_tokens})

        set(mirror_lock_file "${BPM_CACHE_DIR}/${PKG_NAME}/mirror.lock")
        set(mirror_dir "${BPM_CACHE_DIR}/${PKG_NAME}/mirror")

        # load solution data from master solution
        foreach(line IN LISTS master_solution)

            # clear name to prevent accidental reuse if parsing fails
            set(SOL_NAME)
            set(SOL_VERSION)

            set(options "")
            set(oneValueArgs NAME VERSION)
            set(multiValueArgs "")
            separate_arguments(line_tokens UNIX_COMMAND "${line}")
            cmake_parse_arguments(SOL "${options}" "${oneValueArgs}" "${multiValueArgs}" ${line_tokens})
            if("${SOL_NAME}" STREQUAL "${PKG_NAME}")
                if(BPM_VERBOSE)
                    message(STATUS "BPM [${PROJECT_NAME}:${PKG_NAME}]: Found solution: ${line}")
                else()
                    message(STATUS "BPM [${PROJECT_NAME}:${PKG_NAME}]: Found solution: NAME ${SOL_NAME} VERSION ${SOL_VERSION} ... ")
                endif()
                # build solution list 
                if(NOT in_solution_${PKG_NAME})
                    set(in_solution_${PKG_NAME} TRUE)
                    list(APPEND solution "${line}")
                endif()

                break()
            endif()
        endforeach()
        
        # load registry and add to the todo list, but only if not already part of the solution or the todo
        if(NOT DEFINED metadata_${SOL_NAME}_${SOL_VERSION})
            file(LOCK "${mirror_lock_file}")
                execute_process(COMMAND git --git-dir "${mirror_dir}" cat-file blob "${SOL_VERSION}:.bpm-registry" RESULT_VARIABLE res OUTPUT_VARIABLE metadata_tmp ERROR_QUIET)
            file(LOCK "${mirror_lock_file}" RELEASE)

            string(REPLACE "\r\n" "\n" metadata_tmp "${metadata_tmp}") # replace new lines windows to unix style
            string(REPLACE "\n" ";" metadata_${PKG_NAME}_${SOL_VERSION} "${metadata_tmp}") # replace new lines with ; for list seperators
        endif()

        # load the metadata for the package and add entries to the todo list
        foreach(line IN LISTS metadata_${PKG_NAME}_${SOL_VERSION})
            if(line)
                string(STRIP "${line}" line)

                # clear name to prevent accidental reuse if parsing fails
                set(LINE_NAME)

                # parse the line to check if it is not already in the solution or the todo list before adding it to the todo list
                set(options "")
                set(oneValueArgs NAME)
                set(multiValueArgs "")
                separate_arguments(line_tokens UNIX_COMMAND "${line}")
                cmake_parse_arguments(LINE "${options}" "${oneValueArgs}" "${multiValueArgs}" ${line_tokens})

                if((NOT line STREQUAL "") AND (NOT in_todo_${LINE_NAME}))
                    set(in_todo_${LINE_NAME} TRUE)
                    list(APPEND todo_list "${line}")
                endif()
            endif()
        endforeach()

    endwhile()

    set(${out_solution} "${solution}" PARENT_SCOPE)
endfunction()

function(bpm_check_for_bpm_updates BPM_VERSION BPM_REPO)
    message(STATUS "BPM [${PROJECT_NAME}]: Checking for BPM updates...")

    execute_process(
        COMMAND git ls-remote --tags --sort=-version:refname "${BPM_REPO}" "refs/tags/*"
        RESULT_VARIABLE res
        OUTPUT_VARIABLE bpm_version_tags
        ERROR_VARIABLE err
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    if(bpm_version_tags MATCHES "v([0-9]+\\.[0-9]+\\.[0-9]+)")
    set(newest_bpm_version "${CMAKE_MATCH_1}")
        if(newest_bpm_version VERSION_GREATER BPM_VERSION)
            message(STATUS "BPM [${PROJECT_NAME}]:   A newer version of BPM is available: ${newest_bpm_version} (current: ${BPM_VERSION})")
            if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
                message(STATUS "BPM [${PROJECT_NAME}]:     Update with: curl -o cmake/BPM.cmake \"https://github.com/TobiasWallner/BPM.cmake/releases/download/${newest_bpm_version}/BPM.cmake\" -L")
            elseif(CMAKE_SYSTEM_NAME STREQUAL "Windows")
                message(STATUS "BPM [${PROJECT_NAME}]:     Update with: Invoke-WebRequest -Uri \"https://github.com/TobiasWallner/BPM.cmake/releases/download/${newest_bpm_version}/BPM.cmake\" -OutFile \"cmake/BPM.cmake\"")
            endif()
        endif()
    endif()        
    
endfunction()

#
function(BPMMakeAvailable)

    set(BPM_VERSION "v0.5.2")
    set(BPM_REPO "https://github.com/TobiasWallner/BPM.cmake")

    message(STATUS "BPM [${PROJECT_NAME}]: BPM version: ${BPM_VERSION}")

    if(NOT PROJECT_IS_TOP_LEVEL)
        return()
    endif()

    message("")

    set(options VERBOSE NO_DOWNLOAD NO_DOWNLOAD_UPDATES)
    set(oneValueArgs)
    set(multiValueArgs)
    cmake_parse_arguments(_BPM "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    if(_BPM_VERBOSE)
        set(BPM_VERBOSE TRUE)
    endif()
    if(_BPM_NO_DOWNLOAD)
        set(BPM_NO_DOWNLOAD TRUE)
    endif()
    if(_BPM_NO_DOWNLOAD_UPDATES)
        set(BPM_NO_DOWNLOAD_UPDATES TRUE)
    endif()

    if((NOT BPM_NO_DOWNLOAD) AND (NOT BPM_NO_DOWNLOAD_UPDATES))
        bpm_check_for_bpm_updates("${BPM_VERSION}" "${BPM_REPO}")
    endif()

    bpm_get_cache_dir(BPM_CACHE_DIR)

    # -------------------------------
    # Write local repository to file
    # -------------------------------

    get_property(BPM_REGISTRY_ GLOBAL PROPERTY BPM_REGISTRY)
    foreach(PKG_NAME IN LISTS BPM_REGISTRY_)
        # write local registry
        get_property(VERSION_RANGE GLOBAL PROPERTY "BPM_REGISTRY_${PKG_NAME}_VERSION_RANGE")
        get_property(GIT_TAG GLOBAL PROPERTY "BPM_REGISTRY_${PKG_NAME}_GIT_TAG")
        get_property(GIT_REPOSITORY GLOBAL PROPERTY "BPM_REGISTRY_${PKG_NAME}_GIT_REPOSITORY")
        get_property(TYPE GLOBAL PROPERTY "BPM_REGISTRY_${PKG_NAME}_TYPE")
        get_property(OPTIONS GLOBAL PROPERTY "BPM_REGISTRY_${PKG_NAME}_OPTIONS")
        get_property(PKG_PACKAGES GLOBAL PROPERTY "BPM_REGISTRY_${PKG_NAME}_PACKAGES")
        get_property(PKG_PRIVATE GLOBAL PROPERTY "BPM_REGISTRY_${PKG_NAME}_PRIVATE")

        string(REPLACE ";" "-" SAFE_VERSION_RANGE "${VERSION_RANGE}")

        set(options_entry)
        if(OPTIONS)
            string(REPLACE ";" " " options_entry "${OPTIONS}")
            set(options_entry "OPTIONS ${options_entry}")
        endif()

        set(packages_entry)
        if(PKG_PACKAGES)
            string(REPLACE ";" " " packages_entry "${PKG_PACKAGES}")
            set(packages_entry "PACKAGES ${packages_entry}")
        endif()

        set(entry "NAME ${PKG_NAME} VERSION_RANGE ${SAFE_VERSION_RANGE} GIT_TAG ${GIT_TAG} GIT_REPOSITORY ${GIT_REPOSITORY} TYPE ${TYPE} ${options_entry} ${packages_entry}")
        string(APPEND registry_content "${entry};")
        
        # Do not add the entry to the external registry file if the package is marked as private, to prevent leaking private dependencies to other projects that use this project as a dependency
        if(NOT PKG_PRIVATE)
            string(APPEND registry_file_content "${entry}\n")
        endif()

    endforeach()
    
    # TODO: sort fily by package names before writing for improved robustness
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.bpm-registry")
        file(READ "${CMAKE_CURRENT_SOURCE_DIR}/.bpm-registry" old_registry_file_content)
        if(NOT "${registry_file_content}\n" STREQUAL "${old_registry_file_content}")
            file(WRITE "${CMAKE_CURRENT_SOURCE_DIR}/.bpm-registry" "${registry_file_content}\n")
        endif()
    else()
        file(WRITE "${CMAKE_CURRENT_SOURCE_DIR}/.bpm-registry" "${registry_file_content}\n")
    endif()

    # -------------------------------
    # Solve Dependency Graph
    # -------------------------------

    if(NOT DEFINED BPM_DEPENDENCY_SOLUTION)
        message("")
        message(STATUS "BPM [${PROJECT_NAME}]: Solving dependency graph")
        bpm_solve_dependencies("${BPM_CACHE_DIR}" "${registry_content}" solution)

        # write/update solution on change
        string(REPLACE ";" "\n" solution_file_content "${solution}")
        file(WRITE "${CMAKE_BINARY_DIR}/bpm-dependency-solution.txt" "${solution_file_content}")

        message("")
        message(STATUS "BPM [${PROJECT_NAME}]: Dependency Graph Solution ")

        foreach(pkg IN LISTS solution)
            # clear variables to prevent accidental reuse in the loop
            set(PKG_NAME)
            set(PKG_VERSION)
            set(PKG_GIT_REPO)
            set(PKG_OPTIONS)
            set(PKG_PACKAGES)

            # parse solution entry
            set(options)
            set(oneValueArgs NAME VERSION GIT_REPO)
            set(multiValueArgs OPTIONS PACKAGES DEPENDENCIES)
            separate_arguments(pkg_tokens UNIX_COMMAND "${pkg}")
            cmake_parse_arguments(PKG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${pkg_tokens})

            message(STATUS "  + Resolved: ${PKG_NAME}#${PKG_VERSION}")
            if(BPM_VERBOSE)
                message(STATUS "    > GIT_REPO: ${PKG_GIT_REPO}")
                if(PKG_OPTIONS)
                    message(STATUS "    > OPTIONS:")
                    foreach(opt IN LISTS PKG_OPTIONS)
                        message(STATUS "      - ${opt}")
                    endforeach()
                endif()

                if(PKG_PACKAGES)
                    message(STATUS "    > PACKAGES:")
                    foreach(p IN LISTS PKG_PACKAGES)
                        message(STATUS "      - ${p}")
                    endforeach()
                endif()
            endif()
        endforeach()

    else()
        file(READ "${BPM_DEPENDENCY_SOLUTION}" master_solution)
        string(REPLACE "\r\n" "\n" master_solution "${master_solution}")
        string(REPLACE "\n" ";" master_solution "${master_solution}")
        bpm_load_dependencies("${BPM_CACHE_DIR}" "${registry_content}" "${master_solution}" solution)

    endif()

    # ------------------------------------------------------------------------------
    # Add package with `add_subdirectory` or install and add with `find_package`
    # ------------------------------------------------------------------------------
    
    # Provide key value pairs for solution entries
    foreach(pkg IN LISTS solution) 
        # clear variables to prevent accidental reuse in the loop
        set(PKG_NAME)
        set(PKG_VERSION)
        set(PKG_GIT_REPO)
        set(PKG_TYPE)
        set(PKG_OPTIONS)
        set(PKG_PACKAGES)

        # parse solution entry
        set(options)
        set(oneValueArgs NAME VERSION GIT_REPO TYPE)
        set(multiValueArgs OPTIONS PACKAGES DEPENDENCIES)
        separate_arguments(pkg_tokens UNIX_COMMAND "${pkg}")
        cmake_parse_arguments(PKG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${pkg_tokens})

        set(BPM_PKG_${PKG_NAME}_VERSION "${PKG_VERSION}")
        set(BPM_PKG_${PKG_NAME}_GIT_REPO "${PKG_GIT_REPO}")
        set(BPM_PKG_${PKG_NAME}_OPTIONS "${PKG_OPTIONS}")
        set(BPM_PKG_${PKG_NAME}_PACKAGES "${PKG_PACKAGES}")

    endforeach()

    message("")
    message(STATUS "BPM [${PROJECT_NAME}]: Making packages available")
    list(REVERSE solution) # reverse for correct order of installation (dependencies first) - should have leaves first
    foreach(pkg IN LISTS solution)
        # clear variables to prevent accidental reuse in the loop
        set(PKG_NAME)
        set(PKG_VERSION)
        set(PKG_GIT_REPO)
        set(PKG_TYPE)
        set(PKG_OPTIONS)
        set(PKG_PACKAGES)

        # parse solution entry
        set(options)
        set(oneValueArgs NAME VERSION GIT_REPO TYPE)
        set(multiValueArgs OPTIONS PACKAGES DEPENDENCIES)
        separate_arguments(pkg_tokens UNIX_COMMAND "${pkg}")
        cmake_parse_arguments(PKG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${pkg_tokens})

        # Prevent double adding through subdirectory packages
        get_property(PKG_MADE_AVAILABLE GLOBAL PROPERTY "BPM_REGISTRY_${PKG_NAME}_MADE_AVAILABLE")
        if(PKG_MADE_AVAILABLE)
            continue()
        endif()
        set_property(GLOBAL PROPERTY "BPM_REGISTRY_${PKG_NAME}_MADE_AVAILABLE" TRUE)

        set(lib_mirror_dir "${BPM_CACHE_DIR}/${PKG_NAME}/mirror")
        set(lib_mirror_lock_file "${BPM_CACHE_DIR}/${PKG_NAME}/mirror.lock")

        # -------------------------------
        # create manifest
        # -------------------------------

        set(C_COMPILER_HASH)
        set(C_COMPILER_HASH)
        if(CMAKE_C_COMPILER)
            if(IS_ABSOLUTE "${CMAKE_C_COMPILER}")
                # Already an absolute path
                if(EXISTS "${CMAKE_C_COMPILER}")
                    file(SHA256 "${CMAKE_C_COMPILER}" C_COMPILER_HASH)
                endif()
            else()
                # Command name only - resolve it via PATH
                find_program(C_COMPILER_FULL_PATH "${CMAKE_C_COMPILER}")
                if(C_COMPILER_FULL_PATH AND EXISTS "${C_COMPILER_FULL_PATH}")
                    file(SHA256 "${C_COMPILER_FULL_PATH}" C_COMPILER_HASH)
                endif()
            endif()
        endif()

        set(CXX_COMPILER_HASH)
        if(CMAKE_CXX_COMPILER)
            if(IS_ABSOLUTE "${CMAKE_CXX_COMPILER}")
                # Already an absolute path
                if(EXISTS "${CMAKE_CXX_COMPILER}")
                    file(SHA256 "${CMAKE_CXX_COMPILER}" CXX_COMPILER_HASH)
                endif()
            else()
                # Command name only - resolve it via PATH
                find_program(CXX_COMPILER_FULL_PATH "${CMAKE_CXX_COMPILER}")
                if(CXX_COMPILER_FULL_PATH AND EXISTS "${CXX_COMPILER_FULL_PATH}")
                    file(SHA256 "${CXX_COMPILER_FULL_PATH}" CXX_COMPILER_HASH)
                endif()
            endif()
        endif()

        file(LOCK "${lib_mirror_lock_file}")
            execute_process(COMMAND git --git-dir "${lib_mirror_dir}" rev-parse "${PKG_VERSION}^{commit}" RESULT_VARIABLE res OUTPUT_VARIABLE PKG_GIT_COMMIT OUTPUT_STRIP_TRAILING_WHITESPACE)
        file(LOCK "${lib_mirror_lock_file}" RELEASE)

        if(NOT res EQUAL 0)
            message(FATAL_ERROR "BPM [${PROJECT_NAME}:${PKG_NAME}]: Could not convert '${PKG_VERSION}' to commit-hash in mirror '${lib_mirror_dir}'")
        endif()

        # get a list of all dependency (deep information) solutions for the manifest
        set(DEPENDENCIES_MANIFESTS)
        list(SORT PKG_DEPENDENCIES)
        foreach(dep IN LISTS PKG_DEPENDENCIES)
            if(BPM_${dep}_FOUND)
                string(APPEND DEPENDENCIES_MANIFESTS "\n  DEPENDENCY ${dep} VERSION ${BPM_PKG_${PKG_NAME}_VERSION} REPO ${BPM_PKG_${PKG_NAME}_GIT_REPO} OPTIONS ${BPM_PKG_${PKG_NAME}_OPTIONS} PACKAGES ${BPM_PKG_${PKG_NAME}_PACKAGES} MANIFEST ${BPM_${dep}_MANIFEST_HASH}")
            else()
                message(FATAL_ERROR "BPM [${PROJECT_NAME}:${PKG_NAME}]: Dependency '${dep}' referenced but no manifest has been created yet. This should not happen. Please report this to the developers.")
            endif()
        endforeach()

        # sort options list before creating the manifest
        list(SORT PKG_OPTIONS)
        list(SORT PKG_PACKAGES)
        
        
        # turn tag into commit hash
        bpm_create_manifest(manifest
            CMAKE_C_COMPILER_ID
            C_COMPILER_HASH
            CMAKE_C_COMPILER_VERSION
            CMAKE_CXX_COMPILER_ID
            CXX_COMPILER_HASH
            CMAKE_CXX_COMPILER_VERSION
            CMAKE_SYSTEM_NAME
            CMAKE_SYSTEM_PROCESSOR
            CMAKE_VERSION
            CMAKE_GENERATOR
            CMAKE_GENERATOR_PLATFORM
            CMAKE_GENERATOR_TOOLSET
            BUILD_SHARED_LIBS
            CMAKE_POSITION_INDEPENDENT_CODE
            CMAKE_INTERPROCEDURAL_OPTIMIZATION
            CMAKE_C_FLAGS
            CMAKE_CXX_FLAGS
            CMAKE_EXE_LINKER_FLAGS
            CMAKE_SHARED_LINKER_FLAGS
            TOOLCHAIN_HASH
            PKG_NAME
            PKG_VERSION
            PKG_GIT_COMMIT
            PKG_GIT_REPO
            PKG_OPTIONS
            PKG_TYPE
            DEPENDENCIES_MANIFESTS
        ) 

        string(SHA256 manifest_hash "${manifest}")
        string(SUBSTRING "${manifest_hash}" 0 16 SHORT_MANIFEST_HASH)
        string(SUBSTRING "${PKG_GIT_COMMIT}" 0 16 PKG_GIT_COMMIT_HASH)

        # -------------------------------
        # Define hashed directories
        # -------------------------------

        set(lib_src_dir "${BPM_CACHE_DIR}/${PKG_NAME}/src/${PKG_GIT_COMMIT_HASH}")
        set(lib_src_lock_file "${BPM_CACHE_DIR}/${PKG_NAME}/src/${PKG_GIT_COMMIT_HASH}.lock")

        set(lib_build_dir "${BPM_CACHE_DIR}/${PKG_NAME}/build/${SHORT_MANIFEST_HASH}")
        set(lib_build_lock_file "${BPM_CACHE_DIR}/${PKG_NAME}/build/${SHORT_MANIFEST_HASH}.lock")

        set(lib_install_dir "${BPM_CACHE_DIR}/${PKG_NAME}/install/${SHORT_MANIFEST_HASH}")
        set(lib_install_lock_file "${BPM_CACHE_DIR}/${PKG_NAME}/install/${SHORT_MANIFEST_HASH}.lock")

        set(manifest_dir "${BPM_CACHE_DIR}/${PKG_NAME}/manifest")
        set(manifest_file_path "${BPM_CACHE_DIR}/${PKG_NAME}/manifest/${SHORT_MANIFEST_HASH}.manifest")   

        set(lib_src_dir "${BPM_CACHE_DIR}/${PKG_NAME}/src/${PKG_GIT_COMMIT_HASH}")
        set(lib_src_lock_file "${BPM_CACHE_DIR}/${PKG_NAME}/src/${PKG_GIT_COMMIT_HASH}.lock")

        set(lib_build_dir "${BPM_CACHE_DIR}/${PKG_NAME}/build/${SHORT_MANIFEST_HASH}")
        set(lib_build_lock_file "${BPM_CACHE_DIR}/${PKG_NAME}/build/${SHORT_MANIFEST_HASH}.lock")

        set(lib_install_dir "${BPM_CACHE_DIR}/${PKG_NAME}/install/${SHORT_MANIFEST_HASH}")
        set(lib_install_lock_file "${BPM_CACHE_DIR}/${PKG_NAME}/install/${SHORT_MANIFEST_HASH}.lock")

        set(manifest_dir "${BPM_CACHE_DIR}/${PKG_NAME}/manifest")
        set(manifest_file_name "${SHORT_MANIFEST_HASH}.manifest")
        set(manifest_file_path "${BPM_CACHE_DIR}/${PKG_NAME}/manifest/${SHORT_MANIFEST_HASH}.manifest")

        set(BPM_${PKG_NAME}_FOUND TRUE PARENT_SCOPE)
        set(BPM_${PKG_NAME}_FOUND TRUE) # set in current scope too
        set(BPM_${PKG_NAME}_SRC_DIR "${lib_src_dir}" PARENT_SCOPE)
        set(BPM_${PKG_NAME}_SRC_DIR "${lib_src_dir}")
        set(BPM_${PKG_NAME}_BUILD_DIR "${lib_build_dir}" PARENT_SCOPE)
        set(BPM_${PKG_NAME}_BUILD_DIR "${lib_build_dir}")
        set(BPM_${PKG_NAME}_INSTALL_DIR "${lib_install_dir}" PARENT_SCOPE)
        set(BPM_${PKG_NAME}_INSTALL_DIR "${lib_install_dir}")
        set(BPM_${PKG_NAME}_MANIFEST_HASH "${manifest_hash}" PARENT_SCOPE)
        set(BPM_${PKG_NAME}_MANIFEST_HASH "${manifest_hash}")
        set(BPM_${PKG_NAME}_MANIFEST_SHORT_HASH "${SHORT_MANIFEST_HASH}" PARENT_SCOPE)
        set(BPM_${PKG_NAME}_MANIFEST_SHORT_HASH "${SHORT_MANIFEST_HASH}")
        set(BPM_${PKG_NAME}_MANIFEST_FILE_NAME "${manifest_file_name}" PARENT_SCOPE)
        set(BPM_${PKG_NAME}_MANIFEST_FILE_NAME "${manifest_file_name}")
        set(BPM_${PKG_NAME}_MANIFEST_FILE_PATH "${manifest_file_path}" PARENT_SCOPE)
        set(BPM_${PKG_NAME}_MANIFEST_FILE_PATH "${manifest_file_path}")

        list(PREPEND CMAKE_PREFIX_PATH "${lib_install_dir}")
        list(REMOVE_DUPLICATES CMAKE_PREFIX_PATH)
        set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH}" PARENT_SCOPE)

        # write manifest
        if(NOT EXISTS ${manifest_dir})
            file(MAKE_DIRECTORY "${manifest_dir}")
        endif()

        if(NOT EXISTS ${manifest_file_path})
            file(WRITE "${manifest_file_path}" "${manifest}")
        endif()

        if("${PKG_TYPE}" STREQUAL "INSTALL")
            if(BPM_VERBOSE)
                message(STATUS "BPM [${PROJECT_NAME}:${PKG_NAME}]: Make Available (INSTALL): ${PKG_NAME}#${PKG_VERSION} : ${PKG_GIT_REPO}")
            endif()

            
            if(NOT PKG_PACKAGES)
                if(BPM_VERBOSE)
                    message(STATUS "BPM [${PROJECT_NAME}:${PKG_NAME}]: No PACKAGES provided for INSTALL type - assuming package name from library name: PACKAGES=${PKG_NAME}")
                endif()
                set(PKG_PACKAGES "${PKG_NAME}")
            endif()

            bpm_try_find_packages("${PKG_NAME}" "${PKG_PACKAGES}" "${lib_install_dir}" all_packages_found)
            if(all_packages_found)
                foreach(package IN LISTS PKG_PACKAGES)
                    # export variables to parent scope for later use (e.g. for generating manifest)
                    foreach(var IN ITEMS FOUND DIR VERSION CONFIG)
                        if(DEFINED ${package}_${var})
                            set(${package}_${var} "${${package}_${var}}" PARENT_SCOPE)
                        endif()
                    endforeach()
                endforeach()
            else()
                
                set(packages_string)
                foreach(p IN LISTS PKG_PACKAGES)
                    if(NOT packages_string)
                        set(packages_string "${p}")
                    else()
                        string(APPEND packages_string ", ${p}")
                    endif()
                endforeach()
                message(STATUS "BPM [${PROJECT_NAME}:${PKG_NAME}]: Attempt installation: ${packages_string}")

                # clone mirror into source dir
                bpm_clone_from_mirror("${PKG_NAME}" "${lib_mirror_dir}" "${lib_src_dir}" "${PKG_VERSION}" "${lib_mirror_lock_file}" "${lib_src_lock_file}")

                # configure the project
                bpm_configure_library("${BPM_CACHE_DIR}" "${PKG_NAME}" "${lib_src_dir}" "${lib_build_dir}" "${PKG_OPTIONS}" "${solution}" "${lib_src_lock_file}" "${lib_build_lock_file}")

                # build the library
                bpm_build_library("${PKG_NAME}" "${lib_build_dir}" "${lib_build_lock_file}")

                # install the library
                bpm_install_library("${PKG_NAME}" "${lib_build_dir}" "${lib_install_dir}" "${lib_build_lock_file}" "${lib_install_lock_file}")

                # print which packages have been installed
                bpm_show_installed_packages("${PKG_NAME}" "${lib_install_dir}")

                # try to make the packagse available
                bpm_try_find_packages("${PKG_NAME}" "${PKG_PACKAGES}" "${lib_install_dir}" all_packages_found)

                # export variables to parent scope for later use (e.g. for generating manifest)
                if(all_packages_found)
                    foreach(package IN LISTS PKG_PACKAGES)
                        foreach(var IN ITEMS FOUND DIR VERSION CONFIG)
                            if(DEFINED ${package}_${var})
                                set(${package}_${var} "${${package}_${var}}" PARENT_SCOPE)
                            endif()
                        endforeach()
                    endforeach()
                endif()

                if(NOT all_packages_found)
                    bpm_find_installed_packages("${PKG_NAME}" "${lib_install_dir}" installed_packages)
                    message(STATUS "BPM [${PROJECT_NAME}:${PKG_NAME}]: Available packages after installation:")
                    foreach(p IN LISTS installed_packages)
                        message(STATUS "  - ${p}")
                    endforeach()
                    message(FATAL_ERROR "BPM [${PROJECT_NAME}:${PKG_NAME}]: Find packages '${PKG_PACKAGES}' - failed after (re-)install.")
                endif()

                # delete source and build directory
                bpm_load_env_var(BPM_CLEAN_SOURCE_AFTER_INSTALL TRUE)
                if(BPM_CLEAN_SOURCE_AFTER_INSTALL)
                    file(LOCK "${lib_src_lock_file}")
                        set(relative_source_dir "\${BPM_CACHE}/${PKG_NAME}/src/${PKG_GIT_COMMIT_HASH}")
                        message(STATUS "BPM [${PROJECT_NAME}:${PKG_NAME}]: Cleaning source dir: '${relative_source_dir}'")
                        file(REMOVE_RECURSE "${lib_src_dir}")
                    file(LOCK "${lib_src_lock_file}" RELEASE)
                    file(REMOVE "${lib_src_dir}")
                    if(BPM_VERBOSE)
                        message(STATUS "BPM [${PROJECT_NAME}:${PKG_NAME}]: Cleaning source dir - done")
                    endif()
                endif()

                bpm_load_env_var(BPM_CLEAN_BUILD_AFTER_INSTALL TRUE)
                if(BPM_CLEAN_BUILD_AFTER_INSTALL)
                    file(LOCK "${lib_build_lock_file}")
                        set(relative_build_dir "\${BPM_CACHE}/${PKG_NAME}/build/${SHORT_MANIFEST_HASH}")
                        message(STATUS "BPM [${PROJECT_NAME}:${PKG_NAME}]: Cleaning build dir: '${relative_build_dir}'")
                        file(REMOVE_RECURSE "${lib_build_dir}")
                    file(LOCK "${lib_build_lock_file}" RELEASE)
                    file(REMOVE "${lib_build_dir}")
                    if(BPM_VERBOSE)
                        message(STATUS "BPM [${PROJECT_NAME}:${PKG_NAME}]: Cleaning build dir - done")
                    endif()
                endif()
                
            endif()

        elseif("${PKG_TYPE}" STREQUAL "ADD_SUBDIR")
            if(BPM_VERBOSE)
                message(STATUS "BPM [${PROJECT_NAME}:${PKG_NAME}]: Make Available (SUBDIRECTORY): ${PKG_NAME}#${PKG_VERSION} : ${PKG_GIT_REPO}")
            endif()

            # clone mirror into source dir
            bpm_clone_from_mirror("${PKG_NAME}" "${lib_mirror_dir}" "${lib_src_dir}" "${PKG_VERSION}" "${lib_mirror_lock_file}" "${lib_src_lock_file}")

            # provide options
            foreach(opt IN LISTS PKG_OPTIONS)
                string(REPLACE "=" ";" opt_list "${opt}")
                list(GET opt_list 0 opt_name)
                list(GET opt_list 1 opt_value)
                set("${opt_name}" "${opt_value}")
            endforeach()

            # provide solution
            set(BPM_DEPENDENCY_SOLUTION "${CMAKE_BINARY_DIR}/bpm-dependency-solution.txt")

            # provide cache dir
            set(BPM_CACHE "${BPM_CACHE_DIR}")

            # disable test and example build targets if possible
            bpm_find_test_example_options_r("${lib_src_dir}" test_example_options)
            foreach(flag IN LISTS test_example_options)
                set(${flag} OFF)
            endforeach()

            add_subdirectory("${lib_src_dir}" "${lib_build_dir}")
        else()
            message(FATAL_ERROR "BPM [${PROJECT_NAME}:${PKG_NAME}]: Internal error. Unknown TYPE '${PKG_TYPE}'. Shoule be 'INSTALL' or 'ADD_SUBDIR'")
        endif()
    endforeach()

endfunction()


# -----------------------------------------------------
#                   Install
# -----------------------------------------------------

function(BPMCreateInstallPackage)
    include(CMakePackageConfigHelpers)

    if(ARGC EQUAL 1)
        # infere everything from the passed library target
        set(ARG_PACKAGE_NAME ${ARGV0})
        set(ARG_EXPORT_NAMESPACE ${ARGV0})
        set(ARG_TARGETS ${ARGV0})
    else()
        # provide specific arguments
        set(options "")
        set(oneValueArgs PACKAGE_NAME NAMESPACE)
        set(multiValueArgs TARGETS PUBLIC_INCLUDE_DIRS HEADER_FILES_MATCHING)
        cmake_parse_arguments(ARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
    endif()

    if(NOT ARG_PACKAGE_NAME)
        message(FATAL_ERROR "BPMInstall: No NAME provided for the package")
    endif()
    
    if(NOT ARG_TARGETS)
        message(FATAL_ERROR "BPMInstall [${ARG_PACKAGE_NAME}]: No TARGETS provided for the package")
    endif()
    
    if(NOT ARG_EXPORT_NAMESPACE)
        set(ARG_EXPORT_NAMESPACE ${PROJECT_NAME})
    endif()

    if(NOT ARG_HEADER_FILES_MATCHING)
        set(ARG_HEADER_FILES_MATCHING "*.h" "*.hh" "*.hpp" "*.hxx")
    endif()

    if(NOT ARG_PUBLIC_INCLUDE_DIRS)
        set(ARG_PUBLIC_INCLUDE_DIRS "include")
    endif()

    install(TARGETS ${ARG_TARGETS}
        EXPORT ${ARG_PACKAGE_NAME}_export_set
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin
        INCLUDES DESTINATION include
    )
    
    set(files_matching "")
    foreach(m IN LISTS ARG_HEADER_FILES_MATCHING)
        LIST(APPEND files_matching "PATTERN" )
        LIST(APPEND files_matching "${m}")
    endforeach()
    
    foreach(dir IN LISTS ARG_PUBLIC_INCLUDE_DIRS)
        # add '/' to the end of the string if it does not have one already
        if(NOT dir MATCHES "/$")
            set(dir "${dir}/")
        endif()

        # add the include directory
        install(
            DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${dir}"
            DESTINATION "${dir}"
            FILES_MATCHING ${files_matching}
        )
    endforeach()

    install(EXPORT ${ARG_PACKAGE_NAME}_export_set
        FILE "${ARG_PACKAGE_NAME}Targets.cmake"
        NAMESPACE "${ARG_EXPORT_NAMESPACE}::"
        DESTINATION "lib/cmake/${ARG_PACKAGE_NAME}"
    )

    set(config_file_in "${CMAKE_CURRENT_BINARY_DIR}/${ARG_PACKAGE_NAME}Config.cmake.in")
    set(config_file "${CMAKE_CURRENT_BINARY_DIR}/${ARG_PACKAGE_NAME}Config.cmake")
    set(version_file "${CMAKE_CURRENT_BINARY_DIR}/${ARG_PACKAGE_NAME}ConfigVersion.cmake")
    file(WRITE "${config_file_in}" "@PACKAGE_INIT@\n\ninclude(\"\${CMAKE_CURRENT_LIST_DIR}/${ARG_PACKAGE_NAME}Targets.cmake\")\n")

    configure_package_config_file(
        "${config_file_in}"
        "${config_file}"
        INSTALL_DESTINATION "lib/cmake/${ARG_PACKAGE_NAME}"
    )

    install(FILES
        "${config_file}"
        DESTINATION "lib/cmake/${ARG_PACKAGE_NAME}"
    )

    if(PROJECT_VERSION)
        write_basic_package_version_file(
            "${version_file}"
            VERSION "${PROJECT_VERSION}"
            COMPATIBILITY SameMajorVersion
        )

        install(FILES
            "${version_file}"
            DESTINATION "lib/cmake/${ARG_PACKAGE_NAME}"
        )
    endif()

    
endfunction()

