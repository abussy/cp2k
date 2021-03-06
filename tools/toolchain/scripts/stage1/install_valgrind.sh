#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=SC1003,SC1035,SC1083,SC1090
# shellcheck disable=SC2001,SC2002,SC2005,SC2016,SC2091,SC2034,SC2046,SC2086,SC2089,SC2090
# shellcheck disable=SC2124,SC2129,SC2144,SC2153,SC2154,SC2155,SC2163,SC2164,SC2166
# shellcheck disable=SC2235,SC2237

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

valgrind_ver="3.16.1"
valgrind_sha256="c91f3a2f7b02db0f3bc99479861656154d241d2fdb265614ba918cc6720a33ca"
source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_valgrind" ] && rm "${BUILDDIR}/setup_valgrind"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_valgrind" in
  __INSTALL__)
    echo "==================== Installing Valgrind ===================="
    pkg_install_dir="${INSTALLDIR}/valgrind-${valgrind_ver}"
    install_lock_file="$pkg_install_dir/install_successful"
    if verify_checksums "${install_lock_file}"; then
      echo "valgrind-${valgrind_ver} is already installed, skipping it."
    else
      if [ -f valgrind-${valgrind_ver}.tar.bz2 ]; then
        echo "valgrind-${valgrind_ver}.tar.bz2 is found"
      else
        download_pkg ${DOWNLOADER_FLAGS} ${valgrind_sha256} \
          https://www.cp2k.org/static/downloads/valgrind-${valgrind_ver}.tar.bz2
      fi
      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d valgrind-${valgrind_ver} ] && rm -rf valgrind-${valgrind_ver}
      tar -xjf valgrind-${valgrind_ver}.tar.bz2
      cd valgrind-${valgrind_ver}
      ./configure --prefix="${pkg_install_dir}" --libdir="${pkg_install_dir}/lib" > configure.log 2>&1
      make -j $(get_nprocs) > make.log 2>&1
      make -j $(get_nprocs) install > install.log 2>&1
      cd ..
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage1/$(basename ${SCRIPT_NAME})"
    fi
    ;;
  __SYSTEM__)
    echo "==================== Finding Valgrind from system paths ===================="
    check_command valgrind "valgrind"
    ;;
  __DONTUSE__) ;;

  *)
    echo "==================== Linking Valgrind to user paths ===================="
    pkg_install_dir="$with_valgrind"
    check_dir "${with_valgrind}/bin"
    check_dir "${with_valgrind}/lib"
    check_dir "${with_valgrind}/include"
    ;;
esac
if [ "$with_valgrind" != "__DONTUSE__" ]; then
  if [ "$with_valgrind" != "__SYSTEM__" ]; then
    cat << EOF > "${BUILDDIR}/setup_valgrind"
prepend_path PATH "$pkg_install_dir/bin"
prepend_path PATH "$pkg_install_dir/lib"
prepend_path PATH "$pkg_install_dir/include"
EOF
    cat "${BUILDDIR}/setup_valgrind" >> ${SETUPFILE}
  fi
fi
cd "${ROOTDIR}"

# ----------------------------------------------------------------------
# Suppress reporting of known leaks
# ----------------------------------------------------------------------
cat << EOF >> ${INSTALLDIR}/valgrind.supp
{
   BuggySUPERLU
   Memcheck:Cond
   ...
   fun:SymbolicFactorize
}
{
   BuggyMPICH32
   Memcheck:Cond
   ...
   fun:MPIR_Process_status
}
{
   BuggyLD
   Memcheck:Cond
   ...
   fun:expand_dynamic_string_token
}
EOF
# also need to give links to the .supp file in setup file
cat << EOF >> ${SETUPFILE}
export VALGRIND_OPTIONS="--suppressions=${INSTALLDIR}/valgrind.supp --max-stackframe=2168152 --error-exitcode=42"
EOF

load "${BUILDDIR}/setup_valgrind"
write_toolchain_env "${INSTALLDIR}"

report_timing "valgrind"
