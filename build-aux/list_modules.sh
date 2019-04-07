#
# Determine colon-separated list of loaded environment modules

if test "x$MODULESHOME" != "x" ; then

  cs_env_modules=""
  try_modules=""
  try_modules_p=""

  outfile=$1

(
    oldIFS=$IFS; IFS=:
    for m in $LOADEDMODULES; do try_modules="$try_modules $m"; done
    IFS=$oldIFS

    # If using LMOD, re-load associated scripts if possible
    # as top-level (dash) shell may have removed functions...

    if test "x$LMOD_PKG" != "x" ; then
      if test -f "$LMOD_PKG"/init/profile ; then
        source "$LMOD_PKG"/init/profile
      fi
    fi

    module purge

    while test "x$try_modules" != "x$try_modules_p" ;
    do
      try_modules_p=$try_modules
      try_modules=""
      for m in $try_modules_p ; do
        prv_LOADED=$LOADEDMODULES
        module load $m > /dev/null 2>&1
        if test "$prv_LOADED" != "$LOADEDMODULES" ; then
          cs_env_modules="$cs_env_modules $m"
        else
          try_modules="$retry_modules $m"
        fi
      done
    done
    echo "$cs_env_modules" > $outfile
    module list
)

fi
