#!/usr/bin/env bash
# Extract common file formats

# Display usage if no parameters given
if [[ -z "$@" ]]; then
  echo " ${0##*/} <archive> - extract common file formats)"
  exit
fi

# Required program(s)
req_progs=(7z unrar unzip)
for p in ${req_progs[@]}; do
  hash "$p" 2>&- || \
  { echo " Required program \"$p\" not installed."  >&2; exit 1; }
done

# Test if file exists
if [ ! -f "$@" ]; then
  echo "File "$@" doesn't exist"
  exit
fi

# Extract file by using extension as reference
case "$@" in
  *.7z      ) 7z  x       "$@"                ;;
  *.tar.bz2 ) tar xvjf    "$@"                ;;
  *.bz2     ) bunzip2     "$@"                ;;
  *.deb     ) ar  vx      "$@"                ;;
  *.tar.gz  ) tar xzvf    "$@"                ;;
  *.gz      ) gunzip      "$@"                ;;
  *.tar     ) tar xvf     "$@"                ;;
  *.tbz2    ) tar xvjf    "$@"                ;;
  *.tar.xz  ) tar xvf     "$@"                ;;
  *.tgz     ) tar xvzf    "$@"                ;;
  *.rar     ) unrar x     "$@"                ;;
  #*.rar     ) 7z  x       "$@"                ;;
  *.zip     ) unzip       "$@"                ;;
  *.Z       ) uncompress  "$@"                ;;
  *         )             "$@"                ;;
esac