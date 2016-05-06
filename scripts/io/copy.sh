#!/usr/bin/env bash

for i in "$@"
do
case ${i} in
    # Files to move
    -s=*)
    SOURCE="${i#*=}"
    shift
    ;;

    -d=*)
    TARGET="${i#*=}"
    shift
    ;;

    # Unknown
    *)
    ;;
esac
done

cp -rf "${SOURCE}" "${TARGET}"