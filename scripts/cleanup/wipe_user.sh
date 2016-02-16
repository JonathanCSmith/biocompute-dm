#!/usr/bin/env bash
for i in "$@"
do
case ${i} in
    # Username
    -u=*|-user=*)
    USERNAME="biocompute-DM_user_${i#*=}"
    shift
    ;;

    # Unknown
    *)
    ;;
esac
done

# Delete the user
deluser "${USERNAME}"
exit
