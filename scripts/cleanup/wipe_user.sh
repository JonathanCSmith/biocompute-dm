#!/usr/bin/env bash
for i in "$@"
do
case ${i} in
    # Username
    -u=*|-user=*)
    NAME="${i#*=}"
    shift
    ;;

    # Type
    -t=*)
    TYPE="${i#*=}"
    shift
    ;;

    # Unknown
    *)
    ;;
esac
done

USER=""
if [ "${TYPE}" == "user" ]; then
    USER="biocompute-DM_user_${NAME}"
else
    USER="biocompute-DM_group_${NAME}"
fi

# Delete the user
deluser "${USER}"
exit
