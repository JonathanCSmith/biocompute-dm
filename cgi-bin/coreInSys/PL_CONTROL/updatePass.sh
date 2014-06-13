#!/bin/sh

#change the password of sftp chrooted with a restricted shell account

sftp_user=$1
sftp_pass=$2


echo $sftp_pass | passwd $sftp_user --stdin

#getent passwd $sftp_user > $home_dir/etc/passwd


