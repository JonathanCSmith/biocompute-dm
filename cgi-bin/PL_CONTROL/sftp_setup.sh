#!/bin/sh

#script to add a chrooted sftp account with a restricted shell

sftp_user=$1
sftp_pass=$2
home_dir=$3 

useradd --home=$home_dir $sftp_user


echo $sftp_pass | passwd $sftp_user --stdin

usermod -s /usr/bin/rssh $sftp_user


usermod --home $home_dir/data $sftp_user



/home/biocis/www/cgi-bin/coreInSys/PL_CONTROL/mkchroot.sh $home_dir

getent passwd $sftp_user > $home_dir/etc/passwd

echo "user="$sftp_user":022:00011:"$home_dir >>/etc/rssh.conf

mkdir $home_dir/data

chown $sftp_user:$sftp_user $home_dir/data

chmod a+rx $home_dir
