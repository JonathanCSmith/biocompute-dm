#!/bin/sh

username=$1

userdel -r $username

# remove the line with the username from the /etc/rssh.conf file
sed -i".bak" '/'$username'/d' /etc/rssh.conf

# remove the chroot part
rm -rf /home/transfer/$username



