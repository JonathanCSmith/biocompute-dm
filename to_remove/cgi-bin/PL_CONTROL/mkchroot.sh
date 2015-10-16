#!/bin/sh




chroot="${1}"

if [ "${chroot}" == "" ]; then
    echo "FATAL: I need a location to create the chroot!"
    exit 1
fi

#if [ -e ${chroot} ]; then
#    echo "FATAL: ${chroot} already exists!"
#    exit 1
#fi

mkdir -p ${chroot}/{usr/bin,dev,etc,lib64,usr/lib64,usr/libexec/openssh,usr/libexec}

for bin in /usr/bin/scp /usr/bin/rssh /usr/libexec/rssh_chroot_helper /usr/libexec/openssh/sftp-server ;

do

    cp ${bin} ${chroot}${bin}
 

    for lib in `ldd ${bin} | awk '{print $3}'`;
    do
        if [ -f ${lib} ]; then
            cp ${lib} ${chroot}/${lib}
        fi
    done

done

#cp /usr/libexec/rssh_chroot_helper ${chroot}/usr/libexec/
cp /lib64/ld-linux-x86-64.so.2 ${chroot}/lib64/
cp /lib64/libc.so.6 ${chroot}/lib64/
cp /lib64/libnss_compat.so.2 ${chroot}/lib64/
cp /lib64/libcrypt.so.1 ${chroot}/lib64/

#cp /lib/libcrypt.so.1 ${chroot}/lib/
#cp /lib/libnss_compat.so.2 ${chroot}/lib/

mknod -m 0666 ${chroot}/dev/null c 1 3

