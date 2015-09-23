#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <string>

char * sftp_user;

int main(int argc, char** argv)
{
	sftp_user=argv[1];

	std::string const command=std::string("/home/biocis/www/cgi-bin/PL_CONTROL/userDel.sh ")+sftp_user;
	setuid(0);
	system(command.c_str());

	return 0;
}
