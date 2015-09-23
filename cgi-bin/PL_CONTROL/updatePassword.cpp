#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <string>


char * sftp_user;
char * sftp_pass;

int main(int argc, char** argv)
{
	sftp_user=argv[1];
	sftp_pass=argv[2];
	std::string const command=std::string("/home/biocis/www/cgi-bin/coreInSys/PL_CONTROL/updatePass.sh ")+sftp_user+" "+sftp_pass;
	setuid(0);
	system(command.c_str());

	return 0;
}
