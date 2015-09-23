#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <string>

char * sftp_user;
char * sftp_pass;
char * home_dir;

int main(int argc, char** argv)
{
	sftp_user=argv[1];
	sftp_pass=argv[2];
	home_dir=argv[3];

	std::string const command = std::string( "/home/biocis/www/cgi-bin/coreInSys/PL_CONTROL/sftp_setup.sh " ) + argv[1] +" "+argv[2]+" "+argv[3];
	setuid( 0 );
	system( command.c_str() );

	return 0;
}
