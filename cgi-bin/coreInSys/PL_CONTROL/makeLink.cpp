#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <string>

char * targetDir;
char * linkDir;

int main(int argc, char** argv)
{
	targetDir=argv[1];
	linkDir=argv[2];

	std::string const command = std::string("ssh root@159.92.115.5 'ln -s ")+argv[1] +" "+argv[2]+"'";
	setuid( 0 );
	system( command.c_str() );

	return 0;
}
