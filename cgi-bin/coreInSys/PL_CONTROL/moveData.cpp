#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <string>

char * sourceDir;
char * destDir;

int main(int argc, char** argv)
{
	sourceDir=argv[1];
	destDir=argv[2];

	std::string const command = std::string("ssh root@nsd01 'mv ")+argv[1] +" "+argv[2]+"'";
	setuid( 0 );
	system( command.c_str() );

	return 0;
}
