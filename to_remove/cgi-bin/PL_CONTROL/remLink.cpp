#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <string>

char * linkName;

int main(int argc, char** argv)
{
	linkName=argv[1];

	std::string const command = std::string("ssh root@159.92.115.5 'rm  ")+argv[1] +"'";
	setuid( 0 );
	system( command.c_str() );

	return 0;
}