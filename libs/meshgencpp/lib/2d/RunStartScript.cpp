#include <cstdlib>
#include <cstdio>
#include "assert.h"
using namespace std;

void RunStartScript(char*& outputdir)
{
  // run startscript
  // to create output directory if it does not exist
  // and copy data files to output directory
  char command_str[1024];
  int numchars = snprintf(command_str,1024,
			  "if test -f startscript && test -x startscript;\n"
			  "then ./startscript %s %d\n"
			  "else ${MESHGENCPP}/lib/2d/startscript %s %d\n"
			  "fi",outputdir,2,outputdir,2);
  assert_lt(numchars,1023);
  assert_gt(numchars,0);
  int exit_status = 0;
  exit_status = system(command_str);
}
