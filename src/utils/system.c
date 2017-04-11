#include <getopt.h>
#include <limits.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

int system_c(const char* input)
{
  return system(input);
}

bool pwd_c(const int* cwd_size, char* cwd)
{
  char buffer[*cwd_size + 1];
  char * result = getcwd(buffer, *cwd_size + 1);
  if (result==NULL)
  {
    return false;
  }
  else
  {
    strcpy(cwd,result);
    return true;
  }
}

bool get_terminal_width_c(int * result)
{
}

bool get_flag_c
(
  const int  * argc,
        char * argvs,
  const char * options,
        char * flag,
        char * output
)
{
  char * argv[*argc+1];
  argv[0] = argvs;
  
  for (int i=1;i<*argc;++i)
  {
    argv[i] = argv[i-1] + strlen(argv[i-1]) + 1;
  }
  argv[*argc] = NULL;
  
  const int result = getopt(*argc, argv, options);
  if (result==-1)
  {
    return false;
  }
  else
  {
    *flag = result;
    if (optarg==NULL)
    {
      *output = '\0';
    }
    else
    {
      strcpy(output,optarg);
    }
    return true;
  }
}
