#define _POSIX_C_SOURCE 200809L

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
  FILE * output_file = popen("tput cols","r");
  if (!output_file)
  {
    return false;
  }
  const int success2 = fscanf(output_file, "%i", result);
  pclose(output_file);
  if (success2 != 1)
  {
    return false;
  }
  return true;
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
  
  opterr = 0;
  const int result = getopt(*argc, argv, options);
  if (result==-1)
  {
    return false;
  }
  else
  {
    if (result==1)
    {
      *flag = ' ';
    }
    else
    {
      *flag = result;
    }
    
    if (*flag=='?' || *flag==':')
    {
      *output = optopt;
      *(output+1) = '\0';
    }
    else if (optarg==NULL)
    {
      strcpy(output,"");
    }
    else
    {
      strcpy(output,optarg);
    }
    return true;
  }
}
