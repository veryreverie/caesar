#include <limits.h>
#include <stdbool.h>
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
