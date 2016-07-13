#include <stdio.h>
#include <stdlib.h>

extern void userdef(int argc, char **argv);
extern void model2(int argc, char **argv);

int main(int argc, char **argv)
{
  userdef(argc, argv);
  model2(argc, argv);
  return 0;
}
