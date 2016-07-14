#include <stdio.h>
#include <stdlib.h>

extern void model2(int argc, char **argv);
extern void model2(int argc, char **argv);

int main(int argc, char **argv)
{
  model1(argc, argv);
  model2(argc, argv);
  return 0;
}
