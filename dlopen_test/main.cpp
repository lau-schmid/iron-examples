#include <iostream>
#include <dlfcn.h>

using namespace std;

int main(int argc, char *argv[])
{
   if (argc == 1)
     exit(-1);

   cout << "start" << endl;

   

   void *mHandle = dlopen(argv[1], RTLD_LOCAL | RTLD_LAZY);

   if (mHandle == NULL)
   {
     cout << "handle is NULL!" << endl;
     exit(0);
   }
   else
     cout << "handle is not NULL" << endl;
   
   void (*func)();
   func = (void (*)()) dlsym(mHandle, "func");
   
   cout << "got func pointer" << endl;

   func();

   cout << "end" << endl;
}

