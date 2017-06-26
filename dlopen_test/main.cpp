#include <iostream>
#include <dlfcn.h>
#include <cstdlib>

using namespace std;

int main(int argc, char *argv[])
{
   if (argc == 1)
   { 
     cout << "No argument passed!" << endl;
     exit(-1);
   }

   cout << "start" << endl;

  char *filename = argv[1];
  cout << "filename: \"" << filename << "\"" << endl;

   void *mHandle = dlopen(filename, RTLD_LOCAL | RTLD_LAZY);

   if (mHandle == NULL)
   {
     cout << "handle is NULL!" << endl;
     cout << "dlerror: " << dlerror() << endl;
     exit(0);
   }
   else
   {
     cout << "handle is not NULL" << endl;
   }
   
   const char *function_name = "OC_CellML_RHS_routine";
   
   void (*func)(double,double*,double*,double*,double*);
//   func = (void (*)()) dlsym(mHandle, function_name);

   func = (void (*)(double,double*,double*,double*,double*)) dlsym(mHandle,"OC_CellML_RHS_routine");   


   cout << "got func pointer" << endl;
   if( func == NULL)
   {
     cout << "func is NULL" << endl;
   }
   else
   {
    double d1,d2,d3,d4;
    d1 = d2 = d3 = d4 = 0.0;
    func(1.0,&d1,&d2,&d3,&d4);
   }
   
   int *a = (int *)dlsym(mHandle, "a");
   if(a == NULL)
   {
    cout<<"a is NULL"<<endl;
   }
   else
   {
      cout<<"a = "<<*a<<endl;
   }
   cout << "end" << endl;
}

