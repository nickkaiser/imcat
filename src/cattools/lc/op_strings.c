/*
 * op_strings.c
 */

#include <stdio.h>
#include "../../catlib/cat.h"
#include "stack.h"
#include "error.h"
#include <string.h>


void streqinit(item *operator)
{
  item *result, *op1, *op2;
  int tempres, test;
  double *res;
  char *str1, *str2;

  op2 = pop();
  op1 = pop();
  if ((op1->itype != TEXT_TYPE) || (op2->itype != TEXT_TYPE)){
    error_exit("eq: non-string arguments\n");
  }
  result = operator->next = newitem("STREQ_TEMP", NUM_TYPE, 1, 1);
  res = (double *) calloc(1 , sizeof(double));
  result->addr = (void *) res;
  str1 = *((char **) (op1->addr));
  str2 = *((char **) (op2->addr));
  
 /*   /* actually do the op here */ 
/*    printf("no fault in init\n"); */
/*    test = strlen(str1); */
/*    printf("str1 is %d long\n", test); */
/*    tempres = strcmp(str1, str2); */
/*    printf("is it here?\n"); */
/*    /*tempres = strcmp( (char *)(*(str1)), (char *)(*(str2)) );*/ 
/*   /*   printf("str1 is \"%s\" str2 is \"%s\" \n", *(str1), *(str2)); */ 
/*  /*    printf("tempres is %lf\n", tempres); */ 
/*     if(tempres == 0){ */
/*      *res = 1.0; */
/*    } */
/*    else *res = 0.0; */
 
/*  /*   *res = 1.0; */ 
/*  /*    printf("*res is %lf\n", *res); */ 

}
  
  




void streqdoit(item *operator)
{
  item *op1, *op2;
  char *str1, *str2;
  int tempres;
  double *res;
  
  op2 = pop();
  op1 = pop();
  res = (double *) ((operator->next)->addr);
  str1 = *((char **) (op1->addr));
  str2 = *((char **) (op2->addr));
 
  /* do operation here */
  
  tempres = strcmp(str1, str2);
  if(tempres == 0){
    *res = 1.0;
  }

  else *res = 0.0;
  


}



