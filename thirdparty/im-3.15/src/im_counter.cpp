/** \file
 * \brief Processing Counter
 *
 * See Copyright Notice in im_lib.h
 */

#include "im_counter.h"

#include <stdlib.h>
#include <memory.h>


static imCounterCallback iCounterFunc = NULL;
static void* iCounterUserData = NULL;

imCounterCallback imCounterSetCallback(void* user_data, imCounterCallback counter_func)
{
  imCounterCallback old_counter_func = iCounterFunc;
  iCounterFunc = counter_func;
  if (user_data)
    iCounterUserData = user_data;
  return old_counter_func;
}

int imCounterHasCallback(void)
{
  return iCounterFunc!=NULL;
}

struct iCounter
{
  int total;
  int current;
  int has_begin;
  const char* message;
  void* userdata;
};

#define MAX_COUNTERS 50
static iCounter iCounterList[MAX_COUNTERS];

int imCounterBegin(const char* title)
{
  static int first = 1;
  if (first)
  {
    memset(iCounterList, 0, MAX_COUNTERS*sizeof(iCounter));
    first = 0;
  }

  if (!iCounterFunc) 
    return -1;             // counter management is useless

  int counter = -1;
  for (int i = 0; i < MAX_COUNTERS; i++)
  {
    if (iCounterList[i].has_begin == 0)  // the counter is free
    {
      counter = i;
      break;
    }
  }

  if (counter == -1) 
    return -1;             // too many counters

  iCounter *ct = &iCounterList[counter];

  ct->has_begin = 1;

  iCounterFunc(counter, iCounterUserData, title, -1);

  return counter;
}

void imCounterEnd(int counter)
{
  if (counter < 0 || counter >= MAX_COUNTERS || !iCounterFunc)  // invalid counter
    return;

  iCounter *ct = &iCounterList[counter];

  if (ct->has_begin == 0) // counter with no begin
    return;

  iCounterFunc(counter, iCounterUserData, NULL, 1001);
  memset(ct, 0, sizeof(iCounter));
}

void* imCounterGetUserData(int counter)
{
  if (counter < 0 || counter >= MAX_COUNTERS || !iCounterFunc)  // invalid counter
    return NULL;

  iCounter *ct = &iCounterList[counter];

  return ct->userdata;
}

void imCounterSetUserData(int counter, void* userdata)
{
  if (counter < 0 || counter >= MAX_COUNTERS || !iCounterFunc)  // invalid counter
    return;

  iCounter *ct = &iCounterList[counter];

  ct->userdata = userdata;
}

int imCounterInc(int counter)
{
  if (counter < 0 || counter >= MAX_COUNTERS || !iCounterFunc)  // invalid counter
    return 1;

  iCounter *ct = &iCounterList[counter];

  if (ct->has_begin == 0 ||  // counter with no begin or no total
      ct->total == 0)
    return 1;

  const char* msg = NULL;
  if (ct->current == 0)
    msg = ct->message;

  ct->current++;

  int progress = (int)((ct->current * 1000.0f)/ct->total);

  if (ct->current == ct->total)
    ct->current = 0;

  return iCounterFunc(counter, iCounterUserData, msg, progress);
}

int imCounterIncTo(int counter, int count)
{
  if (counter < 0 || counter >= MAX_COUNTERS || !iCounterFunc)  // invalid counter
    return 1;

  iCounter *ct = &iCounterList[counter];

  if (ct->has_begin == 0 ||  // counter with no begin or no total
      ct->total == 0)
    return 1;

  if (count <= 0) count = 0;
  if (count >= ct->total) count = ct->total;

  ct->current = count;

  const char* msg = NULL;
  if (ct->current == 0)
    msg = ct->message;

  int progress = (int)((ct->current * 1000.0f)/ct->total);

  if (ct->current == ct->total)
    ct->current = 0;

  return iCounterFunc(counter, iCounterUserData, msg, progress);
}

void imCounterTotal(int counter, int total, const char* message)
{
  if (counter < 0 || counter >= MAX_COUNTERS || !iCounterFunc)  // invalid counter
    return;

  iCounter *ct = &iCounterList[counter];

  if (ct->has_begin == 0)  // counter with no begin 
    return;

  ct->message = message;
  ct->total = total;
  ct->current = 0;
}
