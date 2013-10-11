/*
 * Error handling
 * Copyright (C) 2011 Wray Buntine 
 * All rights reserved.
 *
 * The contents of this file are subject to the Mozilla Public License
 * Version 1.1 (the "License"); you may not use this file except in
 * compliance with the License. You may obtain a copy of the License at
 * http://www.mozilla.org/MPL/
 *
 * Software distributed under the License is distributed on an "AS IS"
 * basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
 * License for the specific language governing rights and limitations
 * under the License.
 *
 * Author: Wray Buntine (wray.buntine@nicta.com.au)
 *     
 */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "yaps.h"

/*
 *    setup so user can define their own yapper
 */
static void (*yaps)(const char *format, va_list ap) = NULL;
static void call_yaps(const char *fmt, ...) {
  va_list ap;
  if ( yaps==NULL )
    yaps_quit("Cannot call_yaps on NULL yapper\n");
  va_start(ap, fmt);
  yaps(fmt, ap);
  va_end(ap);
}
void yaps_yapper(void (*yapper)(const char *format, va_list ap)) {
  yaps = yapper;
}

/*
 * Nonfatal message 
 */
void yaps_message(const char *fmt, ...)
{
  va_list ap;
  va_start(ap, fmt);
  if ( yaps ) 
    yaps(fmt, ap);
  else
    vfprintf(stderr, fmt, ap);
  va_end(ap);
}
/*
 * Fatal message 
 */
void yaps_quit(const char *fmt, ...)
{
  va_list ap;
  va_start(ap, fmt);
  if ( yaps ) 
    yaps(fmt, ap);
  else
    vfprintf(stderr, fmt, ap);
  va_end(ap);
  exit(1);
}
/*
 * Fatal error related to a system call.
 */
void yaps_sysquit(const char *fmt, ...)
{
  va_list ap;
  if ( yaps ) 
    call_yaps("%s: ", strerror(errno));
  else
    fprintf(stderr, "%s: ", strerror(errno));
  va_start(ap, fmt);
  if ( yaps ) 
    yaps(fmt, ap);
  else
    vfprintf(stderr, fmt, ap);
  va_end(ap);
  exit(1);
}
