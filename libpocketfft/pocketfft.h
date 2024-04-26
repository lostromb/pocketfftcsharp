/*
 * This file is part of pocketfft.
 * Licensed under a 3-clause BSD style license - see LICENSE.md
 */

 /*! \file pocketfft.h
  *  Public interface of the pocketfft library
  *
  *  Copyright (C) 2008-2018 Max-Planck-Society
  *  \author Martin Reinecke
  */

#ifndef POCKETFFT_H
#define POCKETFFT_H

#include <stdlib.h>

struct cfft_plan_i;
typedef struct cfft_plan_i* cfft_plan;
__declspec(dllexport) cfft_plan make_cfft_plan(int length);
__declspec(dllexport) void destroy_cfft_plan(cfft_plan plan);
__declspec(dllexport) int cfft_backward(cfft_plan plan, double c[], double fct);
__declspec(dllexport) int cfft_forward(cfft_plan plan, double c[], double fct);
__declspec(dllexport) int cfft_length(cfft_plan plan);

struct rfft_plan_i;
typedef struct rfft_plan_i* rfft_plan;
__declspec(dllexport) rfft_plan make_rfft_plan(int length);
__declspec(dllexport) void destroy_rfft_plan(rfft_plan plan);
__declspec(dllexport) int rfft_backward(rfft_plan plan, double c[], double fct);
__declspec(dllexport) int rfft_forward(rfft_plan plan, double c[], double fct);
__declspec(dllexport) int rfft_length(rfft_plan plan);

#endif
