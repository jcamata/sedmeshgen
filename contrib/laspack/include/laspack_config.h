/****************************************************************************/
/*                                 elcmp.h                                  */
/****************************************************************************/
/*                                                                          */
/* includes definitions of libMesh                                          */
/*                                                                          */
/* Copyright (C) 1992-1996 Tomas Skalicky. All rights reserved.             */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*        ANY USE OF THIS CODE CONSTITUTES ACCEPTANCE OF THE TERMS          */
/*              OF THE COPYRIGHT NOTICE (SEE FILE COPYRGHT.H)               */
/*                                                                          */
/****************************************************************************/

#ifndef LASPACK_CONFIG_H
#define LASPACK_CONFIG_H

#ifdef __cplusplus
      /* someone included us from C++ */
#define _LP_INCLUDED_FROM_CPLUSPLUS 1
#else
/* compile LASPACK with real arithmetic in C */
#undef _LP_INCLUDED_FROM_CPLUSPLUS
#  endif /* __cplusplus */


#endif /* LASPACK_CONFIG_H */
