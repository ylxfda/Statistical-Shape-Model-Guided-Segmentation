/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: ClickedPointEvent.h,v $
  Language:  C++
  Date:      $Date: 2004/07/05 21:37:24 $
  Version:   $Revision: 1.1 $

  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkClickedPointEvent_h__
#define __itkClickedPointEvent_h__

#include "itkEventObject.h"     

/** 
  \class ClickedPointEvent
  \brief Class for defining ClickedPointEvent 

  This defines a ClickedEventEvent class which is derived from the
  itk::AnyEvent class.

*/

itkEventMacro( ClickedPointEvent, itk::AnyEvent );


#endif

