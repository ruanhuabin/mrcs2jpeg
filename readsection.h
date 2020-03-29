/*******************************************************************
 *       Filename:  readsection.h                                     
 *                                                                 
 *    Description:                                         
 *                                                                 
 *        Version:  1.0                                            
 *        Created:  2020年03月29日 11时03分52秒                                 
 *       Revision:  none                                           
 *       Compiler:  gcc                                           
 *                                                                 
 *         Author:  Ruan Huabin                                      
 *          Email:  ruanhuabin@tsinghua.edu.cn                                        
 *        Company:  Dep. of CS, Tsinghua Unversity                                      
 *                                                                 
 *******************************************************************/

#ifndef  READSECTION_H
#define  READSECTION_H
#include "mrc.h"
bool readMRCSection(float **im, int &dimx, int &dimy, const char* filename, int axis, int slice);
bool readMRCSection2(float **im, int &dimx, int &dimy, MRC &mrc, int axis, int slice);
#endif

