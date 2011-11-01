/********************************************
! Copyright Li Jiahui@ifts <jiahuili@zju.edu.cn>
!
! This code is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details:
! <http://www.gnu.org/licenses/>.
*********************************************/

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

void createdir_(char *path);
void rmdir_(char *path);

void createdir_(char *path)
{
//  char pathname[4];
  int iflag=0;
//  sprintf(pathname, "%04d",1);
  if(access(path,F_OK)<0) iflag=mkdir(path,0766);
//  return iflag;
}

void rmdir_(char *path)
{
  remove(path);
}
