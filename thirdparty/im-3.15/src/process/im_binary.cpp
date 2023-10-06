/** \file
 * \brief Operations for Binary Images
 *
 * See Copyright Notice in im_lib.h
 */


#include <im.h>
#include <im_util.h>

#include "im_process_counter.h"
#include "im_process_loc.h"
#include "im_process_pnt.h"

#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <string.h>
#include <math.h>


/* Direction masks:      */
/*   N     S   W     E    */
static int masks[] = { 0200, 0002, 0040, 0010 };

/*  True if pixel neighbor map indicates the pixel is 8-simple and  */
/*  not an end point and thus can be deleted.  The neighborhood  */
/*  map is defined as an integer of bits abcdefghi with a non-zero  */
/*  bit representing a non-zero pixel.  The bit assignment for the  */
/*  neighborhood is:            */
/*                  */
/*        a b c          */
/*        d e f          */
/*        g h i          */

static unsigned char isdelete[512] = 
{
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
};

static int DoThinImage(imbyte *map, int xsize, int ysize, int counter)
{
  int    x, y;    /* Pixel location    */
  int    i;    /* Pass index      */
  int    pc  = 0;  /* Pass count      */
  int    count  = 1;  /* Deleted pixel count    */
  int    p, q;    /* Neighborhood maps of adjacent cells      */
  imbyte    *qb;    /* Neighborhood maps of previous scanline      */
  int    m;    /* Deletion direction mask  */
  
  qb = (imbyte *) malloc(xsize);
  qb[xsize-1] = 0;    /* Used for lower-right pixel  */

  while ( count ) 
  {    
    /* Scan src_image while deletions  */
    pc++;
    count = 0;

    imCounterTotal(counter, ysize+1, "Processing... (undef.)");

    for ( i = 0 ; i < 4 ; i++ ) 
    {
      m = masks[i];
      
      /* Build initial previous scan buffer.      */
      
      p = map[0] != 0;
      for (x = 0 ; x < xsize-1 ; x++)
      {
        p = ((p<<1)&0006) | (map[x+1] != 0);
        qb[x] = (imbyte)p;
      }

      if (!imCounterInc(counter)) 
      {
        free(qb);
        return 0;
      }

      /* Scan src_image for pixel deletion candidates.    */
      
      for ( y = 0 ; y < ysize-1 ; y++ ) 
      {
        q = qb[0];
        p = ((q<<3)&0110) | (map[(y+1)*xsize] != 0);
        
        for ( x = 0 ; x < xsize-1 ; x++ ) 
        {
          q = qb[x];
          p = ((p<<1)&0666) | ((q<<3)&0110) | (map[(y+1)*xsize + x+1] != 0);
          qb[x] = (imbyte)p;

          if  (((p&m) == 0) && isdelete[p] ) 
          {
            count++;
            map[y*xsize + x] = 0;
          }
        }
        
        /* Process right edge pixel.      */
       
        p = (p<<1)&0666;
        if  ( (p&m) == 0 && isdelete[p] ) 
        {
          count++;
          map[y*xsize + xsize-1] = 0;
        }

        if (!imCounterInc(counter))
        {
          free(qb);
          return 0;
        }
      }
      
      /* Process bottom scan line.        */
      
      for ( x = 0 ; x < xsize ; x++ ) 
      {
        q = qb[x];
        p = ((p<<1)&0666) | ((q<<3)&0110);

        if  ( (p&m) == 0 && isdelete[p] ) 
        {
          count++;
          map[(ysize-1)*xsize + x] = 0;
        }
      }

      if (!imCounterInc(counter))
      {
        free(qb);
        return 0;
      }
    }
  }
  
  free (qb);
  return 1;
}

int imProcessBinThinNhMaps(const imImage* src_image, imImage* dst_image)
{
  int counter = imCounterBegin("BinThinNhMaps");
  imImageCopyData(src_image, dst_image);
  int ret = DoThinImage((imbyte*)dst_image->data[0], dst_image->width, dst_image->height, counter);
  imCounterEnd(counter);
  return ret;
}


/*****************************************************************************/


static inline int zsCountNeighbours(unsigned char* map, int w, int y, int x)
{
  int i, j, count = 0;

  for (i = -1; i <= 1; i++) 
  {
    for (j = -1; j <= 1; j++) 
    {
      if (i != 0 || j != 0)
        count += (map[(y + i)*w + x + j] == 1);
    }
  }

  return count;
}

static inline int zsCountTransitions(unsigned char* map, int w, int y, int x)
{
  return 	((map[(y-1)*w + x] == 0 && map[(y-1)*w + x+1] == 1)
           + (map[(y-1)*w + x+1] == 0 && map[(y)*w + x+1] == 1)
           + (map[(y)*w + x+1] == 0 && map[(y+1)*w + x+1] == 1)
           + (map[(y+1)*w + x+1] == 0 && map[(y+1)*w + x] == 1)
           + (map[(y+1)*w + x] == 0 && map[(y+1)*w + x-1] == 1)
           + (map[(y+1)*w + x-1] == 0 && map[(y)*w + x-1] == 1)
           + (map[(y)*w + x-1] == 0 && map[(y-1)*w + x-1] == 1)
           + (map[(y-1)*w + x-1] == 0 && map[(y-1)*w + x] == 1));
}

static inline int zhangSuenTest1(unsigned char* map, int w, int y, int x)
{
  int neighbours = zsCountNeighbours(map, w, y, x);

  return ((neighbours >= 2 && neighbours <= 6)
          && (zsCountTransitions(map, w, y, x) == 1)
          && (map[(y-1)*w + x] == 0 || map[(y)*w + x+1] == 0 || map[(y+1)*w + x] == 0)
          && (map[(y)*w + x+1] == 0 || map[(y+1)*w + x] == 0 || map[(y)*w + x-1] == 0));
}

static inline int zhangSuenTest2(unsigned char* map, int w, int y, int x)
{
  int neighbours = zsCountNeighbours(map, w, y, x);

  return ((neighbours >= 2 && neighbours <= 6)
          && (zsCountTransitions(map, w, y, x) == 1)
          && (map[(y-1)*w + x] == 0 || map[(y)*w + x+1] == 0 || map[(y)*w + x-1] == 0)
          && (map[(y-1)*w + x] == 0 || map[(y+1)*w + x] == 0 || map[(y)*w + x+1] == 0));
}

int imProcessBinThinZhangSuen(imImage* src_image, imImage* dst_image)
{
  int start_y = 1, start_x = 1, end_y, end_x, i, marker_count, h, w, processed, markers_size;
  unsigned char* src_map = (unsigned char*)src_image->data[0];
  unsigned char* map = (unsigned char*)dst_image->data[0];
  int* markers, x, y;

  w = src_image->width;
  h = src_image->height;
  end_y = h - 2;
  end_x = w - 2;

  int counter = imCounterBegin("BinThinNhMaps");

  if (src_map != map)
    memcpy(map, src_map, src_image->count);

  markers_size = (end_y - start_y+1)*(end_x - start_x+1) * sizeof(int);
  markers = (int*)malloc(markers_size);

  do 
  {
    marker_count = 0;

    imCounterTotal(counter, 2 * (end_y - start_y + 1), "Processing... (undef.)");

    for (y = start_y; y <= end_y; y++) 
    {
      for (x = start_x; x <= end_x; x++) 
      {
        i = y*w + x;
        if (map[i] == 1 && zhangSuenTest1(map, w, y, x) == 1)
        {
          markers[marker_count] = i;
          marker_count++;
        }
      }

      if (!imCounterInc(counter))
      {
        free(markers);
        imCounterEnd(counter);
        return 0;
      }
    }

    processed = (marker_count > 0);

    for (i = 0; i < marker_count; i++) 
    {
      map[markers[i]] = 0;
    }

    marker_count = 0;

    for (y = start_y; y <= end_y; y++) 
    {
      for (x = start_x; x <= end_x; x++) 
      {
        i = y*w + x;
        if (map[i] == 1 && zhangSuenTest2(map, w, y, x) == 1) 
        {
          markers[marker_count] = i;
          marker_count++;
        }
      }

      if (!imCounterInc(counter))
      {
        free(markers);
        imCounterEnd(counter);
        return 0;
      }
    }

    if (processed == 0)
      processed = (marker_count > 0);

    for (i = 0; i < marker_count; i++) 
    {
      map[markers[i]] = 0;
    }

  } while (processed == 1);

  free(markers);
  imCounterEnd(counter);
  return 1;
}
