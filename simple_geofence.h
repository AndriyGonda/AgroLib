/*
 * simple_geofence.h
 * 
 * Copyright 2016 by Andriy Gonda <andriy.gonda@gmail.com>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */

#ifndef SIMPLE_GEOFENCE_H
#define SIMPLE GEOFENCE_H
/*A main geometric point structure*/
typedef struct 
  {
	double x;
	double y;
  } GPOINT;
  
/*A main  geometric line structure*/
 typedef struct
  {
	double A;
	double B;
	double C;
  } GLINE;

#ifdef __cplusplus
extern "C"
{
#endif

/*
 * Get GPS coordinates with source string in format (longtitude; latitude)
 * name: get_coords_geofence
 * @param char *s  -  source string
 * @param double *longtitude  -  a pointer on variable of longtitude
 * @param double *latitude    -  a  pointer on variable of latitude
 */
extern void get_coords_geofence( char *s,  double *longtitude, double *latitude) ;

/*
 * Checking two double numbers on equality with 1e-308 precision
 * name: double_equal
 * @param double a
 * @param double b
 * @return result of comparation a and b 
 */
extern  int double_equal(double a, double b);

/*
 * Checking two double numbers on majority with 1e-308 precision
 * name: double_more
 * @param double a 
 * @param double b
 * @return result of comparation a and b
 */
extern  int double_more(double a, double b);

/*
 * Checking two double numbers on minority with 1e-308 precision
 * name: double_less
 * @param double a 
 * @param double b
 * @return result of comparation a and b 
 */
extern int double_less(double a, double b);

/*
 * Checking two double numbers on majority or equality with 1e-308 precision
 * name: double_more_equal
 * @param double a 
 * @param double b
 * @return result of comparation a and b
 */
extern int double_more_equal(double a, double b);
 
/*
 * Checking two double numbers on maximum
 * name: double_max
 * @param double a
 * @param double b
 * @return maximum number
 */
extern double double_max(double a, double b);

/*
 * Checking two double numbers on minimum
 * name: double_min
 * @param double a
 * @param double b
 * @return minimum number
 */
extern double double_min(double a, double b);
 
/*
 * Checking two geometric points on equality
 * name: equal_point
 * @param GPOINT A
 * @param GPOINT B
 * @return result of comparation two geometric points
 */
extern int equal_point(GPOINT A, GPOINT B);

/*
 * Calculation distance between two geometric points
 * name: distance
 * @param GPOINT A
 * @param GPOINT B
 * @return distance between A and B
 */
extern double distance(GPOINT A, GPOINT B);
 
/*
 * Converting coordinates two points in equation of line
 * name: point_2_to_line
 * @param const GPOINT v
 * @param const GPOINT w
 * @param GLINE *L 
 */
extern void point_2_to_line(const GPOINT v, const GPOINT w, GLINE *L);

/*
 * Checking crossing two lines
 * name: cross_line
 * @param GLINE L1
 * @param GLINE L2
 * @return 1 or 0 if lines L1 and L2 crossed/not crossed
 */
extern int cross_line(GLINE L1, GLINE L2);

/*
 * Finding the point of intersection of two lines
 * name: line2_point
 * @param GLINE L1
 * @param GLINE L2
 * @param GPOINT *P
 * @return 1 if line crossed
 */
extern void line2_to_point(GLINE L1, GLINE L2, GPOINT *P);

/*
 * Checking the intersection of two lines that contains  segments f and s
 * name: segment_line_cross
 * @param GPOINT fL
 * @param GPOINT fR
 * @param GPOINT sL
 * @param GPOINT sR
 * @return 1 if line crossed
 * 
 */
extern int segment_line_cross(GPOINT fL, GPOINT fR, GPOINT sL, GPOINT sR);

/*
 * Checking the intersection of two segments (f and s)
 * name: segment_cross
 * @param GPOINT fL
 * @param GPOINT fR
 * @param GPOINT sL
 * @param GPOINT sR
 * @return 1 if segments crossed
 */
extern int segment_cross(GPOINT fL, GPOINT fR, GPOINT sL, GPOINT sR);

/*
 * Сhecking point belonging to the segment that created by two points A and B
 * name: at_segment
 * @param GPOINT A
 * @param GPOINT B
 * @param GPOINT P
 * @return 1 if point P at segment , that created by A and B points 
 */
extern int at_segment(GPOINT A, GPOINT B , GPOINT P);

/*
 * Checking polygon at simplicity
 * name: is_polygon_simple
 * @param GPOINT (*point)[]
 * @param size_t N_point
 * @return 1 if polygon is simple
 */
extern int is_polygon_simple(GPOINT (*point)[], size_t N_point);

/*
 * Checking geofence at simplicity
 * name: is_geofence_simple
 * @param char *geofence_path
 * @return 1 if polygon is simple
 */
extern int is_geofence_simple(char *geofence_path);

#ifdef __cplusplus
}
#endif
#endif
