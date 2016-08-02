/*
 * agro_lib.h
 * 
 * Copyright 2016 Andriy Gonda <andriy.gonda@gmail.com>
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

#ifndef AGRO_LIB_H
#define AGRO_LIB_H
#ifdef __cplusplus
extern "C"
{
#endif
#include "simple_geofence.h"

typedef struct
 {
  long double longtitude;
  long double latitude;
  unsigned char marker;
 } POINT;


/*Checking point on belonging to geofence*/
extern unsigned point_in_geofence(size_t npol , POINT (*p)[], long double x, long double y);

/*Distance from the pair of points on the grid point*/
extern double  dist_to_point(POINT A,POINT B,POINT C);

/*converting coordinates to x ,y system of coordinates in meters*/
extern int convert_coordinates(long double *longt,long  double *lat);

/*square of geofence with overlay grid ( grid size in meters) */
extern double size_area(POINT  (*sourse_geofence)[],size_t N_geofence, double grid_size);

/*analyse squares of processed, not processed and poly processed ground fields by the GPS-track */
extern void area_analyse(POINT (*geofence)[],size_t N_geofence,POINT (*area)[],size_t N_area,double area_size,double width_size,double grid_size,double *used_size,double *not_used_size, double *poly_used_size);

/*get coordinates latitude and longtitude with the source line*/
extern void get_coordinates( char *line, long double *longtitude, long double *latitude);

/*get a gps-track particle that included in geofence*/
extern void get_area(const char *track_path,const char *geofence_path, const char *area_path); 

/*Reducing points in source geofence points file*/
extern void reducing_points(const char *src_path, const char *dst_path,const double reducing_k);

/*Calculation square of geofence by overlay grid with .fld-geofence file*/
extern double field_size(char *field_path , double grid_size);

/*Analyse track with .tmsg - file and geofence with .fld - file and return results of calculations*/
extern void area_analyser(char *field_path, char *area_path,double area_size, double width, double  grid_size, double *used_size, double *not_used_size, double *poly_used_size);

extern void format_data(char *);
extern double side(POINT a, POINT b);
extern double determinant(POINT p1, POINT p2);

/*Finding perimeter of geofence with .fld file */
extern double field_perimeter(char *fld_sourse);

/*Finding polygon size with array of points*/
extern double poly_size(POINT (*poly)[], size_t N_poly);

/*Convert json field format to field format .fld*/
extern void json_fields_to_fld(char *sourse_filename, int field_id,char *dst_filename);

/* Convert json array messages to .tmsg format file*/
extern void json_messages_to_tmsg(char *sourse_filename,char *destination_filename);

/*Get field perimeter in kilometers and field size in hectares with json field file and id of field by grids algoritm*/
extern void field_size_and_perimeter(char *json_field_file,int field_id ,double grid_size,double *area_size, double *perimeter);

/*Analysis of the area processed, unprocessed and re-processed land by gps-track and grid with json files*/
extern void json_field_analyse(char *json_field_file,int field_id,char *json_message_file,
              double area_size, double width, double grid_size,double *used_field_size,double *not_used_field_size, double *poly_used_field_size);

/*Get size of field in hectares by polygonals algoritm*/
extern double json_polygon_size(char *json_field_sourse,int field_id);

/*Get field square with json field file and id by grids algoritm.*/
extern double json_field_size(char *json_sourse, int field_id, double grid_size);

/*Calculation of field perimeter by json file*/
extern double json_field_perimeter(char *json_field_sourse, int field_id);

/*Calculation the  square of circle geofence*/
extern double circle_geofence_square(double radius);

/*Calculation the perimeter of circle geofence*/
extern double circle_geofence_perimeter(double radius);

extern double radial_distance(POINT c, POINT p);

/*Checking point on including in circle geofence*/
extern int point_in_circle_geofence(POINT center, POINT current,double radius);

/* Analyse track in circle geofence with .tmsg-file */
extern int circle_geofence_analyser_tmsg(double center_longtitude,double center_latitude, double radius,
            double width , double grid_size,char *track_path, double *processed, double *not_processed, double *poly_processed);
/*Analyse track in circle geofence with .json-file*/
extern int json_circle_geofence_analyser(double center_longtitude, double center_latitude,double radius, double width ,
              double grid_size, char *json_path, double *processed, double *not_processed , double *poly_processed);

/*Checking geofence on simplicity with json file*/
extern int json_is_simple_geofence(char *json_source_path , int id);

#ifdef __cplusplus
extern "C"
}
#endif

#endif
