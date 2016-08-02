#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "cJSON.h"
#include "agro_lib.h"

 unsigned point_in_geofence(size_t npol , POINT (*p)[], long double x, long double y)
   {
    unsigned c=0;
    size_t i=0,j=npol-1; 
      for(i=0; i<npol;i++)
        {
          if(( ( ((*p)[i].latitude<=y)&&(y<(*p)[j].latitude) ) || ( ((*p)[j].latitude<=y)&&(y<(*p)[i].latitude) ))&&
            (x>((*p)[j].longtitude-(*p)[i].longtitude)*(y-(*p)[i].latitude)/((*p)[j].latitude-(*p)[i].latitude) + (*p)[i].longtitude) )
              c=!c;
            j=i;
        }
     return c;
    }

int convert_coordinates( long double *longt, long  double *lat)
  {
	 
	 double longtitude = *longt;
	 double latitude   =  *lat;
	 int n = (int)((6+longtitude)/6);
	  double l = ( longtitude  - (3+ 6*(n-1)) ) / 57.29577951;
	  latitude*=M_PI/180.;
	  if(!longt||!lat) return 0;
	  double sin_lat = sin(latitude);
	  double square_sin = sin_lat*sin_lat;
	  double quadro_sin = square_sin*square_sin;
	  double sixty_sin  = quadro_sin*square_sin;
	  *longt = 6367558.4968*latitude - sin(2*latitude)*(16002.8900 +  66.9607*square_sin + 0.3515*quadro_sin -
	       l*l*( 1594561.25+ 5336.535*square_sin + 29.790*quadro_sin  + 0.149*sixty_sin +
	       l*l*(672483.4 - 811219.9*square_sin +  5420.0*quadro_sin  - 10.6*sixty_sin +
	       l*l*(278194. -  830174.*square_sin + 572434.*quadro_sin  - 16010.*sixty_sin +
	       l*l*(109500. - 574700.*square_sin +  863700.*quadro_sin - 398600.*sixty_sin ) ) ) ) );
	       
	   *lat = (5+10*n)*1e+5 + l*cos(latitude)*(6378245 + 21346.1415*square_sin + 107.1590*quadro_sin  +
	        0.5977*sixty_sin + l*l*(1070204.16 - 2136826.66*square_sin +17.98*quadro_sin - 11.99*sixty_sin +
	         l*l*(270806.  - 1523417.*square_sin + 1327645.*quadro_sin - 21701.*sixty_sin +
	          l*l*(79690. - 866190.*square_sin + 1730360.*quadro_sin - 945460.*sixty_sin ) ) ) ); 
	 return 1;
  }
   
double size_area( POINT (*sourse_geofence)[],size_t N_geofence, double grid_size)
  {
      size_t i,j;
      POINT (*geofence)[N_geofence] = malloc(sizeof(*geofence));
      for(size_t i=0; i<N_geofence; i++)
        {
            (*geofence)[i] = (*sourse_geofence)[i];
        }
       #pragma omp parallel for private(i) shared(geofence) schedule(dynamic)
      for(size_t i=0; i<N_geofence; i++)
       {
          convert_coordinates(&(*geofence)[i].longtitude,&(*geofence)[i].latitude);
       }
     
     long double max_longtitude=(*geofence)[0].longtitude,
                 min_longtitude=(*geofence)[0].longtitude,
                 max_latitude=(*geofence)[0].latitude,
                 min_latitude=(*geofence)[0].latitude;
               
         for(size_t i=0; i<N_geofence; i++)
           {
               if(max_longtitude<(*geofence)[i].longtitude) max_longtitude=(*geofence)[i].longtitude;
               if(min_longtitude>(*geofence)[i].longtitude) min_longtitude=(*geofence)[i].longtitude;
               if(max_latitude<(*geofence)[i].latitude) max_latitude=(*geofence)[i].latitude;
               if(min_latitude>(*geofence)[i].latitude) min_latitude=(*geofence)[i].latitude;
           }
           
          size_t N =0;
          size_t N_TOTAL =0;
          
          size_t Npx = (size_t)((max_longtitude - min_longtitude)/grid_size);
          size_t Npy = (size_t)((max_latitude-min_latitude)/grid_size);
         
          omp_set_num_threads(omp_get_num_procs());
         #pragma omp parallel for private(i,j)   reduction(+:N) reduction(+:N_TOTAL) schedule (dynamic) 
          for(i=0; i<=Npx; i++)
             for(j=0; j<=Npy; j++)
                {
                  if(point_in_geofence(N_geofence,geofence,min_longtitude + i*grid_size,min_latitude + j*grid_size)) N++;
                  N_TOTAL++;
                }
               free(geofence);
             double area, area_of_rect =(max_longtitude - min_longtitude)*(max_latitude-min_latitude);
            area =  ((double)N/(double)N_TOTAL)*area_of_rect;
            area/=10000.;
            return area ;
  }
  
  double  dist_to_point(POINT A,POINT B,POINT C)
   {
     POINT AC,AB, BA, BC; //vectors
     //AC vector
     AC.longtitude = C.longtitude - A.longtitude;
     AC.latitude  = C.latitude - A.latitude;
     //BC  vector
     BC.longtitude  =  C.longtitude - B.longtitude;
     BC.latitude   =  C.latitude - B.latitude;
     //AB vector
     AB.longtitude = B.longtitude - A.longtitude;
     AB.latitude   = B.latitude - A.latitude;
     // BA vector
     BA.longtitude =  A.longtitude - B.longtitude;
     BA.latitude   =  A.latitude -  B.latitude;
     //
     long double a ,b,c,p;
      a = sqrt(AB.longtitude*AB.longtitude + AB.latitude*AB.latitude);
      b =sqrt(BC.longtitude*BC.longtitude + BC.latitude*BC.latitude);
      c = sqrt(AC.longtitude*AC.longtitude + AC.latitude*AC.latitude);
      p = (a+b+c)*0.5;
      long double h = (p*(p-a)*(p-b)*(p-c));
      double  AC_AB = (AC.longtitude*AB.longtitude + AC.latitude*AB.latitude);
      double  BC_BA = (BC.longtitude*BA.longtitude + BC.latitude*BA.latitude);
      if(AC_AB<0||BC_BA<0)
      return 999999.;
      return sqrt(h)*(2./a); 
 }
 
 void area_analyse(POINT (*sourse_geofence)[],size_t N_geofence,POINT (*sourse_area)[],size_t N_area,double area_size,double width_size,double grid_size,double *used_size,double *not_used_size, double *poly_used_size)
  {
          POINT (*geofence)[N_geofence]=malloc(sizeof(*geofence));
          for(size_t i=0; i<N_geofence;i++) 
                       (*geofence)[i]=(*sourse_geofence)[i];
           
          POINT (*area)[N_area]=malloc(sizeof(*area));
          for(size_t i=0; i<N_area;i++) 
                          (*area)[i]=(*sourse_area)[i];
           size_t i;
          #pragma omp parallel for private(i) shared(area) schedule(dynamic)               
          for( i=0; i<N_area;i++)
           {
			   convert_coordinates(&(*area)[i].longtitude,&(*area)[i].latitude);
               (*area)[i].marker=0;
           }
           
           #pragma omp parallel for private(i) shared(geofence) schedule(dynamic)
           for(i=0;i<N_geofence; i++)
            {
				convert_coordinates(&(*geofence)[i].longtitude,&(*geofence)[i].latitude);
            }
            
             long double max_longtitude=(*geofence)[0].longtitude,
                         min_longtitude=(*geofence)[0].longtitude,
                         max_latitude=(*geofence)[0].latitude,
                         min_latitude=(*geofence)[0].latitude;
                         
             for( size_t j=0; j<N_geofence; j++)
              {
               if(max_longtitude<(*geofence)[j].longtitude) max_longtitude=(*geofence)[j].longtitude;
               if(min_longtitude>(*geofence)[j].longtitude) min_longtitude=(*geofence)[j].longtitude;
               if(max_latitude<(*geofence)[j].latitude) max_latitude=(*geofence)[j].latitude;
               if(min_latitude>(*geofence)[j].latitude) min_latitude=(*geofence)[j].latitude;
              }
           
           size_t N_grid_points=0;
           size_t Npx,Npy;
                  Npx = (size_t)((max_longtitude-min_longtitude)/grid_size)+2;
                 Npy = (size_t)((max_latitude - min_latitude)/grid_size)+2;
                 
          size_t j;
           omp_set_num_threads(omp_get_num_procs());
          #pragma omp parallel for private(i,j)   reduction(+:N_grid_points) schedule (dynamic)  
           for( i=0; i<=Npx; i++)
            for(j=0; j<=Npy; j++)
                 {
                    if(point_in_geofence(N_geofence,geofence,min_longtitude+i*grid_size,min_latitude+j*grid_size))
                    N_grid_points++; 
                 }
                
            POINT (*grid_point)[N_grid_points] = malloc(sizeof(*grid_point));
            
          size_t k=0;
          for( i=0; i<=Npx;i++)
          {
              for( j=0; j<=Npy; j++)
                 {
                      if(point_in_geofence(N_geofence,geofence,min_longtitude+i*grid_size,min_latitude+j*grid_size))
                      {
                        (*grid_point)[k].longtitude =min_longtitude + i*grid_size;
                        (*grid_point)[k].latitude = min_latitude + j*grid_size;
                        (*grid_point)[k].marker=0;
                        k++;
                      }
                 }
          }
          free(geofence);
        width_size*=0.5;
    #pragma omp parallel for private(i,j)  shared(grid_point,area)  schedule (dynamic) 
        for(i=0; i<N_grid_points; ++i)
             {
                 for(j=0;j<N_area; ++j)
                   {
                       if(dist_to_point((*area)[j],(*area)[j+1],(*grid_point)[i])<=width_size
                       &&(*grid_point)[i].marker<2)
                           (*grid_point)[i].marker++; 
                   }
             } 
             
     size_t N_2used=0;
     size_t N_unused=0; 
    #pragma omp parallel for private(i)   reduction(+:N_unused) reduction(+:N_2used) schedule (dynamic)  
      for(i=0; i<N_grid_points; i++)
       {
           if(!(*grid_point)[i].marker) N_unused++;
           if((*grid_point)[i].marker==2) N_2used++;
       }
       
        *used_size = ((double)(N_grid_points - N_unused)/(double) N_grid_points )*area_size;
        *not_used_size = ((double) N_unused/(double) N_grid_points)*area_size;
        *poly_used_size = ((double) N_2used/(double) N_grid_points)*area_size;   
    free(grid_point);    
    free(area);         
  }


void format_data(char *sourse_string)
 {
     size_t comma_counter =0;
    for(size_t i=0; i<strlen(sourse_string); i++)
       {
           
           if(sourse_string[i]==',')
             {
                 if(comma_counter<2)
                  {
                    sourse_string[i]=' ';
                    comma_counter++; 
                  }
                   else
                  {
                   sourse_string[i]='\n'; 
                   comma_counter=0;
                  }
             }
            
       } 
 }

void json_fields_to_fld(char *sourse_filename, int field_id,char *dst_filename)
{
    FILE *f;long len;char *data;
    
    f=fopen(sourse_filename,"rb");fseek(f,0,SEEK_END);len=ftell(f);fseek(f,0,SEEK_SET);
    data=(char*)malloc(len+1);
    if(fread(data,1,len,f));
    fclose(f);
    char field_id_format[255];
    sprintf(field_id_format,"%d",field_id);
    cJSON *json = cJSON_Parse(data);
    cJSON *format = cJSON_GetObjectItem(json,field_id_format);
    format_data(cJSON_GetObjectItem(format,"points")->valuestring);
     FILE *dst=fopen(dst_filename,"w");
    fprintf(dst,"%s",cJSON_GetObjectItem(format,"points")->valuestring);
    cJSON_Delete(json);
    fclose(dst);
    free(data);
}


void json_messages_to_tmsg(char *sourse_filename,char *destination_filename)
{
    FILE *f;long len;char *data;
    
    f=fopen(sourse_filename,"rb");fseek(f,0,SEEK_END);len=ftell(f);fseek(f,0,SEEK_SET);
    data=(char*)malloc(len+1);
    if(fread(data,1,len,f));
    fclose(f);
    cJSON *json = cJSON_Parse(data);
    FILE *dst = fopen(destination_filename,"w");
  for(size_t i=0; i<cJSON_GetArraySize(json); i++)
    fprintf(dst,"%.10lf %.10lf\n",cJSON_GetObjectItem(cJSON_GetArrayItem(json,i),"x")->valuedouble,cJSON_GetObjectItem(cJSON_GetArrayItem(json,i),"y")->valuedouble);
    cJSON_Delete(json);
    fclose(dst);
    free(data);
}

void get_coordinates( char *line, long double *longtitude,long double *latitude) 
 {
     char *p = strtok(line," ");
    if(p) *longtitude  = atof(p);
    p = strtok(NULL," ");
    if(p) *latitude = atof(p);
    p = strtok(NULL," ");
    p=NULL;
 }
 
 
 void get_area(const char *track_path,const char *geofence_path, const char *area_path)
 {
     FILE *geofence_file = fopen(geofence_path,"r");
     size_t N_geofence=0;
     char line[255];
     while(!feof(geofence_file))
      {
          if(fgets(line,sizeof(line),geofence_file)) N_geofence++;
      }
      fclose(geofence_file);
      POINT (*geofence)[N_geofence] = malloc(sizeof(*geofence));
      geofence_file = fopen(geofence_path,"r");
      for(size_t i=0; i<N_geofence&&!feof(geofence_file);i++)
       {
           fgets(line,sizeof(line),geofence_file);
           get_coordinates(line,&(*geofence)[i].latitude,&(*geofence)[i].longtitude);
       }
      fclose(geofence_file);
      FILE *track=fopen(track_path,"r");
      FILE *area = fopen(area_path,"w");
      long double latitude;
      long double longtitude;
          while(!feof(track))
           {
               fgets(line,sizeof(line),track);
               get_coordinates(line,&longtitude,&latitude);
              if(point_in_geofence(N_geofence,geofence,longtitude,latitude))
                {
                    fprintf(area,"%.10Lf  %.10Lf\n",longtitude,latitude);
                }
           }
          fclose(track);
          fclose(area);
      free(geofence);
 }

void reducing_points(const char *src_path, const char *dst_path,const double reducing_k )
  {
      
        double derivative(POINT p1 ,POINT p2 )
         {
            return  (p2.latitude - p1.latitude)/(p2.longtitude-p1.longtitude);
         }
  
     FILE *src  = fopen(src_path,"r");
     if(!src) 
         {
            fprintf(stderr,"Error! File not opened!\n");
            exit(1); 
         }
     size_t N=0;
     char line[255];
     
     while (!feof(src))
      {
        fgets(line,sizeof(line),src);
        if(!feof(src)) N++;
      }
     fclose(src);
     
     POINT (*p)[N] = malloc(sizeof(*p));
     size_t i=0;
     src = fopen(src_path,"r");
     while (!feof(src)&&i<N)
      {
          fgets(line,sizeof(line),src);
          get_coordinates(line,&(*p)[i].latitude,&(*p)[i].longtitude);
          i++;
      }
      fclose(src);
      
        FILE *dst = fopen(dst_path,"w");
        if(!dst) 
         {
            fprintf(stderr,"Error! File not opened!\n");
            exit(1); 
         }
         
        for(i=1; i<N-1; i++)
         {
             if(!((derivative((*p)[i],(*p)[i-1])/derivative((*p)[i],(*p)[i+1]))<=(2-reducing_k)&&
             (derivative((*p)[i],(*p)[i-1])/derivative((*p)[i],(*p)[i+1]))>=reducing_k))
              fprintf(dst,"%.10Lf   %.10Lf\n",(*p)[i].latitude,(*p)[i].longtitude);
         }
         
        fclose(dst);
     free(p);
  }
  
  double field_size(char *field_path , double grid_size)
 {
     FILE *geofence = fopen(field_path,"r");
     size_t N_field = 0;
     char line[255];
     while(!feof(geofence))
      {
          if(fgets(line,sizeof(line),geofence) ) N_field++;
      }
      fclose(geofence);
      POINT (*field)[N_field] = malloc(sizeof(*field));
      geofence = fopen(field_path,"r");
      for(size_t i=0; i<N_field&&!feof(geofence); i++)
        {
            fgets(line,sizeof(line),geofence);
            get_coordinates(line,&(*field)[i].latitude,&(*field)[i].longtitude);
        }
        double size  = (double)size_area(field, N_field,grid_size);
       free(field);
       return size;
 }

void area_analyser(char *field_path, char *area_path,double area_size, double width, double  grid_size, double *used_size, double *not_used_size, double *poly_used_size)
 {
     FILE *field_file = fopen(field_path,"r");
     size_t N_field=0;
     char line[255];
     while(!feof(field_file))
      {
          if(fgets(line,sizeof(line),field_file)) N_field++;
      }
      fclose(field_file);
      field_file = fopen(field_path,"r");
     POINT (*field)[N_field] = malloc(sizeof(*field));
     for(size_t i=0; i<N_field&&!feof(field_file); i++)
      {
          fgets(line,sizeof(line),field_file);
          get_coordinates(line,&(*field)[i].latitude,&(*field)[i].longtitude);
      }
      fclose(field_file);
      FILE *area_file = fopen(area_path,"r");
      size_t N_area=0;
      while(!feof(area_file))
       {
           if(fgets(line,sizeof(line),area_file)) N_area++;
       }
      fclose(area_file);
      POINT (*area)[N_area] = malloc(sizeof(*area));
      area_file = fopen(area_path,"r");
      for(size_t i=0; i<N_area&&!feof(area_file); i++)
        {
            fgets(line,sizeof(line),area_file);
            get_coordinates(line,&(*area)[i].longtitude,&(*area)[i].latitude);
        }
      fclose(area_file);
      double used, unused, p_used;
      area_analyse(field,N_field,area,N_area,area_size,width,grid_size,&used,&unused,&p_used);
      *used_size=used;
      *not_used_size=unused;
      *poly_used_size=p_used;
     free(field);
     free(area);
 }  

 double side(POINT a, POINT b)
  {
      return sqrt( (b.longtitude-a.longtitude)*(b.longtitude-a.longtitude) + (b.latitude-a.latitude)*(b.latitude-a.latitude) );
  }

 double field_perimeter(char *fld_sourse)
 {
     FILE *field_file = fopen(fld_sourse,"r");
     size_t N_field =0;
     char line[255];
     while(!feof(field_file))
      {
          if(fgets(line,sizeof(line),field_file)) N_field++;
      }
      fclose(field_file);
      POINT (*field)[N_field] = malloc(sizeof(*field));
      field_file = fopen(fld_sourse,"r");
      for(size_t i=0; i<N_field&&!feof(field_file); i++)
       {
           fgets(line,sizeof(line),field_file);
           get_coordinates(line,&(*field)[i].latitude,&(*field)[i].longtitude);
       }
       size_t i;
       #pragma omp parallel for private(i) shared(field) schedule(dynamic)
       for(i=0;i<N_field; i++)
        {
			convert_coordinates(&(*field)[i].longtitude,&(*field)[i].latitude);
		}
       fclose(field_file);
       double perimeter = side((*field)[0],(*field)[N_field-1]);
       
       for(i=0; i<N_field-1; i++)
        {
            perimeter+=side((*field)[i],(*field)[i+1]);
        }
        free(field);
       return perimeter*0.001;
 }
 
 void field_size_and_perimeter(char *json_field_file,int field_id ,double grid_size,double *area_size, double *perimeter)
  {    
       char tmp_file_fld[5000];
       strcpy(tmp_file_fld,json_field_file);
       strcat(tmp_file_fld,".fld");
       json_fields_to_fld(json_field_file,field_id,tmp_file_fld);
       *area_size=field_size(tmp_file_fld,grid_size);
       *perimeter=field_perimeter(tmp_file_fld);
       remove(tmp_file_fld);
  }
  
 void json_field_analyse(char *json_field_file,int field_id,char *json_message_file,double area_size, double width, double grid_size,double *used_field_size,double *not_used_field_size, double *poly_used_field_size)
  {
     char tmp_field[6000];
     strcpy(tmp_field,json_field_file);
     strcat(tmp_field,".fld");
     
     char tmp_track[6000];
     strcpy(tmp_track,json_message_file);
     strcat(tmp_track,"_track.tmsg");
     
     char tmp_area[6000];
     strcpy(tmp_area,json_message_file);
     strcat(tmp_area,"_area.tmsg");
     json_fields_to_fld(json_field_file,field_id,tmp_field);
     json_messages_to_tmsg(json_message_file,tmp_track);
     get_area(tmp_track,tmp_field,tmp_area);  
     area_analyser(tmp_field,tmp_area,area_size,width,grid_size,used_field_size,not_used_field_size,poly_used_field_size);
     remove(tmp_field);
     remove(tmp_track);
     remove(tmp_area);
  }

 double determinant(POINT p1, POINT p2)
 {
     return p1.longtitude*p2.latitude - p2.longtitude*p1.latitude;
 }   
 
 double poly_size(POINT (*poly)[], size_t N_poly)
  {
    double S=0;
    for(size_t i=0; i<N_poly-1; i++)
     {
         S+=determinant((*poly)[i],(*poly)[i+1]);
     }
     return abs(S)*0.5;
  }
 
 double json_polygon_size(char *json_field_sourse,int field_id)
  {
      double size;
      char tmp_path[5000];
        strcpy(tmp_path,json_field_sourse);
        strcat(tmp_path,".fld");
        json_fields_to_fld(json_field_sourse,field_id,tmp_path);
       
        size_t N_field=0;
        FILE *src = fopen(tmp_path, "r");
        char line[255];
         while(!feof(src))
          {
              if(fgets(line,sizeof(line),src)) N_field++;
          }
        fclose(src);
        POINT (*field)[N_field] = malloc(sizeof(*field));
        src = fopen(tmp_path, "r");
            for(size_t i=0; i<N_field; i++)
            {
                fgets(line,sizeof(line),src);
                get_coordinates(line,&(*field)[i].latitude,&(*field)[i].longtitude);
            }
            
            size_t i;
            #pragma omp parallel for private(i) shared(field) schedule(dynamic)
            for(i=0; i<N_field; i++)
             {
				 convert_coordinates(&(*field)[i].longtitude,&(*field)[i].latitude);
			 }
			 
        fclose(src);
        remove(tmp_path);
        size= poly_size(field,N_field);
        size/=10000.;
        free(field);
        return size;
  }  
  


 double json_field_size(char *json_sourse, int field_id, double grid_size)
  {
      char tmp[6000];
      strcpy(tmp,json_sourse);
      strcat(tmp,".fld");
      json_fields_to_fld(json_sourse,field_id,tmp);
      double size = field_size(tmp,grid_size);
      remove(tmp);
      return size;
  }
 
 double json_field_perimeter(char *json_field_sourse, int  field_id)
  {
      char tmp[6000];
      strcpy(tmp,json_field_sourse);
      strcat(tmp,".fld");
      json_fields_to_fld(json_field_sourse,field_id,tmp);
      double perimeter = field_perimeter(tmp);
      remove(tmp);
      return perimeter;
  }
  
  double circle_geofence_square(double radius)
  {
	  return M_PI*radius*radius/10000.;
  }
 
 double circle_geofence_perimeter(double radius)
  {
	  return 2.*M_PI*radius/1000.;
  }

double radial_distance(POINT c, POINT p)
 {
	 return sqrt( (c.longtitude - p.longtitude)*(c.longtitude - p.longtitude) + (c.latitude - p.latitude)*(c.latitude - p.latitude) );
 }

int point_in_circle_geofence(POINT center, POINT current, double radius)
 {
	 return (radial_distance(center,current)<=radius)?1:0;
 }

int circle_geofence_analyser_tmsg(double center_longtitude,double center_latitude, double radius,
 double width , double grid_size,char *track_path,double *processed, double *not_processed, double *poly_processed)
 {
	size_t N_source_track =0;
	char line[256];
	FILE *track_file  = fopen(track_path,"r");
	 if(!track_file) return 0;
	 while( !feof(track_file) )
	  {
		if(fgets(line,sizeof(line),track_file) ) N_source_track++;  
	  }
	 fclose(track_file);
	 
	 POINT (*source_track) [N_source_track] = malloc(sizeof(*source_track) );
	 
	 if(!source_track) return 0;
	 if(!track_file) return 0;
	 
	 track_file  = fopen(track_path,"r");
	  for(size_t i=0; i<N_source_track; i++ )
	      {
			if(fgets(line,sizeof(line),track_file) )
			 {
				 get_coordinates(line,&(*source_track)[i].longtitude, &(*source_track)[i].latitude);
				 (*source_track)[i].marker =0;
		     }  
		  }
	  fclose(track_file);  
	  
	  size_t i;
	  #pragma omp parallel for private(i) shared(source_track) schedule(dynamic)
      for( i=0; i<N_source_track; i++)
        {
			convert_coordinates(&(*source_track)[i].longtitude,&(*source_track)[i].latitude);
		}
		
	 POINT center ={center_longtitude,center_latitude,0};
	  convert_coordinates(&center.longtitude,&center.latitude);
	  
	 size_t N_dst_track=0;
	  for(i=0; i<N_source_track; i++)
	    {
		 if(point_in_circle_geofence(center,(*source_track)[i],radius) ) N_dst_track++; 
		}
		 if(!N_dst_track)  return 0;
	  POINT (*dst_track)[N_dst_track] = malloc( sizeof(*dst_track) );
	  if(!dst_track) return 0;
	  size_t j=0;
	   for(i=0; i<N_source_track; i++)
	    {
			if(point_in_circle_geofence(center,(*source_track)[i],radius) )
			  {
				  (*dst_track)[j] = (*source_track)[i];
				  j++;
			  }
		}
	 free(source_track);
	 
	 size_t Npx = (size_t) (2*radius/grid_size);
     size_t Npy = (size_t) (2*radius/grid_size);
     
     POINT grid_p;
     size_t N_grid=0;
        for(size_t i=0; i<Npx; i++)
		    for(size_t j=0; j<Npy; j++)
		         {
					 grid_p.longtitude = (center.longtitude - radius) +i*grid_size;
					 grid_p.latitude  = (center.latitude    - radius) + j*grid_size;
					 if(point_in_circle_geofence(center,grid_p,radius) )  N_grid++;
				 }
	 
	 POINT (*grid_point)[N_grid] = malloc(sizeof(*grid_point) );
	 if(!grid_point) return  0;
	 size_t k=0;
	 for(size_t i=0; i<Npx; i++)
		    for(size_t j=0; j<Npy; j++)
		         {
					 grid_p.longtitude = (center.longtitude - radius) +i*grid_size;
					 grid_p.latitude  = (center.latitude    - radius) + j*grid_size;
					 if(point_in_circle_geofence(center,grid_p,radius) ) 
					   {
						   (*grid_point)[k] = grid_p;
						   (*grid_point)[k].marker=0;
						   k++;
					   }
				 }
				 
	 width*=0.5;
	 
	 #pragma omp parallel for private(i,j)  shared(dst_track,grid_p) schedule(dynamic)
	 for(i=0; i<N_grid; i++)
	   for(j=0; j<N_dst_track-1; j++)
	         {
				if(dist_to_point((*dst_track)[j],(*dst_track)[j+1],(*grid_point)[i])<=width
				   &&(*grid_point)[i].marker<2)
				       (*grid_point)[i].marker++;  
			 }
	 free(dst_track);
	 size_t N_used=0, N_not_used=0, N_poly_used=0;
	     for(size_t i=0; i<N_grid; i++)
			    {
				 if(!(*grid_point)[i].marker ) N_not_used++;
				 if((*grid_point) [i].marker ) N_used++;
				 if((*grid_point)[i].marker>1) N_poly_used++;
				}
	 *processed = circle_geofence_square(radius)*((double) N_used / (double) N_grid);
	 *not_processed = circle_geofence_square(radius)*( (double) N_not_used / (double) N_grid);
	 *poly_processed = circle_geofence_square(radius)*( (double) N_poly_used / (double) N_grid);
	 free(grid_point);
	 return 1;
 }

 int json_circle_geofence_analyser(double center_longtitude, 
 double center_latitude,double radius, double width , double grid_size, char *json_path, double *processed, double *not_processed , double *poly_processed)
    {
		if(!json_path ) return 0;
		char (*tmp_path)[strlen(json_path)+strlen(".tmsg")+2] = malloc(sizeof(*tmp_path));
		strcpy((*tmp_path),json_path);
		 strcat((*tmp_path),".tmsg");
		 json_messages_to_tmsg(json_path,(*tmp_path));
		 int status =  circle_geofence_analyser_tmsg(center_longtitude,center_latitude, radius,width , grid_size,(*tmp_path),processed,not_processed,poly_processed);
		  remove((*tmp_path) );
		 if(!status) return 0;
		
		 free(tmp_path);
		 return status;
	}
	
	
 int json_is_simple_geofence(char *json_source_path , int id)
   {
	   char (*tmp_path) [strlen(json_source_path)+10] = malloc(sizeof(*tmp_path) );
	   strcpy((*tmp_path),json_source_path);
	   strcat((*tmp_path),".fld");
	   json_fields_to_fld(json_source_path, id, (*tmp_path));
	   int status = is_geofence_simple((*tmp_path));
	   remove((*tmp_path));
	   free((*tmp_path));
	   return status;
   }
