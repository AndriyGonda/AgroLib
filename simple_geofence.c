#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include "simple_geofence.h" 

 int double_equal(double a, double b){return (a-b)<=1e-308&&(b-a)<=1e-308;}
 int double_more(double a, double b) {return (a-b)>1e-308;}
 int double_less(double a, double b) { return (b-a)>1e-308;} 
 int double_more_equal(double a, double b){return (b-a)<=1e-308;}
 double double_max(double a, double b) { return (double_more(a,b))? a:b;}
 double double_min(double a, double b) { return (double_less(a,b))?a:b;}
 int equal_point(GPOINT A, GPOINT B)  {return double_equal(A.x,B.x)&&double_equal(A.y,B.y);}
 double distance(GPOINT A, GPOINT B) { return sqrt( (A.x-B.x)*(A.x-B.x) + (A.y-B.y)*(A.y-B.y) );}
 
 void get_coords_geofence( char *s,  double *longtitude, double *latitude) 
     {
	  char *p = strtok(s," ");
	  if(p) *longtitude  = atof(p);
	       p = strtok(NULL," ");
	  if(p) *latitude = atof(p);
	  p = strtok(NULL," ");
	  p=NULL;
   }
   
 void point_2_to_line(const GPOINT p1, const GPOINT p2 , GLINE *L)
	{
	   L->A = p2.y -p1.y;
	   L->B = p1.x -p2.x;
	   L->C = p1.y*(p2.x -p1.x)- p1.x*(p2.y -p1.y);
	}
 
 int cross_line(GLINE L1, GLINE L2)
	{
		double status = L1.A*L2.B - L2.A*L1.B;
		return !double_equal(status,0);
	}

 void line2_to_point(GLINE L1, GLINE L2, GPOINT *P)
  {
	 double st = L1.A*L2.B - L2.A*L1.B;
	  P->x = -(L1.C*L2.B - L2.C*L1.B) / st;
	  P->y = -(L1.A*L2.C - L1.C*L2.A) / st;
  }
  
  int at_segment(GPOINT A, GPOINT B , GPOINT P)
   {
	  if(equal_point(A,B) ) return equal_point(A,P);
	  return (double_equal( (B.x-A.x)*(P.y-A.y) , (B.y - A.y)*(P.x-A.x) )  && 
	  ( (double_more_equal(P.x,A.x) && double_more_equal(B.x,P.x) ) || ( double_more_equal(P.x,B.x) && double_more_equal(A.x,P.x) ) ) );
   }
   
 int segment_line_cross(GPOINT fL, GPOINT fR, GPOINT sL, GPOINT sR)
  {
	  double minf,maxf, mins,maxs;
	  GPOINT ZeroPoint={0.,0.};
	  minf = double_min(distance(fL,ZeroPoint),distance(fR,ZeroPoint));
	  maxf = double_max(distance(fL,ZeroPoint),distance(fR,ZeroPoint));
	  mins = double_min(distance(sL,ZeroPoint),distance(sR,ZeroPoint));
	  maxs = double_max(distance(sL,ZeroPoint),distance(sR,ZeroPoint));
	  if(double_equal(minf,maxs)||double_equal(maxf,mins)) return 0;
	  if(double_more(mins,maxf)||double_more(minf,maxs) ) return 1;
	  else return 2;
  }
  
  int segment_cross(GPOINT fL, GPOINT fR, GPOINT sL, GPOINT sR)
     {
		 GPOINT rs;
		 GLINE L1,L2;
	  point_2_to_line(fL,fR,&L1);
	  point_2_to_line(sL,sR,&L2);
	  if( cross_line(L1,L2))
	    {
			line2_to_point(L1,L2,&rs);
		   if(equal_point(rs,fL)||equal_point(rs,fR)||equal_point(rs,sL)||equal_point(rs,sR)) return 5;
			else
			 if(at_segment(fL,fR,rs) && at_segment(sL,sR,rs)) return 7;
			  else 
			    if(at_segment(fL,fR,rs) && at_segment(sL,sR,rs)) return 6;
			      else return 4; 
		}
		else
		  if(double_equal(L1.A*L2.B,L2.A*L1.B) &&!double_equal(L1.C,L2.C) ) return 3;
		   else
		   return segment_line_cross(fL,fR,sL,sR);
	 }
  
  int is_polygon_simple(GPOINT (*point)[], size_t N_point)
   {
	   int status=0;
	   for(size_t i=0; i<N_point-1; i++)
	     for(size_t j=i+1;j<N_point;j++)
	         {
				 switch(segment_cross((*point)[i],(*point)[i+1],(*point)[j],(*point)[j+1])%N_point)
				  {
					  case 0: status=0; break;
					  case 2: status=0; break;
					  case 6: status=0; break;
					  case 7: status=0; break;
					  default: status =1;
				  }
			 }
	return status;     
   }

	int is_geofence_simple(char *geofence_path)
	 {
		 size_t N_geofence =0;
		char line[256];
		
		FILE *src = fopen(geofence_path,"r");
		while(!feof(src))
		 {
			 if(fgets(line,sizeof(line),src) ) N_geofence++;
		 }
		 fclose(src);
		 
		 GPOINT (*geofence)[N_geofence] = malloc(sizeof(*geofence));
		 src= fopen(geofence_path,"r");
		 for(size_t i = 0 ; i<N_geofence&&!feof(src); i++ )
		   {
			   if(fgets(line, sizeof(line),src) )
			    get_coords_geofence(line,&(*geofence)[i].x,&(*geofence)[i].y);
		   }
		   fclose(src);
		   int status =is_polygon_simple(geofence,N_geofence);
		   free(geofence);
		   return status;
		   
	 }
	 
  /*Created by Andriy Gonda*/
