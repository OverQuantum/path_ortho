/*
path_ortho

Orthogonalization of paths (polygones) with simplification

Copyright © 2013 OverQuantum

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.


Author contacts:
http://overquantum.livejournal.com
https://github.com/OverQuantum

Project homepage:
https://github.com/OverQuantum/path_ortho




History
2013.03.13 idea
2013.03.15 first probe on VB6
2013.03.16 C++ version started
2013.03.17 working
2013.03.17 SVG path-d syntax removed
2013.03.17 works fine, RC1

TODO:
- handle unclosed paths
- consider multiple paths
- external base vector




*/

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS  //This is to stop MSVC complaining about sprintf() and so on
#endif

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

typedef double _float;
typedef __int32 i32;

//class for storing path
class Path
{
public:
	_float* x; //node coordinates
	_float* y;
	i32 num; //number of active node
	i32 alloc; //allocated size

	Path();//constructor
	~Path();//destructor

	void Alloc(i32 size);//(re)allocate arrays with cleaning
};

//constructor
Path::Path()
{
	x = NULL;
	y = NULL;
	num = 0;
	alloc = 0;
}

//destructor
Path::~Path()
{
	if (x) delete x;
	if (y) delete y;
}

//(re)allocate arrays with cleaning
void Path::Alloc(i32 size)
{
	if (x) delete x;
	if (y) delete y;
	x = new _float[size];
	y = new _float[size];
	alloc = size;
	num=0;
}

//Orthogonalize path
//dest - for output
//src - input
//collapseLen - minimal vector length in result, smaller will be collapsed
//return: 0 - OK, other - error
i32 OrthogonalizePath(Path *dest, Path *src, _float collapseLen)
{
	i32 i,j,num, num1, dirprev, istart, num2;
	_float x1,y1,x2,y2,x4,y4,c;
	_float xbase,ybase,xbase2,ybase2;
	_float d,v2,v4,sumv2,sumv4;
	i32 *dir;//flags of direction for all path vectors, 0 - along, 1 - perpendicular
	i32 *dirgroup;//group of directions for result path

	num = src->num;

	//check number of nodes and closing state
	if (num<5)
	{
		printf("ERROR: closed path with less than 5 nodes will just collapse\n");
		return 1;
	}
	if (src->x[0]!=src->x[num-1] || src->y[0]!=src->y[num-1])
	{
		printf("ERROR: path must be closed\n");
		return 2;
	}

	//calculate two versions of base direction vector
	//1st on vectors in range (0;90) degrees
	//2nd on vectors in range (-45;45) degrees
	xbase=0;ybase=0;xbase2=0;ybase2=0;
    for (i=1;i<num;i++)
	{
		//vector from one node to another
        x1 = src->x[i]-src->x[i-1];
        y1 = src->y[i]-src->y[i-1];
        d = x1*x1+y1*y1;//and its length
        
        x1*=d;//increase length for increasing effect of longer vectors
        y1*=d;
        
		if (x1<0)
		{
			x1=-x1; y1=-y1;
		}
        //(x1;y1) now in range (-90;90) degr
        
        if (y1<0)
		{
            x2=-y1; y2=x1;
		}
		else
		{
            x2=x1; y2=y1;
		}
        //(x2,y2) now in range (0;90) degr

        xbase+=x2;//data goes into 1st base vector
        ybase+=y2;
        
        if (x1>abs(y1))
		{
            x4=x1;y4=y2; //in range (-45;45) degr
		}
		else
		{
            if (y1<0)
			{
                x4=-y1; y4=x1; //was in range (-90;-45) degr
			}
			else
			{
                x4=y1; y4=-x1; //was in range (45;90) degr
			}
		}
        //(x4,y4) now in range (-45;45) degr

        xbase2+=x4;//data goes into 2nd base vector
        ybase2+=y4;

	}

	//normalize both base vectors
    d=1.0/sqrt(xbase*xbase+ybase*ybase);
    xbase*=d;
    ybase*=d;
    
    d=1.0/sqrt(xbase2*xbase2+ybase2*ybase2);
    xbase2*=d;
    ybase2*=d;

	//calculate for both base vector square error
    sumv2=0;
    sumv4=0;
    for (i=1;i<num;i++)
	{
		//path vector from one node to another
        x1 = src->x[i]-src->x[i-1];
        y1 = src->y[i]-src->y[i-1];

		//for 1st base vector
        v2 = x1 * xbase + y1 * ybase;//length of path vector along base vector
        v4 = x1 * ybase - y1 * xbase;//length of path vector perpendicular to base vector
		v2*=v2;//square
		v4*=v4;
        if (v2>v4)
            sumv2+=v4; //path vector is along base vector, square error is defined by perpendicular
		else
            sumv2+=v2; //path vector is perpendicular to base vector, square error is defined by along length
    
		//for 1st base vector
        v2 = x1 * xbase2 + y1 * ybase2;//length of path vector along base vector
        v4 = x1 * ybase2 - y1 * xbase2;//length of path vector perpendicular to base vector
		v2*=v2;//square
		v4*=v4;
        if (v2>v4)
            sumv4+=v4; //path vector is along base vector, square error is defined by perpendicular
		else
            sumv4+=v2; //path vector is perpendicular to base vector, square error is defined by along length
	}

	if (sumv2>sumv4)
	{
		//square error of 1st base vector is larger, so we will use 2nd vector as base
		xbase=xbase2;
		ybase=ybase2;
	}

	//for now on xbase,ybase - base vector

	//allocate arrays
	dir = new i32[num+5];
	dirgroup = new i32[num+5];
	dest->Alloc(num);//will not be larger than source path

	//fill dir array with directions
	for (i=1;i<num;i++)
	{
		//path vector from one node to another
        x1 = src->x[i]-src->x[i-1];
        y1 = src->y[i]-src->y[i-1];

        v2 = abs(x1*xbase+y1*ybase);//length of path vector along base vector
        v4 = abs(x1*ybase-y1*xbase);//length of path vector perpendicular to base vector
        if (v2>v4)
			dir[i-1]=0;//path vector is along base vector
		else
            dir[i-1]=1;//path vector is perpendicular to base vector
	}

RecalcOrtho:

	//calculate C-parameter for group of nodes
	//group is a sequental vectors having equal dir
	//C-param is from line equation A*x+B*y = C, where (A,B) - vector perpendicular to line
    
	num1 = num-1;//number of vectors in closed source path
	istart=-1;//group not started
	num2=0;//number of groups
    dirprev=dir[num1-1];//first vector will be in group with last vectors in case of same dir

	for(i=0;i<num*2;i++)//we must cycle more than one time to group last vectors with first ones
	{
        if (dir[i%num1]!=dirprev) // %num1 is to handle more than one cycle
		{
			//dir has changed - group finished

            if (istart>=0)
			{
				c=0;
                if (dirprev==0)
				{
                    //path vector along base vector, using (B,-A) - perpendicular to base vector
                    for(j=istart;j<=i;j++)
                        c+=src->x[j%num1]*ybase-src->y[j%num1]*xbase;
				}
                else
				{
                    //path vector perpendicular to base vector, using base vector (A,B)
                    for(j=istart;j<=i;j++)
                        c+=src->x[j%num1]*xbase+src->y[j%num1]*ybase;
				}
                c=c/(1.0+i-istart); //get average C-param of line

                if (dirprev==0)
				{
                    dest->x[num2]=c; //x of dest path is used to store C-param of along-direction
					dest->x[num2+1]=c;//next node also have same C-param
				}
				else
				{
                    dest->y[num2]=c; //y of dest path is used to store C-param of ortho-direction
					dest->y[num2+1]=c;
				}
				dirgroup[num2]=istart;//keep group starting node
                num2++;//group finished
                if (i>=num1)
				{
					//we have cycled over 0 and so may finish

					//store C-param of last group into first node
                    if (dirprev==0)
                        dest->x[0]=c;
					else
                        dest->y[0]=c;

                    break;//exit from for() loop
				}
			}
            dirprev = dir[i];//next group
            istart = i;//group start node
		}
	}

	//calculate nodes of new path by equations of two lines
	//A*x+B*y=C1 - ortho-direction
	//B*x-A*y=C2 - along-direction
	//=>
	//x=(A*C1+B*C2)/(A*A+B*B)
	//y=(-A*C2+B*C1)/(A*A+B*B)
	//where (A*A+B*B) = 1
    for (i=0;i<num2;i++)
	{
        x1 = xbase*dest->y[i]+ybase*dest->x[i];
		y1 = -xbase*dest->x[i]+ybase*dest->y[i];
		dest->x[i]=x1;
		dest->y[i]=y1;
	}
	
	//last node repeat first for closed path
	dest->x[num2]=dest->x[0];
	dest->y[num2]=dest->y[0];
	dest->num = num2+1;

	//check length of vector for collapsing 
	if (collapseLen>0)
	{
		_float clen2=collapseLen*collapseLen;//checking is on squares, to not waste time for square-root
		i32 iend;
		i32 flag=0;
		//num2=dest->num;
		for (i=1;i<dest->num;i++)
		{
			//path vector from one node to another
			x1 = dest->x[i]-dest->x[i-1];
			y1 = dest->y[i]-dest->y[i-1];
			d = x1*x1+y1*y1;//and its length

			if (d<clen2)
			{
				//smaller vector found
				//we should collapse it by joining dir-group of source path, which makes this vector, with near dir-group(s)

				istart=dirgroup[i-1];//start of dir-group, which makes this vector

				//end of dir-group, could be cycled
				if (i==(dest->num-1))
					iend=dirgroup[0];
				else
					iend=dirgroup[i];

				//get dirprev from previous group for erasing this group
				//simple reversing 1<->0 is not good, as two or more sequental vectors should be collapsed as a whole
				if (istart==0)
					dirprev=dir[src->num-2];//cycled
				else
					dirprev=dir[istart-1];
				
				for (j=istart;j<=iend;j++)
					dir[j]=dirprev;//change dir of this group
				flag=1;
			}
		}
		if (flag>0)
			goto RecalcOrtho;//anything found - goto calc of C-params
	}

	//done, free arrays
	delete dir;
	delete dirgroup;
	return 0;
}

//print path by printf with specified accuracy
void PrintPath(Path *src, i32 accuracy)
{
	if (src->num<2) //zero or 1 node - nothing to print
		return;
	
	//prepare format for printf-ing coordinates
	char format[20];
	sprintf(format,"%%.%if %%.%if",accuracy,accuracy);

	//printf first node
	printf(format,src->x[0],src->y[0]);
	for (i32 i=1;i<src->num;i++)
	{
		//printf other nodes
		printf(", ");
		printf(format,src->x[i],src->y[i]);
	}
	printf("\n");
}

//Parse path from string into object
//return: 0 - OK, other - error (none yet)
i32 ParsePath(Path *dest, char* text)
{
	i32 i,len;
	i32 nodesNumEst;
	i32 node;
	char coord[30];
	_float x,y;
	i32 found_coords, start_coord;

	len = (i32)strlen(text);
	nodesNumEst = 0;

	//Estimate number of nodes by number of separators
	for (i=0;i<len;i++)
		if (text[i]==',') nodesNumEst++;

	nodesNumEst=10+nodesNumEst;//10 just in case
	
	//allocate nodes in path
	dest->Alloc(nodesNumEst);

	//now parse string
	node=0;
	found_coords=0;
	start_coord=0;
	
	for (i=0;i<=len;i++) //with null-terminator for finishing last number
	{
		if (text[i]==',' || text[i]==' ' || text[i]==0) //any delimiter, zero is null-terminator
		{
			if (i>start_coord)
			{
				//copy text of one number into separate array to call ato* function
				strncpy(coord,text+start_coord,i-start_coord);
				coord[i-start_coord]=0;//add null-terminator

				if (found_coords==0)
					x = atof(coord);
				else if (found_coords==1)
				{
					y = atof(coord);

					//x and y found - node complete
					dest->x[dest->num] = x;
					dest->y[dest->num] = y;
					dest->num++;
				}
				//3rd and further coordinates are ignored
				found_coords++;

				if (text[i]==',')
					found_coords=0;//if ',' then go to next node
			}
			start_coord=i+1;//sequental delimiters treated as one
		}
	}

	return 0;
}

//Usage: path_ortho [params] "<data>"
//Example: path_ortho a2 c10 "6218 8805, 6295 8675, 6501 8798, 6425 8927, 6218 8805"
//
//params:
//  aN - accuracy of output, N - integer from 0 to 20, 0 by default
//  cF - collapse len, F - float >=0, 0 by default
//data - source path
//  format: x y (,x y)
//  space delimits coordinates from each other, comma delimits nodes
//  one or zero coordinates between commas are ignored (node not created)
//  third and further coordinates are also ignored
//output:
//  result path, format same as data
//  or error text, starts with ERROR: 
//  or nothing

int main(int argumentCount, char **arguments)
{
	i32 accuracy;
	_float collapseLen;

	Path *sourcePath;
	Path *destPath;

	if (argumentCount<2)
		return 0;//no arguments - exit
	
	accuracy = 0;
	collapseLen = 0;

	//check arguments for params
	for (i32 i=1;i<(argumentCount-1);i++)
	{
		//simple params, defined by first char
		switch(arguments[i][0])
		{
		case 'a'://accuracy
			accuracy = atoi(arguments[i]+1);
			if (accuracy<0)
				accuracy=0;
			else if (accuracy>20)
				accuracy=20;
			break;
		case 'c'://collapse len
			collapseLen = atof(arguments[i]+1);
			if (collapseLen<0)
				collapseLen = 0;
			break;
		default://unknowns - ignore
			break;
		}
	}

	//create objects
	sourcePath = new Path();
	destPath = new Path();

	while(true)
	{
		//process data
		if (ParsePath(sourcePath,arguments[argumentCount-1])!=0)
			break;
		if (OrthogonalizePath(destPath,sourcePath,collapseLen)!=0)
			break;

		//print to output
		PrintPath(destPath,accuracy);
		break;
	}

	//free memory and exit
	delete sourcePath;
	delete destPath;
	return 0;
}