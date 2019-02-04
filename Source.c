/*Written by Nitzan Cohen in C programming course - BGU 2017.*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <GL/glut.h>
// declaration on a global pointer that holds the address of a matrix that will iclude values of energy and time for task 4.
double ** PTR_To_Matrix;
typedef struct picdimensions { // using this struct to move information to other tasks.
	unsigned int wpix;
	unsigned int hpix;
}picdimensions;
typedef struct GAGLIST { // structure type : linked list
	float avgpower;
	int xcord, ycord;
	int size;
	struct GAGLIST *next;
}GAGLIST;
typedef struct details { // structure used inside recursion
	int highestx;
	int lowesty;
	int counter;
	long int sum;
	int lowestx;
	int highesty;
} details;
typedef struct dataforprint
{
	int x;
	int gap;
}dataforprint;
typedef struct LC {
	int x;
	int y;
	struct LC *next;
}LC;
typedef struct GAG // include all the veriably that we need for GAG in all the progrem
{
	int x;
	int y;
	int size;
	double intensity;
	struct GAG *ptr1;
	double short_lane1;
	struct GAG *ptr2;
	double short_lane2;
	struct GAG *ptr3;
	double short_lane3;
	int counter;
	char visit;

}GAG;
typedef struct Shortest_routh
{
	double sum;
	int count;
	int flag;
}Shortest_route;

int array_size(); // prototypes
GAG* Readingfromdust(int);
void Top_5(GAG*, int);
GAG*Map(GAG*, int);
int chack_input(int);
GAG *Tree(GAG*, GAG*);
GAG *Best_routh(GAG*, Shortest_route*, double, GAG*);
LC *Map_route(GAG*, Shortest_route*, double, GAG*, LC*, LC*);
void Create_BMP();
void mark_the_shorter_path_in_matrix(LC*, int**, double, double);
void print_path_by_matrix(int**,int,int);
void prints_pluses_on_the_GAGs(GAG*, int, int**,int,int);
void deallocate_tree(GAG*);
void ChancesForGAG(GAG*, int, int, int);
void minMaxGAGs(GAG*, int, dataforprint*);
void Create_BMP_max();
void print_area_with_max_GAGS(dataforprint*, int, int);
double Calculatonofthetraillength(GAG*, int);
double startLocationOfaGAG(int, GAG*);
void energyCAlculation(double, GAG*, int);
void myKey(unsigned char);
void WriteText(double, double, const char *);
void Warning_text(double, double, char *);
void display();
void Glut();
unsigned char **convertbmptomatrix (FILE *, unsigned char**, unsigned int, unsigned int);
GAGLIST *convertmatrixtolist(unsigned char**, int, int,unsigned int,unsigned int);
GAGLIST *findgag(unsigned char**, int, int, unsigned int, unsigned int, details*,GAGLIST*);
void printdustfile(FILE *, errno_t,GAGLIST*);


void main()
{
	dataforprint *ptrdata=NULL;//use in "print_area_with_max_GAGS" function for extract information from "minMaxGAGs" function.
	GAG *GAGarray = NULL;// use to point on the array of GAGs
	GAG *Root;//point on the root of the tree.
	Shortest_route *pointer;
	LC *link;
	LC *bighead = NULL;//link list of the shortest way.
	LC *head;// copy the link list of the shrtest way.
	LC *free_link_list;// use to free the dynemic memory;
	LC *tmp20;
	int **matrix2 = NULL;// matrix that simulate the bmp stracture(the pixels area).
	int Number_of_GAGs = 0;
	int tree = 0;// read the inputs from the users and use it to be the root of the tree.
	double m = 0, c = 0;// use for the line equation we use in "mark_the_shorter_path_in_matrix" function
	///////////////////////////////////////////////////////////////////
	//Extracting data and preserving the data for the whole functions
	struct picdimensions dims;
	dims.wpix = 0; //intiallizing
	dims.hpix = 0; 
	errno_t err;
	FILE *p1;
	err = fopen_s(&p1, "d2.bmp", "rb");
	if (err != 0) //check if the file can be open
	{
		printf("Problem open the file d2.bmp");
	}
	else // File open successfuly
	{
		fseek(p1, 18, SEEK_SET);
		for (int i = 0; i < 4; i++) // bitwise operation to recieve the info inside binary file (required to scan 4 bytes)
			dims.wpix |= fgetc(p1) << (i * 8);
		fseek(p1, 22, SEEK_SET);
		for (int i = 0; i < 4; i++)
			dims.hpix |= fgetc(p1) << (i * 8);
		fclose(p1);
	}
	////////////////////////////////////////////////////////////////////////
	int MenuChoice = 0;
	int userChoice = 0;
	Number_of_GAGs = array_size();
	if (Number_of_GAGs != -2)
	{
		GAGarray = Readingfromdust(Number_of_GAGs);
	}
	while (MenuChoice != 9)
	{

		printf_s("--------------------------\nNano mechine services\n--------------------------\nMenu:\n");
		printf_s("1. Scan dust.\n2. Top 5.\n3. Route it.\n4. Energy cosumption report.\n5. Students addition\n9. Exit.\n");
		printf_s("Enter choice: ");
		scanf_s("%d", &MenuChoice);
		printf_s("\n");
		if (MenuChoice == 1)
		{
			errno_t err;
			FILE *p1;
			int count = 0;
			unsigned int wpix = 0, hpix = 0;
			unsigned char **matrix;
			err = fopen_s(&p1, "d2.bmp", "rb");
			if (err != 0) //check if the file can be open
				printf("Problem open the file d2.bmp");
			else // File open successfuly
			{
				fseek(p1, 18, SEEK_SET);
				for (int i = 0; i < 4; i++) // bitwise operation to recieve the info inside binary file (required to scan 4 bytes)
					wpix |= fgetc(p1) << (i * 8);
				fseek(p1, 22, SEEK_SET);
				for (int i = 0; i < 4; i++)
					hpix |= fgetc(p1) << (i * 8);
				int superpadding = 0;
				superpadding = (4 - (wpix * 3) % 4) % 4; // allocating memory with padding bytes
				matrix = (unsigned char**)malloc((sizeof(unsigned char*)*hpix)); // allocating dynamic memory to 2D matrix
				for (int i = 0; i < hpix; i++)
					matrix[i] = (unsigned char*)malloc((sizeof(unsigned char)*wpix) + superpadding); // 2D matrix with dimensions: wpix*hpix
				matrix = convertbmptomatrix(p1, matrix, wpix, hpix); // converting the picture to matrix
				GAGLIST *head = convertmatrixtolist(matrix, 0, 0, wpix, hpix); // converting the matrix to linked list
				err = fclose(p1); // closing the file
				GAGLIST *tmp = head;
				if (head != NULL)
				{
					printf("Coordinate (%d,%d)\n", head->xcord, head->ycord);
					printf("Size %d\n", head->size);
					printf("Average intensity %0.1f\n", head->avgpower);
				}
				else
					printf("Total of 0 big dust particles.\n\n");
				printdustfile(p1, err, tmp); // creating the dust.txt file
				tmp = head;
				while (head != NULL) // free the memory allocated for the linked list
				{
					count++;
					tmp = head;
					head = head->next;
					free(tmp);
				}
				if (count > 0)
				{
					printf("Total of %d big dust particles.\n", count);
					printf("\n");
				}
				for (int i = 0; i < hpix; i++)
				{// free the memory allocated for the matrix
					matrix[i];
				}
				free(matrix);
				if (GAGarray == NULL)
				Number_of_GAGs = array_size();
				GAGarray = Readingfromdust(Number_of_GAGs);
			}
		}
		else if (MenuChoice == 2) 
		{
			// printing the 5 biggest dust areas
			if (GAGarray != NULL)
			{
				Top_5(GAGarray, Number_of_GAGs);
			}
			else
				printf("Option 1 wasn't selected, therefore option 2 can not be applied.\n");
		}
		else if (MenuChoice == 3)
		{
			// printing the map route and information
			if (Number_of_GAGs == 0);
			else if (GAGarray != NULL)
			{
				Map(GAGarray, Number_of_GAGs);
				tree = chack_input(Number_of_GAGs);
				Root = malloc(sizeof(GAG));
				Root = Tree(&GAGarray[tree - 1], Root);
				pointer = malloc(sizeof(Shortest_route));
				link = malloc(sizeof(LC));
				bighead = link;
				pointer->sum = 9999999999; // finding minimum value
				pointer->count = 0;
				pointer->flag = 0;
				Best_routh(Root, pointer, 0, GAGarray);
				if (pointer->sum != 9999999999)
				{
					pointer->count = 0; // using count as flag now and not counter
					Map_route(Root, pointer, 0, GAGarray, link, bighead);
					head = bighead;// to save the data find in bighead.
					free_link_list = bighead;
					printf_s("Best route found:\n");
					while (bighead != NULL) // printing the linked list 
					{
						tmp20 = bighead;
						printf("(%d,%d)", tmp20->x, tmp20->y);
						if (bighead->next != NULL)
							printf("->");
						bighead = bighead->next;
					}
					Create_BMP();
					matrix2 = malloc(sizeof(int*)*dims.hpix);
					if (matrix2 == NULL) exit(1);//chack if the memory allocation doesn't work.
					for (int i = 0; i < dims.hpix; i++)//build the matrix to simulete the bmp pixel stracture.
					{
						matrix2[i] = malloc(sizeof(int)*dims.wpix);
						if (matrix2 == NULL) exit(2);//chack if the memory allocation doesn't work.
					}
					for (int i = 0; i < dims.hpix; i++)// run over all the matrix in nested loop and puts zeros
					{
						for (int j = 0; j < dims.wpix; j++)
						{
							matrix2[i][j] = 0;
						}
					}
					while (head->next != NULL)
					{
						m = (double)((head->next->x) - (head->x));//chack if the slope in the line equation is undefined.
						if (m == 0)
							m = 9999999999;
						else
							m = ((double)((head->next->y) - (head->y)) / (double)((head->next->x) - (head->x)));// the slope i line equation.
						c = -m*(double)head->x + (double)head->y;// the constant value line equation.
						mark_the_shorter_path_in_matrix(head, matrix2, m, c);
						head = head->next;
					}
					print_path_by_matrix(matrix2, dims.wpix, dims.hpix);
					prints_pluses_on_the_GAGs(GAGarray, Number_of_GAGs, matrix2, dims.wpix, dims.hpix);
					while (free_link_list != NULL) // deallocate link list.
					{
						tmp20 = free_link_list;
						free_link_list = free_link_list->next;
						free(tmp20);
					}
				}
				else
				{
					printf_s("There is no route from (%d,%d) to (%d,%d)\n",&GAGarray[tree-1].x, &GAGarray[tree - 1].y,GAGarray->x,GAGarray->y);
					Create_BMP();
					matrix2 = malloc(sizeof(int*)*dims.hpix);
					if (matrix2 == NULL) exit(1);//chack if the memory allocation doesn't work.
					for (int i = 0; i < dims.hpix; i++)//build the matrix to simulete the bmp pixel stracture.
					{
						matrix2[i] = malloc(sizeof(int)*dims.wpix);
						if (matrix2 == NULL) exit(2);//chack if the memory allocation doesn't work.
					}
					for (int i = 0; i < dims.hpix; i++)// run over all the matrix in nested loop and puts zeros
					{
						for (int j = 0; j < dims.wpix; j++)
						{
							matrix2[i][j] = 0;
						}
					}
					prints_pluses_on_the_GAGs(GAGarray, Number_of_GAGs, matrix2, dims.wpix, dims.hpix);
				}
				for (int i = 0; i < dims.hpix; i++) // free the memory allocated for the matrix
					free(matrix2[i]);
				free(matrix2);
				deallocate_tree(Root);
			}
			else
			{
				printf("Option 1 wasn't selected, therefore option 3 can not be applied.\n");
			}
		}
		else if (MenuChoice == 4)
		{
			//printing the energy consumption and plot
			if (Number_of_GAGs == 0)
				printf_s("There is no valid content in dust.txt\n");
			else if (GAGarray != NULL)
			{
				printf_s("Enter length of route up to 9: ");
				scanf_s("%d", &userChoice);
				while (userChoice < 1 || userChoice>9)
				{
					printf_s("Enter length of route up to 9: ");
					scanf_s("%d", &userChoice);
				}
				double length = Calculatonofthetraillength(GAGarray, userChoice);
				energyCAlculation(length, GAGarray, userChoice);
				Glut();
			}
			else
				printf("Option 1 wasn't selected, therefore option 4 can not be applied.\n");
		}
		else if (MenuChoice == 5)
		{
			// Choice 5
			if (GAGarray != NULL)
			{
				ptrdata = malloc(sizeof(dataforprint)); // allocating memory for updating outside-data.
				ChancesForGAG(GAGarray, Number_of_GAGs, dims.wpix, dims.hpix);
				minMaxGAGs(GAGarray, Number_of_GAGs, ptrdata);
				Create_BMP_max();
				print_area_with_max_GAGS(ptrdata, dims.wpix, dims.hpix);
			}
			else
				printf("Option 1 wasn't selected, therefore option 5 can not be applied.\n");
		}
		else if (MenuChoice == 9)
		{
			printf_s("Good bye!\n");
			if (GAGarray != NULL)
			{
				free(GAGarray);
			}
			free(ptrdata);
			exit(1);
		}
		else
		{
			printf_s("Bad input, try again\n\n");
		}
	}
}

unsigned char ** convertbmptomatrix(FILE *p1, unsigned char **matrix, unsigned int wpix, unsigned int hpix)
{

	//Main idea: Scan the whole PixelArray(of bmp file) and extract only one of the bytes in each pixel to a matrix.
	unsigned int offset = 0;
	unsigned char R, G, B;
	int padding = 0;
	padding = (4 - (wpix * 3) % 4) % 4; //Equation to calculate number of padding bytes in the end of a row.
	fseek(p1, 10, SEEK_SET);
	for (int i = 0; i < 4; i++) // bitwise operation to recieve the info inside binary file (required to scan 4 bytes)
		offset |= fgetc(p1) << (i * 8); //offset to reach the array
	fseek(p1, offset, SEEK_SET);
	for (int i = 0; i < hpix; i++) // the height equal to i
	{
		for (int j = 0; j < wpix; j++) // the width equal to j
		{
			R = getc(p1); // moving the pointer 1 byte
			G = getc(p1); // moving the pointer 1 byte
			B = getc(p1); // moving to the next pixel
			matrix[i][j] = R; // R,G,B have the same values
		}
		fseek(p1, padding, SEEK_CUR); // Skip padding bytes if there are.
	}
	return matrix;
}
GAGLIST *convertmatrixtolist(unsigned char **matrix, int x, int y, unsigned int wpix, unsigned int hpix)
{
	
	//Main idea: Scan the whole matrix, if there is pixel equal/high of 80 intensity, activate findgag function.
	GAGLIST *head = NULL;
	GAGLIST *tail = NULL;
	GAGLIST *temp;
	// temp2: type:GAG* pointer, use: update data inside recursion to create a temporary member of the type of linked list.
	GAGLIST *temp2;
	temp2 = (GAGLIST*)malloc(sizeof(GAGLIST)); 
	// ptr: type: details pointer, use: storing data inside recursion.
	details *ptr;
	ptr = (details*)malloc(sizeof(details));
	ptr->counter = 1; // The loops doesn't count the first pixel in the GAG, so the counter is intialized to 1.
	// Values that will easily replaced to get minimum and maximums.
	ptr->lowesty = 999999; 
	ptr->highestx = 0; 
	ptr->sum = 0;
	ptr->lowestx = 9999999; 
	ptr->highesty = 0; 

	for (int i = 0; i < hpix; i++) // the search begins
	{
		for (int j = 0; j < wpix; j++)
		{
			if (matrix[i][j] >= 80)
			{
				temp2=findgag(matrix, i, j, wpix, hpix, ptr,temp2); // Found a pixel that is high then 80, call findgag function.
				//reintialize
				ptr->counter = 1;
				ptr->lowesty = 999999;
				ptr->highestx = 0;
				ptr->sum = 0;
				ptr->lowestx = 999999;
				ptr->highesty = 0;
				if (temp2 != NULL)
				{
					if (temp2->size >= 30) // add a member to the linked list only if its high then 30.
					{
						temp=(GAGLIST*)malloc(sizeof(GAGLIST)); // adding a member to the linked list by allocating memory
						temp->avgpower = temp2->avgpower; // copy the details inside the temporary member to a solid member of the list.
						temp->size = temp2->size;
						temp->xcord = temp2->xcord;
						temp->ycord = temp2->ycord;
						temp->next = NULL;
						if (tail)
							tail->next = temp;
						else
							head = temp;
						tail = temp;
						tail->next = NULL;
					}
				}
				}
			}
		}
	free(ptr); // free allocated memory
	free(temp2);// free allocated memory
	return head; // returning the head of the list
	}



GAGLIST *findgag(unsigned char **matrix, int y, int x, unsigned int wpix, unsigned int hpix,details *ptr,GAGLIST *temp) 
{
	
	//y equal to i(height) and x equal to j(width)
	// Main idea:
	//1. run with recursion to search up,down,left,right until all the calls are closed
	//2. update data to a reference pointer structure
	//3. save the info inside it
	//4. create temporary member of the linked list and return it.
	int flag = 0;
	long int t = matrix[y][x];
	if (ptr->lowestx > x) // finding the minimum x
		ptr->lowestx = x;
	if (ptr->highestx < x) // finding the maximum x
		ptr->highestx = x;
	if (ptr->highesty < y) // finding the maximum y
		ptr->highesty = y;
	if (ptr->lowesty > y) // finding the minimum y
		ptr->lowesty = y;
	if (matrix[y][x] >= 80)
	{
		ptr->sum = ptr->sum + t;
		matrix[y][x] = 0;
		if (y - 1 >= 0 && matrix[y - 1][x] >= 80) // up
		{
			ptr->counter++;
			findgag(matrix, y - 1, x, wpix, hpix, ptr,temp);
			flag = 1;
		}
		if (y + 1 < hpix && matrix[y + 1][x] >= 80) // down
		{
			ptr->counter++;
			findgag(matrix, y + 1, x, wpix, hpix, ptr,temp);
			flag = 1;
		}
		if (x - 1 >= 0 && matrix[y][x - 1] >= 80) // left
		{
			ptr->counter++;
			findgag(matrix, y, x - 1, wpix, hpix, ptr, temp);
			flag = 1;
		}
		if (x + 1 < wpix && matrix[y][x + 1] >= 80) //right 
		{
			ptr->counter++;
			findgag(matrix, y, x + 1, wpix, hpix, ptr,temp);
			flag = 1;
		}
		if (flag == 0)
			// creating a temporary member of the linked list
			// to improve the speed of the recursion , update data only when going reverse step
		{
			
			temp->avgpower = ((float)ptr->sum / (float)ptr->counter); // casting operation
			temp->size = ptr->counter;
			temp->xcord =(ptr->lowestx) + (int)ceil(((ptr->highestx - ptr->lowestx) / 2)+0.5); // rounding operation
			temp->ycord =(ptr->lowesty) + (int)ceil(((ptr->highesty - ptr->lowesty) / 2)+0.5);
			return temp;
		}
		else if(flag==1)
		{
			return temp;
		}
	}
}

void printdustfile(GAGLIST *p1,errno_t err,GAGLIST *tmp)
{
	
	//creating the dust.txt file with the pointer and a temporary pointer that point to the head of the linked list.
	GAGLIST *temprint;
	err = fopen_s(&p1, "dust.txt", "w+");
	if (err != 0) // check if the file can be written
		printf("Problem writing the file dust.txt\n");
	fprintf_s(p1, "coordinate	size	average Intensity\n");
	fprintf_s(p1, "==========	====	=================\n");
	temprint = tmp;
	while (tmp != NULL)
	{
		fprintf_s(p1, "(%d,%d)   	", tmp->xcord, tmp->ycord);
		fprintf_s(p1, "%d\t", tmp->size);
		fprintf_s(p1, "%0.1f\n", tmp->avgpower);
		temprint = tmp;
		tmp = tmp->next;
	}
	err = fclose(p1);
}
double Calculatonofthetraillength(GAG* ptrTogag, int userchoice)
{
	/*This function will return a double type number that contains the value of the distance of a route that goes from the start of
	the first GAG in dust document to the end of the GAG that was chosen by the user.
	The function receieves poiter to an array that was built from the information in the dust form and an integer
	that reflects the number of the GAG that will be in the end of the route.
	*/
	double distance = 0;
	GAG *ptrNext = ptrTogag;

	ptrNext++;
	GAG *ptr = ptrTogag;
	for (int i = 0; i < userchoice; i++)
	{
		if (userchoice == 1)
		{
			/* In case that the route contains only one GAG the route will be defines from it's start to it's end- meaning the length will
			be equal to it's diameter sqrt(size)*/
			distance = (sqrt(ptrTogag->size));
		}
		else
		{
			if (i == 0)
			{
				//In the first GAG the route will also contain the distance from the start to it's middle, meaning radius.
				distance += (sqrt(ptrTogag->size)) / 2;
			}
			else if (i == userchoice - 1)
			{
				//In the last GAG the route will also contain the distance from the middle to the end of the GAG, meaning radius.
				distance += ((sqrt(ptrTogag->size)) / 2);

			}
			if (i <= userchoice - 2)
			{
				distance += sqrt((ptrTogag->x - ptrNext->x)*(ptrTogag->x - ptrNext->x) + (ptrTogag->y - ptrNext->y)*(ptrTogag->y - ptrNext->y));
			}
			ptrTogag++;
			ptrNext++;
		}



		//printf_s("%f \n", distance);	
	}
	return distance;
}


double startLocationOfaGAG(int GAGindex, GAG* ptrToArray)
{
	/*This function calculates the time that a GAG starts and return it. the function will recieve a index of a GAG (based on the order in the dust list)
	and the list of the GAGs, using a pointer to an array. the function uses recursive calls. the function returns 0 for the the start of the first GAG.
	for the rest of the GAGs a calcukation is made by using the start of the GAG before and the distance between the GAGs.
	*/
	double distanceFromLast = 0;
	double timeOfstart = 0;
	//GAG *ptr = ptrToArray;
	if (GAGindex == 0)
	{
		//printf_s("start");
		return 0;
	}
	else
	{
		distanceFromLast = sqrt(pow((ptrToArray[GAGindex].x - ptrToArray[GAGindex - 1].x), 2) + pow((ptrToArray[GAGindex].y - ptrToArray[GAGindex - 1].y), 2));
		distanceFromLast = distanceFromLast - sqrt(ptrToArray[GAGindex].size) / 2;
		timeOfstart = startLocationOfaGAG((GAGindex - 1), ptrToArray) + sqrt(ptrToArray[GAGindex - 1].size) / 2 * 0.01 + distanceFromLast*0.01;

	}
	//printf_s("%f", timeOfstart);
	return timeOfstart;

}
void energyCAlculation(double routelength, GAG* ptrTogag, int userschoice)
{
	/*The function caculates the Energy that the nanometric component will use for cleaning GAGs. the function recieves the length of the route that was calculated before,
	the array of the GAGs from dust and the number of GAGs that the user chose. The energy will be calculated by equation. The results and the time of the calculated results
	will be entered to a matrix which is a global item and can be used for the Glut.
	*/

	//EnergyAndTime = (double***)malloc(sizeof(double**)*3);
	//for (int a = 0; a<3; a++)
	//EnergyAndTime[a] = (double**)malloc(sizeof(double*)*20);
	//double ** EnergyAndTime = malloc(sizeof(double*) *20);
	//for (int i = 0; i < 20; i++)
	//{
	//	EnergyAndTime[i] = malloc(sizeof(double) * 3);
	//}
	double EnergyAndTime[20][3];
	int indexOf_gag = 0;
	//	double Energy = 0;
	int i = 0, j = 0;
	double p, S;
	double time = 0;
	double Energy = 0;
	double EnergyNext = 0;
	int counter = 0;
	int flag = 0;
	//EnergyAndTime[0][0] = 0;
	//EnergyAndTime[1][0] = 0;
	//using the start of a GAG function to define if the component is inside a gag ot out of it.
	double startofGAG = startLocationOfaGAG(indexOf_gag, ptrTogag);
	double startof_nextGAG = startLocationOfaGAG((indexOf_gag + 1), ptrTogag);
	for (i = 0; i < 20; i++)
		for (j = 0; j < 3; j++)
			EnergyAndTime[i][j] = 0;
	i = 0;
	while (time <= routelength*0.01)
	{
		if (userschoice == 1)
		{
			S = sqrt(ptrTogag[0].size);
			Energy = EnergyNext;
			EnergyNext = (Energy * 1 * S*0.75 + 10)*0.001 + Energy;
			if ((EnergyNext - Energy) / 0.001 >= 550 && flag != 1)
			{
				flag = 1;
				EnergyAndTime[1][2] = EnergyNext;
				EnergyAndTime[2][2] = time;
			}
			//printf_s("%f", EnergyNext - Energy));
			time += 0.001;
		}
		//caculation of S and p
		else if (time >= startofGAG && time <= (startofGAG + sqrt(ptrTogag[indexOf_gag].size)*0.01))
		{
			//staps inside of GAG
			S = sqrt(ptrTogag[indexOf_gag].size);
			p = 1;
			time += 0.001;
			Energy = EnergyNext;
			EnergyNext = (Energy*p*S*0.75 + 10)*0.001 + Energy;
			if ((EnergyNext - Energy) / 0.001 >= 550 && flag != 1)
			{
				flag = 1;
				EnergyAndTime[1][2] = EnergyNext;
				EnergyAndTime[2][2] = time;
			}
		}
		else
		{
			EnergyAndTime[i + 1][0] = time;
			EnergyAndTime[i + 1][1] = EnergyNext;
			i++;
			while (time >= (startofGAG + sqrt(ptrTogag[indexOf_gag].size)*0.01) && time <= startof_nextGAG)
			{
				//staps between GAGs
				S = 0;
				p = 0;
				time += 0.001;
				Energy = EnergyNext;
				EnergyNext = Energy + 0.01;
				if ((EnergyNext - Energy) / 0.001 >= 550 && flag != 1)
				{
					flag = 1;
					EnergyAndTime[1][2] = EnergyNext;
					EnergyAndTime[2][2] = time;
				}
				//	counter++;
			}
			EnergyAndTime[i + 1][0] = time;
			EnergyAndTime[i + 1][1] = EnergyNext;
			i++;
			indexOf_gag++;
			startofGAG = startLocationOfaGAG(indexOf_gag, ptrTogag);
			startof_nextGAG = startLocationOfaGAG(indexOf_gag + 1, ptrTogag);
		}


	}

	EnergyAndTime[0][2] = (double)flag;
	EnergyAndTime[i + 1][0] = time;
	EnergyAndTime[i + 1][1] = EnergyNext;
	i++;

	PTR_To_Matrix = malloc(sizeof(double*) * 20); // Copy the data inside local matrix to global pointer to matrix.
	for (int i = 0; i < 20; i++)
	{
		PTR_To_Matrix[i] = malloc(sizeof(double) * 3);
	}
	for (int i = 0; i < 20; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			PTR_To_Matrix[i][j] = EnergyAndTime[i][j];
			
		}
		
	}
	free(ptrTogag); // free the array of GAG.

}
void myKey(unsigned char key)
{
	/*This function defines what happens in every keyboard action to the glut window. Based on the direction in every click on the keyboard
	*/
	switch (key)
	{
	default:

		glutDestroyWindow(1);
		exit(1);
		break;
	}

}
void WriteText(double x, double y, const char *string)
{
	// This function will be used for the glut operation. It is a general function to insert GREEN text to the display of the glutwindow.
	glColor3d(0.0, 0.5, 0.0);
	glRasterPos2d(x, y);
	while (*string)
	{
		glutBitmapCharacter(GLUT_BITMAP_9_BY_15, *string++);
	}
}
void Warning_text(double x, double y, char *string)
{
	// This function will be used for the glut operation. It is a general function to insert RED text to the display of the glutwindow.
	glColor3d(0.7, 0.0, 0.0);
	glRasterPos2d(x, y);
	while (*string)
	{
		glutBitmapCharacter(GLUT_BITMAP_9_BY_15, *string++);
	}

}
void display()
{
	//That function will determine the display of the glut window. all the events including text, lines and points will be in that function.
	//Declaration on the color of the background(white)
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT);



	// Painting X and Y Axis
	glLineWidth(3);
	glColor3f(0.0, 0.0, 0.0);
	glBegin(GL_LINES);
	glVertex2f(-0.7, -0.7);
	glVertex2f(0.8, -0.7);
	glEnd();
	glBegin(GL_LINES);
	glVertex2f(-0.7, -0.7);
	glVertex2f(-0.7, 0.8);
	glEnd();
	//Arrows in the end of the axis
	glBegin(GL_LINES);
	glVertex2f(0.8, -0.7);
	glVertex2f(0.75, -0.65);
	glEnd();
	glBegin(GL_LINES);
	glVertex2f(0.8, -0.7);
	glVertex2f(0.75, -0.75);
	glEnd();
	glBegin(GL_LINES);
	glVertex2f(-0.7, 0.8);
	glVertex2f(-0.73, 0.75);
	glEnd();
	glBegin(GL_LINES);
	glVertex2f(-0.7, 0.8);
	glVertex2f(-0.67, 0.75);
	glEnd();

	//scale marks on x axis
	glBegin(GL_LINES);
	glVertex2f(0.7, -0.7);
	glVertex2f(0.7, -0.65);
	glEnd();
	glBegin(GL_LINES);
	glVertex2f(0.35, -0.7);
	glVertex2f(0.35, -0.65);
	glEnd();
	glBegin(GL_LINES);
	glVertex2f(0.0, -0.7);
	glVertex2f(0.0, -0.65);
	glEnd();
	glBegin(GL_LINES);
	glVertex2f(-0.35, -0.7);
	glVertex2f(-0.35, -0.65);
	glEnd();

	//scale marks on y axis
	glBegin(GL_LINES);
	glVertex2f(-0.7, 0.65);
	glVertex2f(-0.65, 0.65);
	glEnd();
	glBegin(GL_LINES);
	glVertex2f(-0.7, 0.35);
	glVertex2f(-0.65, 0.35);
	glEnd();
	glBegin(GL_LINES);
	glVertex2f(-0.7, 0.0);
	glVertex2f(-0.65, 0.0);
	glEnd();
	glBegin(GL_LINES);
	glVertex2f(-0.7, -0.35);
	glVertex2f(-0.65, -0.35);
	glEnd();


	// Title of the x axis
	char timeTitle[] = "Time [s]";
	WriteText(0.7, -0.9, &timeTitle);

	// Title of the y axis
	char EnergyTitle[] = "Energy [J]";
	WriteText(-0.8, 0.9, &EnergyTitle);

	// the loop will find the highst value of the energy in the matrix. this value will later determine the hightst values on the axis scales. 	
	int i = 1;
	//for (int i = 0; i < 3; i++)
	//{
	//	for (int j = 0; j < 20; j++)
	//	{
	//		printf("%lf ",PTR_To_Matrix[i][j]);
	//	}
	//	printf("\n");
	//}
	while (PTR_To_Matrix[i][1] != 0)
	{
		i++;
	}
	//double gapOftime = (0.8 + 0.7) / EnergyAndTime[i - 1][0];
	//double gapOfEnergy = (0.8 + 0.7) / EnergyAndTime[i - 1][1];


	//Scale text of the Time axis updates dynemicly
	char value[10];
	double last_stap_time = 0;
	last_stap_time = (int)(PTR_To_Matrix[i - 1][0]) + 1;
	if ((int)last_stap_time % 2 != 0)
	{
		last_stap_time = (double)last_stap_time + 1;
	}
	//converting double to string	
	snprintf(value, 10, "%.2lf", last_stap_time);
	WriteText(0.7, -0.8, &value);
	snprintf(value, 10, "%.2lf", last_stap_time / 4);
	WriteText(-0.35, -0.8, &value);
	snprintf(value, 10, "%.2lf", last_stap_time / 4 * 2);
	WriteText(0.0, -0.8, &value);
	snprintf(value, 10, "%.2lf", last_stap_time / 4 * 3);
	WriteText(0.35, -0.8, &value);


	//Scale text of the Energy axis updates dynemicly
	int final_stap_Energy = 0;
	final_stap_Energy = (int)(PTR_To_Matrix[i - 1][1]) + 1;
	while (final_stap_Energy % 10 != 0)
	{
		final_stap_Energy++;
	}
	snprintf(value, 10, "%d", final_stap_Energy);
	WriteText(-0.8, 0.65, &value);
	snprintf(value, 10, "%d", final_stap_Energy / 4);
	WriteText(-0.8, -0.35, &value);
	snprintf(value, 10, "%d", final_stap_Energy / 4 * 2);
	WriteText(-0.8, 0.0, &value);
	snprintf(value, 10, "%d", final_stap_Energy / 4 * 3);
	WriteText(-0.8, 0.35, &value);

	// inserting dots to create the graphing, the dots will be determine by the values in the matrix and converting to the "GLUT Scale".
	for (int j = 0; j < i - 1; j++)
	{
		if (j == 0)
		{
			glLineWidth(3);
			glColor3f(0.0, 0.0, 1.0);
			glBegin(GL_LINES);
			glVertex2f(-0.7, -0.7);
			glVertex2f(1.5*(PTR_To_Matrix[j + 1][0]) / last_stap_time - 0.7, 1.5*(PTR_To_Matrix[j + 1][1]) / final_stap_Energy - 0.7);
			glEnd();
		}
		else
		{
			if (i - 1 != 1)
			{
				glLineWidth(3);
				glColor3f(0.0, 0.0, 1.0);
				glBegin(GL_LINES);
				glVertex2f((PTR_To_Matrix[j][0])*(1.5) / last_stap_time - 0.7, (PTR_To_Matrix[j][1]) *(1.5) / final_stap_Energy - 0.7);
				glVertex2f((PTR_To_Matrix[j + 1][0])*(1.5) / last_stap_time - 0.7, (PTR_To_Matrix[j + 1][1])*1.5 / final_stap_Energy - 0.7);
				glEnd();
			}
		}

	}

	//title in the head of the window
	char EnergyTitleOfPage[50] = "Total energy consumption: ";
	int lengthOfString = strlen(EnergyTitleOfPage);
	char EnergyMax[10] = { '\0' };
	char continueOfsentence[4] = " [J]";
	snprintf(EnergyMax, 10, "%.3lf", PTR_To_Matrix[i - 1][1]);
	int n;
	for (int n = 0; n <= strlen(EnergyMax); n++)
	{
		EnergyTitleOfPage[lengthOfString + n] = EnergyMax[n];
	}
	lengthOfString = strlen(EnergyTitleOfPage);
	for (int m = 0; m < 4; m++)
	{
		EnergyTitleOfPage[lengthOfString + m] = continueOfsentence[m];
	}

	WriteText(-0.4, 0.9, &EnergyTitleOfPage);

	//Excessive energy allert
	// the warning will be displayed when E ̇≥550[J/s] , and will be displayed in the place the exaption occured

	if (PTR_To_Matrix[0][2] != 0)
	{
		double X_coordinate = PTR_To_Matrix[2][2] * 1.5 / last_stap_time - 0.7;
		double Y_coordinate = PTR_To_Matrix[1][2] * 1.5 / final_stap_Energy - 0.7;
		Warning_text(X_coordinate + 0.04, Y_coordinate - 0.02, "Excessive energy");
		Warning_text(X_coordinate + 0.09, Y_coordinate - 0.07, "consumption!");
		glLineWidth(3);
		glColor3f(0.7, 0.0, 0.0);
		glBegin(GL_LINE_LOOP);
		glVertex2f(X_coordinate + 0.01, Y_coordinate + 0.04);
		glVertex2f(X_coordinate + 0.01, Y_coordinate - 0.1);
		glVertex2f(X_coordinate + 0.5, Y_coordinate - 0.1);
		glVertex2f(X_coordinate + 0.5, Y_coordinate + 0.04);
		glEnd();

	}


	/* flush GL buffers */
	glFlush();

	for (int i = 0; i < 20; i++) // free the global matrix.
		free(PTR_To_Matrix[i]);
	free(PTR_To_Matrix);



}
void Glut()
{
	//This function sends events that will determine the glut's window properties.
	int foo = 1;
	char * bar[1] = { " " };
	//	double x = 8.7;
	glutInit(&foo, bar);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	//window size
	glutInitWindowSize(640, 480);
	//window coordinates to open
	glutInitWindowPosition(0, 0);
	//title for the window
	glutCreateWindow("Energy consumption in time");
	//calling for the display function
	glutDisplayFunc(display);
	glutKeyboardFunc(myKey);
	glutMainLoop();
}
/*
The function return the number of GAG we have in the Dust file.
It open the Dust file and every time recognize \n its mean that it finish a line.
The the first line arent GAG so she ruturn the number of line-2 and this is the number of GAG.*/
int array_size()
{
	errno_t checking;
	FILE *dust;
	int ch = 0;//run sign sign in every line.
	int count_line = 0;
	checking=fopen_s(&dust, "dust.txt", "rt");
	if (checking == 0)
	{
		if (!dust)//chack that the file was open.
		{
			return -1;
		}
		ch = fgetc(dust);
		while (ch != EOF)// count how many line in the dust file
		{
			if (ch == '\n')
			{
				count_line += 1;
			}
			ch = fgetc(dust);//promote ch to the next sign.
		}
		fclose(dust);
		return count_line - 2;//we dont need the 2 first line for the array.
	}
	else
		return -2;
}
/*
The function get the number of GAG we have and ruturn pointer to array,
the array is array of GAG struct.
The x,y corrdinete, size of GAG and the intensity of the GAG are reading from the Dust file.
All the other fields in the structher are inicialized for our need later.
*/
GAG* Readingfromdust(int Number_of_GAGs) //put all the information from the dust file in arrey
{
	//	int usersNumber;
	int arrayindex = 0;
	char line[60];// every line reading into the array.
	GAG *arrayofGAG;// point on the begining of the array of GAG we return in the end.
	FILE *mydust;
	char *ptrline;// use to point on the char we chack in the Dust file.
	arrayofGAG = malloc(sizeof(GAG)*(Number_of_GAGs)); //allocate dynemic memory for the arrey
	if (arrayofGAG == NULL) exit(1);//chacking.
	fopen_s(&mydust, "dust.txt", "rt");
	arrayofGAG = malloc(sizeof(GAG)*(Number_of_GAGs));
		for (int i = 0; i < Number_of_GAGs + 2; i++) // read every line and by the sing put every detail in here place in the arrey
		{
			arrayofGAG[arrayindex].ptr1 = NULL;
			arrayofGAG[arrayindex].ptr2 = NULL;
			arrayofGAG[arrayindex].ptr3 = NULL;
			arrayofGAG[arrayindex].short_lane1 = INFINITY;
			arrayofGAG[arrayindex].short_lane2 = INFINITY;
			arrayofGAG[arrayindex].short_lane3 = INFINITY;
			arrayofGAG[arrayindex].counter = 0;
			arrayofGAG[arrayindex].visit = NULL;

			fgets(line, 60, mydust); //read the line
			if (i >= 2)//dont need the first 2 lines.
			{
				for (int j = 0; j < strlen(line); j++) // run on every char in the line
				{
					if (line[j] == '(')// the place were the x found.
					{
						ptrline = &line[j + 1];// the place were the x found.
						arrayofGAG[arrayindex].x = atoi(ptrline);
					}
					else if (line[j] == ',')
					{
						ptrline = &line[j + 1];// the place were the y found.
						arrayofGAG[arrayindex].y = atoi(ptrline);
					}
					else if (line[j] == ')')
					{
						ptrline = &line[j + 5];// the place were the size found.
						arrayofGAG[arrayindex].size = atoi(ptrline);
						for (int k = j + 5; k < strlen(line); k++)
							if (line[k] == '\t')
							{
								ptrline = &line[k];// the place were the intensity found.
								arrayofGAG[arrayindex].intensity = atof(ptrline);
								break;
							}
					}
				}
				arrayindex++;//promote the index that point on the struct we need in the aaray.
			}
		}
		fclose(mydust);
		return arrayofGAG;
}
/*
The function get pointer for the array of GAG that wes built and the number of GAG.
In the end it printing the 5 Bigest GAG and close.
It is running in loop that run as the number of the GAG we have. each loop chacking if the GAG is
on of the Bigest 5 and if is swap im with the smallest in the list.put im in is place on the list
that Bigest1 is the Bigest and Bigest5 is the smallest.
*/
void Top_5(GAG *arrayofGAG, int Number_of_GAGs)
{
	GAG calibratin; //set all the pointers to zero on the relevant veriable.
	calibratin.x = 0;
	calibratin.y = 0;
	calibratin.size = 0;
	calibratin.intensity = 0;
	GAG *Bigest_1 = &calibratin;
	GAG *Bigest_2 = &calibratin;
	GAG *Bigest_3 = &calibratin;
	GAG *Bigest_4 = &calibratin;
	GAG *Bigest_5 = &calibratin;
	GAG *temp = &calibratin;// using for switch pointers.
	for (int i = 0; i < Number_of_GAGs; i++) // loop that run all over the arrayofGAG
	{
		//the loop always put the GAG with the biggest number in Bigest1 and in the begin chack if it need to get inside the top 5 biggest GAG.
		if (arrayofGAG[i].size > Bigest_5->size)
		{
			Bigest_5 = &arrayofGAG[i];
			if (Bigest_5->size > Bigest_4->size)
			{   //switch
				temp = Bigest_5;
				Bigest_5 = Bigest_4;
				Bigest_4 = temp;
				if (Bigest_4->size > Bigest_3->size)
				{  //switch
					temp = Bigest_4;
					Bigest_4 = Bigest_3;
					Bigest_3 = temp;
					if (Bigest_3->size > Bigest_2->size)
					{  //switch
						temp = Bigest_3;
						Bigest_3 = Bigest_2;
						Bigest_2 = temp;
						if (Bigest_2->size > Bigest_1->size)
						{  //switch
							temp = Bigest_2;
							Bigest_2 = Bigest_1;
							Bigest_1 = temp;
						}
					}
				}
			}
		}
	}
	printf_s("Biggest dust particles found:\n");
	printf_s("coordinate\tsize\taverage Intensity\n");
	printf_s("==========\t====\t=================\n");
	if (Number_of_GAGs > 0)
	{
		printf_s("(%d,%d)   \t%d\t%.1lf\n", Bigest_1->x, Bigest_1->y, Bigest_1->size, Bigest_1->intensity);
		if (Number_of_GAGs > 1)
		{
			printf_s("(%d,%d)   \t%d\t%.1lf\n", Bigest_2->x, Bigest_2->y, Bigest_2->size, Bigest_2->intensity);
			if (Number_of_GAGs > 2)
			{
				printf_s("(%d,%d)   \t%d\t%.1lf\n", Bigest_3->x, Bigest_3->y, Bigest_3->size, Bigest_3->intensity);
				if (Number_of_GAGs > 3)
				{
					printf_s("(%d,%d)   \t%d\t%.1lf\n", Bigest_4->x, Bigest_4->y, Bigest_4->size, Bigest_4->intensity);
					if (Number_of_GAGs > 4)
					{
						printf_s("(%d,%d)   \t%d\t%.1lf\n", Bigest_5->x, Bigest_5->y, Bigest_5->size, Bigest_5->intensity);
					}
				}
			}
		}
	}
	printf_s("\n");
	return;
}
/*
The function get pointer for the array of GAG and the number of GAG we have in total.
It strar running from the first GAG in the list, calculete the distance from all the other GAGS, and
point on the closest 3. than move to the next GAG in the list and do the same. there is exaption that whan
GAG have a 3 other GAGS that point on him the function will skeep on him in the calculation and no more GAGS
can point on him.
the function return pointer to the first GAG.
*/
GAG*Map(GAG *arrayOFGAG, int Number_of_GAGs)
{
	GAG *temp_ptr = &arrayOFGAG[0];// all the pointers and the short_lanes in the GAG still didnt change from the begging.
	double distance = 0;//mesure the distance.
	double temp_distance = 0;//use to switch the ditance between 2 diffrent sohrt_line.
	for (int i = 0; i < Number_of_GAGs; i++)//nested loop run all over arrayofGAG.
	{
		for (int j = 0; j < Number_of_GAGs; j++)
		{
			if ((i == j) || (arrayOFGAG[j].counter == 3));//GAG cant point on himself and cant have more than 3 GAG that points on him.
			else
			{
				distance = sqrt(pow((arrayOFGAG[i].x - arrayOFGAG[j].x), 2) + pow((arrayOFGAG[i].y - arrayOFGAG[j].y), 2));//calcaulete the distance.
				if (distance < arrayOFGAG[i].short_lane3)//chack if its one of the 3 shortest path. short_lane3 is the longest in the list.
				{
					if (arrayOFGAG[i].ptr3)
					{
						arrayOFGAG[i].ptr3->counter--;//the GAG we chack is no longer point on this GAG. 
					}
					arrayOFGAG[i].ptr3 = &arrayOFGAG[j];
					arrayOFGAG[j].counter++;//new GAG point on him.
					arrayOFGAG[i].short_lane3 = distance;//save the ditance in the GAG in the same order as the pointer.
					if (distance < arrayOFGAG[i].short_lane2)
					{  //switch pointers and distances.
						temp_ptr = arrayOFGAG[i].ptr3;
						arrayOFGAG[i].ptr3 = arrayOFGAG[i].ptr2;
						arrayOFGAG[i].ptr2 = temp_ptr;
						temp_distance = arrayOFGAG[i].short_lane3;
						arrayOFGAG[i].short_lane3 = arrayOFGAG[i].short_lane2;
						arrayOFGAG[i].short_lane2 = temp_distance;
						if (distance < arrayOFGAG[i].short_lane1)
						{  //switch pointers and distances.
							temp_ptr = arrayOFGAG[i].ptr2;
							arrayOFGAG[i].ptr2 = arrayOFGAG[i].ptr1;
							arrayOFGAG[i].ptr1 = temp_ptr;
							temp_distance = arrayOFGAG[i].short_lane2;
							arrayOFGAG[i].short_lane2 = arrayOFGAG[i].short_lane1;
							arrayOFGAG[i].short_lane1 = temp_distance;
						}
					}
				}
			}
		}
	}
	printf_s("First map connection:\n");
	if (arrayOFGAG[0].ptr1 == NULL)
		printf_s("(%d,%d) no neighbours where found", arrayOFGAG[0].x, arrayOFGAG[0].y);
	else
	{
		printf_s("(%d,%d) neighbours are:", arrayOFGAG[0].x, arrayOFGAG[0].y);
		printf_s("(%d,%d),", arrayOFGAG[0].ptr1->x, arrayOFGAG[0].ptr1->y);
		if (arrayOFGAG[0].ptr2 != NULL)
		{
			printf_s("(%d,%d),", arrayOFGAG[0].ptr2->x, arrayOFGAG[0].ptr2->y);
			if (arrayOFGAG[0].ptr3 != NULL)
				printf_s("(%d,%d)\n", arrayOFGAG[0].ptr3->x, arrayOFGAG[0].ptr3->y);
		}
	}
	return arrayOFGAG;
}
/*
The function get the numbers of GAGs that we have.
Its chacking that the user input is integer number in the renge betwwen 1-number of GAGs.
Its putting the input in string and than chack char by char if they all valid for integer.
Its return the user input(just whan its valid, if not its run in loop until the input is good)
*/
int chack_input(int Number_of_GAGs)
{
	char chack_string[1000];//int is max 10 char.
	char testing = NULL;
	int i = 0;
	int input_length;
	char c;
	while ((c = getchar()) != '\n' && c != EOF);
	int input = 0;//scan number from the user.
	do
	{
		printf_s("\nSelect a number between 1 to %d: ", Number_of_GAGs);
		memset(chack_string, NULL, 1000); //set the string to NULL.
		i = 0;
		gets_s(chack_string, sizeof(chack_string));
		//fputc('\0', stdin);
		//fgets(chack_string, 1000, stdin);
		input_length = strlen(chack_string) - 1;
		input = atoi(chack_string);  //convert the chars to int.
		testing = 'T';  //T symbul for TRUE.
		if ((input_length > 10) || (input<1) || (input>Number_of_GAGs) || (chack_string[i] < '1') || (chack_string[i] > '9'))
			//the input must be between 1 to the number of GAG we have, and int is MAX 10 char so if the string is longer than that its doesnt a int .
			testing = 'F';  //F symbul for false.
		else
		{
			i++;
			while (i < input_length)
			{
				if (chack_string[i] < '0' || chack_string[i] > '9')// each char in int must be beetwen 0-9.
				{
					testing = 'F';
					break;  //if one char in the input is incorrect so its doesnt int for sure.
				}
				i++;
			}
		}
		if (testing == 'T')
			printf_s("The selected number is %d\n", input);// if after all the chacks its still 'T' so its a int
		else
		{
			printf_s("Input Err");// if after all the chacks its 'F' so its not a int.
		}
	} while (testing == 'F');// if the input is doesnt int its will run again until the input will be int.
	return input;
}
/* The function gets pointer to the array of GAGs,and pointer to the routh of the tree we build.
The root is the is the number of GAG the user inputs.
The function transform the map we built in "Map" function to a Tree.
It return pointer to the root of the tree.
*/
GAG *Tree(GAG *arrayOFGAG, GAG *Rott)
{
	char temp = arrayOFGAG->visit;
	GAG *hold;
	//pass the information from the GAG in the array to the GAG in the tree.
	Rott->x = arrayOFGAG->x;
	Rott->y = arrayOFGAG->y;
	Rott->short_lane1 = arrayOFGAG->short_lane1;
	Rott->short_lane2 = arrayOFGAG->short_lane2;
	Rott->short_lane3 = arrayOFGAG->short_lane3;
	Rott->visit = temp;
	arrayOFGAG->visit = 'v'; // Mark as visited.
							 //go shortest way.						
	if (!arrayOFGAG->ptr1 || arrayOFGAG->ptr1->visit == 'v')// end of branch
		Rott->ptr1 = NULL;
	else {
		hold = malloc(sizeof(GAG));//allocate memory to the node i the tree.
		Rott->ptr1 = Tree(arrayOFGAG->ptr1, hold);// recursion call.
												  //pass the information from the GAG in the array to the GAG in the tree.
		Rott->ptr1->x = arrayOFGAG->ptr1->x;
		Rott->ptr1->y = arrayOFGAG->ptr1->y;
		Rott->ptr1->short_lane1 = arrayOFGAG->ptr1->short_lane1;
		Rott->ptr1->short_lane2 = arrayOFGAG->ptr1->short_lane2;
		Rott->ptr1->short_lane3 = arrayOFGAG->ptr1->short_lane3;
	}
	//go middele distance
	if (!arrayOFGAG->ptr2 || arrayOFGAG->ptr2->visit == 'v')// end of branch
		Rott->ptr2 = NULL;
	else {
		hold = malloc(sizeof(GAG));
		Rott->ptr2 = Tree(arrayOFGAG->ptr2, hold);// recursion call.
												  //pass the information from the GAG in the array to the GAG in the tree.
		Rott->ptr2->x = arrayOFGAG->ptr2->x;
		Rott->ptr2->y = arrayOFGAG->ptr2->y;
		Rott->ptr2->short_lane1 = arrayOFGAG->ptr2->short_lane1;
		Rott->ptr2->short_lane2 = arrayOFGAG->ptr2->short_lane2;
		Rott->ptr2->short_lane3 = arrayOFGAG->ptr2->short_lane3;
	}
	//go longest way
	if (!arrayOFGAG->ptr3 || arrayOFGAG->ptr3->visit == 'v')// end of branch
		Rott->ptr3 = NULL;
	else {
		hold = malloc(sizeof(GAG));
		Rott->ptr3 = Tree(arrayOFGAG->ptr3, hold);// recursion call.
												  //pass the information from the GAG in the array to the GAG in the tree.
		Rott->ptr3->x = arrayOFGAG->ptr3->x;
		Rott->ptr3->y = arrayOFGAG->ptr3->y;
		Rott->ptr3->short_lane1 = arrayOFGAG->ptr3->short_lane1;
		Rott->ptr3->short_lane2 = arrayOFGAG->ptr3->short_lane2;
		Rott->ptr3->short_lane3 = arrayOFGAG->ptr3->short_lane3;
	}
	arrayOFGAG->visit = temp;
	return Rott;
}
GAG *Best_routh(GAG *Root, Shortest_route *ptr, double sum, GAG *arrayOFGAG)
{
	//Finding the shortest distance
	// Main idea: Run with 2 sums, outside sum and inside sum, find the minimum sum on the way to 1st GaG and update it to ptr->sum.
	// In order to declare that the distance is the shortest in the tree, the function must visit all the 1st GaGs in the tree.
	// After finishing, the minimum distance is achieved.
	if (Root != NULL)
	{
		if (Root->x == arrayOFGAG->x && Root->y == arrayOFGAG->y && Root->visit != 'b')
		{
			if (sum < ptr->sum)
			{
				ptr->sum = sum;
				ptr->count++;

			}
			Root->visit = 'b'; // Mark as visited.
			return ptr;
		}
		sum = sum + Root->short_lane1; //sum the branch
		Best_routh(Root->ptr1, ptr, sum, arrayOFGAG);
		sum = sum - Root->short_lane1; // Recursion reverse step(going 1 step outside , backwards) -> subtraction of the last branch visited .
		sum = sum + Root->short_lane2;
		Best_routh(Root->ptr2, ptr, sum, arrayOFGAG);
		sum = sum - Root->short_lane2;
		sum = sum + Root->short_lane3;
		Best_routh(Root->ptr3, ptr, sum, arrayOFGAG);
		sum = sum - Root->short_lane3;
	}

}
LC *Map_route(GAG *Root, Shortest_route *ptr, double sum, GAG *arrayOFGAG, LC *link, LC *bighead)
{
	//Main idea: Mapping the route forward by reconstrucing a linked list forward and backward inside-out recursivly until objective has found.
	//Way: Using the same way of function "Best_Route" in order to sum the distances, but this time preserving the data of the route.
	// The function allocate memory when going Inside-Step and free the memory when going Outside-Step (Reversing)
	// Once the shortest distance reached, 2 flags going ON (count and flag) and the function no longer free the memory.
	double epsilon = 0.01;
	if (Root != NULL)
	{
		if (Root->x == arrayOFGAG->x && Root->y == arrayOFGAG->y && Root->visit == 'b')
		{
			if (abs(sum - ptr->sum)<epsilon) // comparing between 2 double numbers.
			{
				ptr->count++; // flag - closing the recursion
				ptr->flag = 1; // flag - closing the recursion
				link->x = Root->x;
				link->y = Root->y;
				link->next = NULL;
			}
			Root->visit = 'v'; // Mark as visited.
		}
		if (ptr->flag == 0)
		{
			if (ptr->count > 0)
				return;
			link->next = malloc(sizeof(LC)); // Allocate memory
			link->x = Root->x; // Update data
			link->y = Root->y;
			sum = sum + Root->short_lane1; // Step inside - Sum.
			Map_route(Root->ptr1, ptr, sum, arrayOFGAG, link->next, bighead);
			sum = sum - Root->short_lane1; // Step Outside - Subtract the branch.
			if (ptr->count == 0) // Free memory only when flag is OFF.
				free(link->next);
			link->x = Root->x;
			link->y = Root->y;
			if (ptr->count > 0)
				return;
			link->next = malloc(sizeof(LC));
			sum = sum + Root->short_lane2;
			Map_route(Root->ptr2, ptr, sum, arrayOFGAG, link->next, bighead);
			sum = sum - Root->short_lane2;
			if (ptr->count == 0) // Free memory only when flag is OFF.
				free(link->next);
			link->x = Root->x;
			link->y = Root->y;
			if (ptr->count > 0)
				return;
			link->next = malloc(sizeof(LC));
			sum = sum + Root->short_lane3;
			Map_route(Root->ptr3, ptr, sum, arrayOFGAG, link->next, bighead);
			sum = sum - Root->short_lane3;
			if (ptr->count == 0) // Free memory only when flag is OFF.
				free(link->next);
		}
	}
}
/*
The function make a copy to the bmp file(d2).
*/
void Create_BMP()
{
	FILE * original;// open the bmp picture.
	FILE * copy;// use to copy the data for a new picture.
	unsigned int get_char;
	fopen_s(&original, "d2.bmp", "rb");// open the original bmp.
	fopen_s(&copy, "dustcopy.bmp", "wb");//open new bmp. 
	while ((get_char = fgetc(original)) != EOF)//run until we finih scan all the original bmp.
	{
		fputc(get_char, copy);// copy the data we scan from the original file to the new one.
	}
	fclose(original);
	fclose(copy);
	printf_s("\n\nCheck for a new copy file dustcopy.bmp\n");
	return;
}
/*
The function get pointer to the shortest way, matrix that simulate the bmp stracture(the pixels area)
and to double, m is the slope line and c is the const(y=mx+c).
Its run from one GAG to another and chacks the relevant pixel in compere to the line eqution.
every step the pixel that is the closest is marked and its move toward him and do it until its
arrive to the second GAG.in the end we have in the matrix all the pixels we need to paint
mark by the number 1.
*/
void mark_the_shorter_path_in_matrix(LC *ptrGAG, int **matrix, double m, double c)
{
	double d1 = 0, d2 = 0, d3 = 0;
	if ((ptrGAG->next->x == ptrGAG->x) && (ptrGAG->next->y == ptrGAG->y))//until we arrive to the second GAG.
		return;
	if (m == 9999999999)// liner equation(x=const)
	{
		// moving just in the y axis.
		if (ptrGAG->next->y > ptrGAG->y)
		{
			//going up in the y axis until reached the GAG.
			ptrGAG->y++;
			matrix[ptrGAG->y][ptrGAG->x] = 1;
			mark_the_shorter_path_in_matrix(ptrGAG, matrix, m, c);
		}
		else
		{
			//going down in the y axis until reached the GAG.
			ptrGAG->y--;
			matrix[ptrGAG->y][ptrGAG->x] = 1;
			mark_the_shorter_path_in_matrix(ptrGAG, matrix, m, c);
		}
	}
	else if (m > 0)// liner equation(y=mx+c)
	{
		if (ptrGAG->next->y > ptrGAG->y)
			// moving up in the y axis and right in the x axis.
		{
			//chack the distance from the 3 relevant pixels.
			d1 = (fabs(m*((double)ptrGAG->x + 1) - (double)ptrGAG->y + c));
			d2 = (fabs(m*((double)ptrGAG->x) - (double)ptrGAG->y - 1 + c));
			d3 = (fabs(m*((double)ptrGAG->x + 1) - (double)ptrGAG->y - 1 + c));
			//chose the shortest way, promote who we need and mark the pixel in the matrix.
			if ((d1 <= d2) && (d1 <= d3))
			{
				ptrGAG->x++;
				matrix[ptrGAG->y][ptrGAG->x] = 1;
				mark_the_shorter_path_in_matrix(ptrGAG, matrix, m, c);
			}
			else if ((d2 <= d1) && (d2 <= d3))
			{
				ptrGAG->y++;
				matrix[ptrGAG->y][ptrGAG->x] = 1;
				mark_the_shorter_path_in_matrix(ptrGAG, matrix, m, c);
			}
			else if ((d3 <= d1) && (d3 <= d2))
			{
				ptrGAG->y++;
				ptrGAG->x++;
				matrix[ptrGAG->y][ptrGAG->x] = 1;
				mark_the_shorter_path_in_matrix(ptrGAG, matrix, m, c);
			}
		}
		else if (ptrGAG->next->y < ptrGAG->y) // moving down in the y axis and left in the x axis.
		{
			//chack the distance from the 3 relevant pixels.
			d1 = (fabs(m*((double)ptrGAG->x - 1) - (double)ptrGAG->y + c));
			d2 = (fabs(m*((double)ptrGAG->x) - (double)ptrGAG->y + 1 + c));
			d3 = (fabs(m*((double)ptrGAG->x - 1) - (double)ptrGAG->y + 1 + c));
			//chose the shortest way, promote who we need and mark the pixel in the matrix.
			if ((d1 <= d2) && (d1 <= d3))
			{
				ptrGAG->x--;
				matrix[ptrGAG->y][ptrGAG->x] = 1;
				mark_the_shorter_path_in_matrix(ptrGAG, matrix, m, c);
			}
			else if ((d2 <= d1) && (d2 <= d3))
			{
				ptrGAG->y--;
				matrix[ptrGAG->y][ptrGAG->x] = 1;
				mark_the_shorter_path_in_matrix(ptrGAG, matrix, m, c);
			}
			else if ((d3 <= d1) && (d3 <= d2))
			{
				ptrGAG->y--;
				ptrGAG->x--;
				matrix[ptrGAG->y][ptrGAG->x] = 1;
				mark_the_shorter_path_in_matrix(ptrGAG, matrix, m, c);
			}
		}
		else if (ptrGAG->next->y == ptrGAG->y)
		{
			if (ptrGAG->next->x > ptrGAG->x)
				//going right in the x axis until reached the GAG.
			{
				ptrGAG->x++;
				matrix[ptrGAG->y][ptrGAG->x] = 1;
				mark_the_shorter_path_in_matrix(ptrGAG, matrix, m, c);
			}
			else
				//going left in the x axis until reached the GAG.
			{
				ptrGAG->x--;
				matrix[ptrGAG->y][ptrGAG->x] = 1;
				mark_the_shorter_path_in_matrix(ptrGAG, matrix, m, c);
			}
		}
	}
	else if (m < 0)// liner equation(y=-mx+c)
	{
		if (ptrGAG->next->y > ptrGAG->y)
			// moving up in the y axis and left in the x axis.
		{
			//chack the distance from the 3 relevant pixels.
			d1 = (fabs(m*((double)ptrGAG->x - 1) - (double)ptrGAG->y + c));
			d2 = (fabs(m*((double)ptrGAG->x) - (double)ptrGAG->y - 1 + c));
			d3 = (fabs(m*((double)ptrGAG->x - 1) - (double)ptrGAG->y - 1 + c));
			//chose the shortest way, promote who we need and mark the pixel in the matrix.
			if ((d1 <= d2) && (d1 <= d3))
			{
				ptrGAG->x--;
				matrix[ptrGAG->y][ptrGAG->x] = 1;
				mark_the_shorter_path_in_matrix(ptrGAG, matrix, m, c);
			}
			else if ((d2 <= d1) && (d2 <= d3))
			{
				ptrGAG->y++;
				matrix[ptrGAG->y][ptrGAG->x] = 1;
				mark_the_shorter_path_in_matrix(ptrGAG, matrix, m, c);
			}
			else if ((d3 <= d1) && (d3 <= d2))
			{
				ptrGAG->y++;
				ptrGAG->x--;
				matrix[ptrGAG->y][ptrGAG->x] = 1;
				mark_the_shorter_path_in_matrix(ptrGAG, matrix, m, c);
			}
		}
		else if (ptrGAG->next->y < ptrGAG->y) // moving down in the y axis and right in the x axis.
		{
			//chack the distance from the 3 relevant pixels.
			d1 = fabs((m*((double)ptrGAG->x + 1) - (double)ptrGAG->y + c));
			d2 = fabs((m*((double)ptrGAG->x) - (double)ptrGAG->y + 1 + c));
			d3 = fabs((m*((double)ptrGAG->x + 1) - (double)ptrGAG->y + 1 + c));
			//chose the shortest way, promote who we need and mark the pixel in the matrix.
			if ((d1 <= d2) && (d1 <= d3))
			{
				ptrGAG->x++;
				matrix[ptrGAG->y][ptrGAG->x] = 1;
				mark_the_shorter_path_in_matrix(ptrGAG, matrix, m, c);
			}
			else if ((d2 <= d1) && (d2 <= d3))
			{
				ptrGAG->y--;
				matrix[ptrGAG->y][ptrGAG->x] = 1;
				mark_the_shorter_path_in_matrix(ptrGAG, matrix, m, c);
			}
			else if ((d3 <= d1) && (d3 <= d2))
			{
				ptrGAG->y--;
				ptrGAG->x++;
				matrix[ptrGAG->y][ptrGAG->x] = 1;
				mark_the_shorter_path_in_matrix(ptrGAG, matrix, m, c);
			}
		}
		else if (ptrGAG->next->y == ptrGAG->y)
		{
			if (ptrGAG->next->x > ptrGAG->x)
				//going right in the x axis until reached the GAG.
			{
				ptrGAG->x++;
				matrix[ptrGAG->y][ptrGAG->x] = 1;
				mark_the_shorter_path_in_matrix(ptrGAG, matrix, m, c);
			}
			else
				//going left in the x axis until reached the GAG.
			{
				ptrGAG->x--;
				matrix[ptrGAG->y][ptrGAG->x] = 1;
				mark_the_shorter_path_in_matrix(ptrGAG, matrix, m, c);
			}
		}
	}
	else if (m == 0)// liner equation(y=const)
					// moving just in the x axis.
	{
		if (ptrGAG->next->x > ptrGAG->x)
			//going right in the x axis until reached the GAG.
		{
			ptrGAG->x++;
			matrix[ptrGAG->y][ptrGAG->x] = 1;
			mark_the_shorter_path_in_matrix(ptrGAG, matrix, m, c);
		}
		else
			//going left in the x axis until reached the GAG.
		{
			ptrGAG->x--;
			matrix[ptrGAG->y][ptrGAG->x] = 1;
			mark_the_shorter_path_in_matrix(ptrGAG, matrix, m, c);
		}
	}
	return;
}
/*
The function get the marked matrix from "mark_the_shorter_path_in_matrix" function.
Its run all over the picture and paint in read every placed that marked in the matrix.
in the end in the picture we see a red stright line that connect between every 2
GAGs that are in the shortest way.
*/
void print_path_by_matrix(int **matrix,int wpix,int hpix)
{
	int padding = 0;
	padding = (4 - (wpix * 3) % 4) % 4;
	unsigned int offset = 0;
	FILE *p2; // points on the first char in the bmp file
	fopen_s(&p2, "dustcopy.bmp", "r+");//opening the copy of d2 picture called dustcopy
	fseek(p2, 10, SEEK_SET);
	for (int i = 0; i < 4; i++)
		offset |= fgetc(p2) << (i * 8);
	fseek(p2, offset, SEEK_SET);// starting from the information about the pixels
	for (int i = 0; i < hpix; i++)
	{// nested loop which scan the all picture and change the pixels information.	
		for (int j = 0; j < wpix; j++)
		{
			if (matrix[i][j] == 1)// marked places in the matrix
			{
				// paint the pixel in red.
				fputc(0, p2); // B
				fputc(0, p2); // G
				fputc(255, p2); // R
			}
			else
			{
				fseek(p2, 3, SEEK_CUR);// starts to advance from pixel to another (3 bytes- r+g+b) from the start of the file (SEEK_CUR)

			}
		}
		fseek(p2, padding, SEEK_CUR);
	}
	fclose(p2);
	printf_s("Done connecting bigest dust particles in dustcopy.bmp\n");
}
/*
The function gets pointer to the array of GAGS,the number of GAGs and pointer to a matrix which simulate the bmp stracture(the pixels area).
It puts a green plus(5x5 pixels) in every GAG cordination according to the matrix.
*/
void prints_pluses_on_the_GAGs(GAG* arrayofGAG, int numOfGAGs, int **matrix,int wpix,int hpix)
{
	int padding = 0;
	padding = (4 - (wpix * 3) % 4) % 4;
	unsigned int offset = 0;
	FILE *p2; // points on the first char in the bmp file
	for (int i = 0; i < hpix; i++)// run over all the matrix in nested loop and puts zeros
	{
		for (int j = 0; j < wpix; j++)
		{
			matrix[i][j] = 0;
		}
	}
	for (int p = 0; p < numOfGAGs; p++)
		/*puts 1 in the matrix in shape of plus(5x5) that the center is the GAG cordination, do that for all the GAGs.  */
	{
		matrix[arrayofGAG[p].y][arrayofGAG[p].x] = 1;//center of plus
													 //const x
		matrix[arrayofGAG[p].y - 1][arrayofGAG[p].x] = 1;
		matrix[arrayofGAG[p].y - 2][arrayofGAG[p].x] = 1;
		matrix[arrayofGAG[p].y + 1][arrayofGAG[p].x] = 1;
		matrix[arrayofGAG[p].y + 2][arrayofGAG[p].x] = 1;
		// const y
		matrix[arrayofGAG[p].y][arrayofGAG[p].x - 1] = 1;
		matrix[arrayofGAG[p].y][arrayofGAG[p].x - 2] = 1;
		matrix[arrayofGAG[p].y][arrayofGAG[p].x + 1] = 1;
		matrix[arrayofGAG[p].y][arrayofGAG[p].x + 2] = 1;
	}
	fopen_s(&p2, "dustcopy.bmp", "r+");//opening the copy of d2 picture called "dustcopy"
	fseek(p2, 10, SEEK_SET);
	for (int i = 0; i < 4; i++)
		offset |= fgetc(p2) << (i * 8);
	fseek(p2, offset, SEEK_SET);// starting from the information about the pixels
	for (int i = 0; i < hpix; i++)//runs on y cordination
	{// nested loop which scan the all picture and change the pixels information (puts in every GAG's cordination a green plus)
	 //according to the matrix.
		for (int j = 0; j < wpix; j++)// runs on x_cordination
		{
			if (matrix[i][j] == 1)
			{
				// changes the color to green in every place that the matrix value is 1.
				fputc(0, p2); // B
				fputc(255, p2); // G
				fputc(0, p2); // R
			}
			else
			{
				fseek(p2, 3, SEEK_CUR);// starts to advance from pixel to another (3 bytes- r+g+b) from the start of the file (SEEK_CUR)

			}
		}
		fseek(p2, padding, SEEK_CUR);
	}
	fclose(p2);
	printf_s("Done marking coordinate of bigest dust particles in dustcopy.bmp\n");
	printf_s("\n");
	return;
}
/*
The function get pointer to the routh of the tree ane deallocate the tree.*/
void deallocate_tree(GAG *Root)
{
	if (Root == NULL)
		return;
	deallocate_tree(Root->ptr1);
	deallocate_tree(Root->ptr2);
	deallocate_tree(Root->ptr3);
	free(Root);
}
/* 
this function calculate the chances that a random pixel found inside a GAG*/
void ChancesForGAG(GAG*  arrayofGAG, int Number_of_GAGs,int wpix,int hpix)
{
	int SumOfSizes, i = 0, j = 0;
	SumOfSizes = 0;
	int firstGAG_x = 0;
	unsigned int SizeOfPic = 0;
	double ChanceFor_GAG = 0;
	for (i = 0; i < Number_of_GAGs; i++)
	{
		SumOfSizes = arrayofGAG[i].size;//the amount of pixels that are part of a GAG
	}
	// using the first part variables for the picture wide and length

	SizeOfPic = wpix * hpix;// this information is in the main so we can use in the function
	ChanceFor_GAG = ((double)SumOfSizes / (double)SizeOfPic) * 100;
	printf_s("\n");
	printf_s("The chances to choose a pixel which belongs to a GAG from the picture is %.3lf%%\n", ChanceFor_GAG);
	printf_s("\n");
}
/*
this function outputs 4 areas, and information about each area - their location, amount of GAGs in each area, their place compare two the others areas.*/  
void minMaxGAGs(GAG* arrayofGAG, int Number_of_GAGs, dataforprint *ptr)
{
	double gap = 0, firstGAG_x = 1000000, lastGAG_x = 0;
	int area1 = 0, area2 = 0, area3 = 0, area4 = 0;
	int total_size1 = 0, total_size2 = 0, total_size3 = 0, total_size4 = 0;
	int temp;
	int *array1;// this array contains the amount of GAGs in each area --> we will sort this array from min to max
	double first_x_cordination;

	for (int i = 0; i < Number_of_GAGs; i++)// finding the first and the last GAG with the smallest and the biggest x cordanation
	{

		if (arrayofGAG[i].x > lastGAG_x)
			lastGAG_x = arrayofGAG[i].x;
		else if (arrayofGAG[i].x < firstGAG_x)
			firstGAG_x = arrayofGAG[i].x;
	}

	gap = (lastGAG_x - firstGAG_x) / 4;// the gap between two neighbor areas
	for (int j = 0; j < Number_of_GAGs; j++)// chacking the x cordination of each GAG between the first GAG to the last one 
	{// chacking for each GAG in which area it found (according to it's x cordination)
		if (arrayofGAG[j].x >= firstGAG_x && arrayofGAG[j].x < (firstGAG_x + gap))
		{
			total_size1 += arrayofGAG[j].size;
			area1++;
		}


		else if (arrayofGAG[j].x >(firstGAG_x + gap) && arrayofGAG[j].x < (firstGAG_x + 2 * gap))
		{
			total_size2 += arrayofGAG[j].size;
			area2++;
		}
		else if (arrayofGAG[j].x >(firstGAG_x + 2 * gap) && arrayofGAG[j].x < (firstGAG_x + 3 * gap))
		{
			total_size3 += arrayofGAG[j].size;
			area3++;
		}
		else if (arrayofGAG[j].x >(firstGAG_x + 3 * gap) && arrayofGAG[j].x <= lastGAG_x)
		{
			total_size4 += arrayofGAG[j].size;
			area4++;
		}
	}

	array1 = (int)(malloc(sizeof(int) * 4));// puts the amount of gags for each area in array
	if (array1 == NULL) exit(2);//chacking if the malloc operation succeed
	array1[0] = total_size1;
	array1[1] = total_size2;
	array1[2] = total_size3;
	array1[3] = total_size4;


	/*   Bubble sorting begins */

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < (3 - i); j++)
		{
			if (array1[j] > array1[j + 1])
			{
				temp = array1[j];
				array1[j] = array1[j + 1];
				array1[j + 1] = temp;//now array1 sorted - from minumum size of GAGs to maximun size 
			}
		}
	}
	printf_s("The sizes of the GAGs (in pixels) in each area organized from the minimum to maximum are: \n");
	for (int i = 0; i < 4; i++)
	{
		printf_s("%d ", array1[i]);
	}
	printf_s("\n\n");

	// prints the areas' cordantion, the amount of GAGs in each one and the amount of pixels that are GAGs.
	printf_s("Area 1 between  %d to %d  x cordination has %d GAGs and size of %d pixels\n", (int)firstGAG_x, (int)(firstGAG_x + gap), area1, total_size1);
	printf_s("Area 2 between  %d to %d  x cordination has %d GAGs and size of %d pixels\n", (int)(firstGAG_x + gap), (int)(firstGAG_x + 2 * gap), area2, total_size2);
	printf_s("Area 3 between  %d to %d  x cordination has %d GAGs and size of %d pixels\n", (int)(firstGAG_x + 2 * gap), (int)(firstGAG_x + 3 * gap), area3, total_size3);
	printf_s("Area 4 between  %d to %d  x cordination has %d GAGs and size of %d pixels\n", (int)(firstGAG_x + 3 * gap), (int)lastGAG_x, area4, total_size4);
	printf_s("\n");
	// chacking which area got the maximum GAGs

	if (array1[3] == total_size1)
	{
		printf_s("The area with the maximum size of GAGs is area 1!\n");
		first_x_cordination = firstGAG_x;
	}
	else if (array1[3] == total_size2)
	{
		printf_s("The area with the maximum size of GAGs is area 2!\n");
		first_x_cordination = firstGAG_x + gap;
	}
	else if (array1[3] == total_size3)
	{
		printf_s("The area with the maximum size of GAGs is area 3!\n");
		first_x_cordination = firstGAG_x + 2 * gap;
	}
	else if (array1[3] == total_size4)
	{
		printf_s("The area with the maximum size of GAGs is area 4!\n");
		first_x_cordination = firstGAG_x + 3 * gap;
	}
	if (total_size1 == 0 && total_size2 == 0 && total_size3 == 0 && total_size4 == 0)
		printf_s("There is no GAG's in the picture.\n");
	printf_s("\n");
	ptr->x = first_x_cordination;// with this two pointers we will be able two use this information in "print_area_with_max_GAGS" function 
	ptr->gap = gap;

	free(array1);
}
/*
this function coping the original picture (d2) for use this copy in function "print_area_with_GAGs". the name of the new pic is "maxGAGarea"*/
void Create_BMP_max()
{
	FILE * bSrc;
	FILE * bCpy;
	unsigned int get_char;
	fopen_s(&bSrc, "d2.bmp", "rb");
	fopen_s(&bCpy, "maxGAGarea.bmp", "wb");
	while ((get_char = fgetc(bSrc)) != EOF)
	{
		fputc(get_char, bCpy);
	}
	fclose(bSrc);
	fclose(bCpy);
	return;
}
/*
this fuction mark the area with the most GAGs with two verticales red lines*/
void print_area_with_max_GAGS(dataforprint *ptr,int wpix,int hpix)
{
	int padding = 0;
	padding = (4 - (wpix * 3) % 4) % 4;
	unsigned int offset = 0;
	FILE *p1; // points on the first char in the bmp file
	fopen_s(&p1, "maxGAGarea.bmp", "r+");//opening the copy of d2 picture called maxGAGarea
	fseek(p1, 10, SEEK_SET);
	for (int i = 0; i < 4; i++)
		offset |= fgetc(p1) << (i * 8);
	fseek(p1, offset, SEEK_SET);// starting from the information about the pixels
	for (int i = 0; i < hpix; i++)
		// nested loop which scan the all picture and change the pixels information (leave only the reds pixels and puts zero in the blue and green pixels)
	{
		for (int j = 0; j < wpix; j++)
		{//gets from "minMaxGAGs" function pointer to stract dataforprint which contain the information about x_cordination and the gap.
			if (j == (ptr->x) || j == ((ptr->x) + (ptr->gap)))
			{
				// changes the color of two columns which their x_cordinations are the cordinations of the area with maximum GAGs
				fputc(0, p1); // B
				fputc(0, p1); // G
				fputc(255, p1); // R

			}
			else
			{
				fseek(p1, 3, SEEK_CUR);// starts to advance from pixel to another (3 bytes- r+g+b) from the start of the file (SEEK_CUR)

			}
		}
		fseek(p1, padding, SEEK_CUR);
	}
	fclose(p1);// close the file
}