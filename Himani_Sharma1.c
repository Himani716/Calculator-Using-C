#include<stdio.h>
#include<conio.h>
#include<math.h>
#include<stdlib.h>
#include"Himani_Sharma4.h" 
#include"Himani_Sharma2.h"
#include"Himani_Sharma3.h"
int Mathematics_Calculator();
int Statistical_Calculator();
int Numerical_calculator();
int main()
{
    int num;
    while(1){
	 
    printf("\n 1) Mathematics Calculator");
    printf("\n 2) Statistical Calculator");
    printf("\n 3) Numerical Calculator");
    printf("\n 4) exit");
    printf("\n enter a number:");
    scanf("%d",&num);
    switch(num)
    {
               case 1:
                    printf("\n Mathematical Calculations :");
                    Mathematics_Calculator();
                    break;
               case 2 :
                    printf("\n Statistical Calculations :");
                    Statistical_Calculator();
                    break;
               case 3 :
                    printf("\n Numerical Calculations :");
                    Numerical_calculator();
                    break;
               case 4 :
                    printf("\n Exiting the calculator" );
                    return 0;
                default:
                	printf("\n Please choose a different option");
                	break;
    }
}
}
int Mathematics_Calculator(){
	int choice = 0;
	double num,result,n,r;
	int x,y,hcf,s,t,lcm;
	while(1)
	{
		printf("\n ---Mathematics Contents---\n");
		printf("\n 1) Addition");
		printf("\n 2) Subraction");
		printf("\n 3) Multiplication");
		printf("\n 4) Division");
		printf("\n 5) Exponent");
		printf("\n 6) Factorial");
		printf("\n 7) Permutation");
		printf("\n 8) Combination");
		printf("\n 9) HCF");
		printf("\n 10) LCM");
		printf("\n 11) Simple Interest");
		printf("\n 12) Compound Interest");
		printf("\n 13) Create AP");
		printf("\n 14) Create HP");
		printf("\n 15) Quadratic Roots");
		printf("\n 16) Exit");
		printf("\n Enter your choice :");
		scanf("%d",&choice);// select any one of the above options
	
	switch(choice){
	case 1: Addition(); break;
	case 2: Subraction(); break;
	case 3: Multiplication(); break;
	case 4: Division(); break;
	case 5: Exponent(); break;
    case 6:
    	printf("\n Enter a number");
    	scanf("%lf",&num);
    	if (num<0)
    	printf("\n Factorial is not defined for negative numbers");
    	else {
    		result = Factorial(num);
    		printf("\n----------------------------");
    		printf("\nFactorial of %.2lf = %.2lf\n", num, result);
    		printf("\n----------------------------");
    }
		break;
    
	case 7: 
		printf("\n Enter n and r :");
		scanf("%lf %lf",&n,&r);
		double n_fact,n_r,n_r_fact,npr;
		n_fact = Factorial(n);
		n_r = n-r;
		n_r_fact = Factorial(n_r);
		npr = n_fact/n_r_fact;
		printf("\n----------------------------");
		printf("\n %lfP%lf is %lf ",n,r,npr);
		printf("\n----------------------------");
		break;
	case 8: 
		printf("\n Enter n and r :");
		scanf("%lf %lf",&n,&r);
		double n_fact1,n_r1,n_r_fact1,ncr,r_fact;
		n_fact1 = Factorial(n);
		r_fact = Factorial(r);
		n_r1 = n-r;
		n_r_fact1 = Factorial(n_r1);
		ncr = n_fact1/(n_r_fact1*r_fact);
		printf("\n----------------------------");
		printf("\n %lfC%lf is %lf ",n,r,ncr);
		printf("\n----------------------------");
		break;
	    case 9:
    	
		printf("\n enter two numbers :");
		scanf("%d %d",&x,&y);
		hcf = HCF(x,y);
		printf("\n----------------------------");
		printf("\n HCF of %d and %d is %d",x,y,hcf);
		printf("\n----------------------------");
		break;
		
	    case 10:
		//int s,t,lcm;
		printf("\n enter two numbers :");
		scanf("%d %d",&s,&t);
		lcm = LCM(s,t);
		printf("\n----------------------------");
		printf("\n LCM of %d and %d is %d",s,t,lcm);
		printf("\n----------------------------"); break;
	    case 11: Simple_Interest(); break;
	    case 12: Compound_interest(); break;
	    case 13: Arithmetic(); break;
	    case 14: Harmonic(); break;
	    case 15 : Quad_roots(); break;
	    case 16: printf("\nExited Mathematical section\n"); return 0;
	    default:
	    printf("\n Invalid choice , Enter a valid number");  
	    break;
	}
}
}

int Numerical_calculator(){
	int pick = 0;
		while(1)
	{
		printf("\n ---Numerical Contents---\n");
		printf("\n 1)Bisection Method");
		printf("\n 2)Secant Method");
		printf("\n 3)Regula Falsi Method");
		printf("\n 4)Newton Raphson Method");
		printf("\n 5)Muller Method");
		printf("\n 6)Chebyshev Method");
		printf("\n 7)Bairstow Method");
		printf("\n 8)Trapezoidal Method");
		printf("\n 9)Simpson1by3 Method");
		printf("\n 10)Simpson3by8 Method");
		printf("\n 11)Birge Vieta Method");
		printf("\n 12)Forward Substitution Method");
		printf("\n 13)Backward Subsitution");
		printf("\n 14)Jacobi Iteration Method");
		printf("\n 15)Lagrange interpolation Method");
		printf("\n 16)Hermite Interpolation Method");
		printf("\n 17) Exit");
		printf("\n Enter your choice :");
		scanf("%d",&pick);// select any one of the above options
		
	switch(pick){
	case 1: Bisection_Method(); break;
	case 2: Secant_Method(); break;
	case 3: Regula_Falsi(); break;
	case 4: Newton_Raphson_Method(); break;
	case 5: Muller_Method(); break;	
	case 6: Chebyshev();  break;
	case 7: Bairstow_Method(); break;
	case 8: Trapezoidal(); break;
	case 9 : Simpson1by3(); break;
	case 10 : Simpson3by8(); break;
	case 11: Birge_Vieta(); break;
	case 12: Forward_Subtitution(); break;
	case 13 : Backward_Subsitution(); break;
	case 14: Jacobi_Iteration(); break;
	case 15 :Lagrange_interpolation();
	case 16 :Hermite_IP();
	case 17 : printf("\nExited Numerical section\n"); return 0;
	default:
	    printf("\n Invalid choice , Enter a valid number");	
	    break;
}
}
}
int Statistical_Calculator(){
	int choose = 0;
		while(1)
	{
		printf("\n ---Statistical Contents---\n");
		printf("\n 1) Mean");
		printf("\n 2) Variance");
		printf("\n 3) Standard Deviation ");
		printf("\n 4) Range");
		printf("\n 5) Mode");
		printf("\n 6) Skewness and kurtosis");
		printf("\n 7) Coefficient Of variation:");
		printf("\n 8) Geometric Mean");
		printf("\n 9) Harmonic Mean");
		printf("\n 10) Exit")	;
		printf("\n Enter your choice :");
		scanf("%d",&choose);// select any one of the above options
		
	switch(choose){
		case 1: Mean();break;
		case 2: Variance(); break;
		case 3: Std_deviation(); break;
		case 4: Range(); break;
		case 5 : Mode(); break;
		case 6 : Skewness_Kurtosis(); break;
		case 7 : CV(); break;
		case 8: GeometricMean(); break;
		case 9:HarmonicMean(); break;
		case 10: printf("\nExited Statistical section\n"); return 0;
		default:
	    printf("\n Invalid choice , Enter a valid number");
	    break;
		}
	}
}


