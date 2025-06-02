




int Addition(){
	    double sum;
	    int n,i;
	    printf("\n enter Size of the array :");
	    scanf("%d",&n);
	    double number[n];
	    for (i=0;i<n;i++){
		    printf("\n enter number:");
		    scanf("%lf",&number[i]);
	    }
	    for (i=0;i<n;i++)
		    sum = sum + number[i];
		printf("\n----------------------------");
	    printf("\n Sum = %lf\n",sum);
	    printf("----------------------------");
}
int Subraction(){
	    double diff;
	    int n,i;
	    printf("\n enter Size of the array :");
	    scanf("%d",&n);
	    double number[n];
	    for (i=0;i<n;i++){
		    printf("\n enter number:");
		    scanf("%lf",&number[i]);
	    }
	    diff = number[0];
	    for (i=1;i<n;i++)//The subtraction loop should start from the second element (i = 1) 
		    diff = diff - number[i];//after assigning the first element to diff.
		printf("\n----------------------------");
	    printf("\n Subraction = %lf\n",diff);
	    printf("----------------------------");
}
int Multiplication(){
	    double prod ;
	    int n,i;
	    printf("\n enter Size of the array :");
	    scanf("%d",&n);
	    double number[n];
	    for (i=0;i<n;i++){
		    printf("\n enter number:");
		    scanf("%lf",&number[i]);
	    }
	    prod = number[0];
	    for (i=1;i<n;i++)
		    prod = prod * number[i];
		printf("\n----------------------------");
	    printf("\n Multiplication = %lf\n",prod);
	    printf("----------------------------");
}
int Division(){
	    double div ;
	    int n,i;
	    printf("\n enter Size of the array :");
	    scanf("%d",&n);
	    double number[n];
	    if (n==0)
	    	printf("\n Cannot be divided by zero; enter another number");
	    else {
	    for (i=0;i<n;i++){
		    printf("\n enter number:");
		    scanf("%lf",&number[i]);
	    }
	    div = number[0];
	    for (i=1;i<n;i++)
		    div = div/number[i];
	}
		printf("\n----------------------------");
	    printf("\n Division = %lf\n",div);
		printf("-----------------------------");
}
int Exponent(){
	double x,a,power=1;
	int i;
	printf("\n Enter tha value of x and a:");
	scanf("%lf %lf",&x,&a);
	for (i=1;i<=a;i++){
		power = power*x;
	}
	printf("\n----------------------------");
	printf("\n x raise to power a is %lf ",power);
	printf("\n-----------------------------");
}
int Factorial(double x) {
    if (x == 0 || x == 1) // Base case
        return 1;
    else 
        return x * Factorial(x - 1); // Recursive step
}
int HCF(int a, int b) {
    if (b == 0)
        return a;
    else
        return HCF(b, a % b);
}
int LCM(int a, int  b) {
    return (a * b) / HCF(a, b);
}
int Harmonic(void)
{
    float a , d , term , sum = 0;
    int i , n;

    printf("\nEnter a , d , n\n");
    scanf("%f %f %d" , &a , &d , &n);

    for (i = 1 ; i <= n ; i++)
    {
        term = 1.0 / (a + ((i - 1) * d));
        printf("%f " , term);

        sum = sum + term;
    }
    printf("\n%f\n" , sum);
}
int Arithmetic(void)
{
    float a , d , term ;
    int i , n;

    printf("\nEnter a , d , n\n");
    scanf("%f %f %d" , &a , &d , &n);

    for (i = 1 ; i <= n ; i++)
    {
        term = (a + ((i - 1) * d));
        printf("%f " , term);
    }
}
int Simple_Interest()
{
	int t;
	float r,si,p;
	printf("\nEnter Principal Amount (P): ");
	scanf("%f",&p);
	printf("\nEnter Rate of Interest (R): ");
	scanf("%f",&r);
	printf("\nEnter Time (T): ");
	scanf("%d",&t);
	si = (p*t*r)/100.0;
	printf("\n-----------------------------");
	printf("\n Simple interest = %f%",si);
	printf("\n-----------------------------");
	getch();
}
int Compound_interest(){
	int t;
	float r,ci,amount,p;
	printf("\nEnter Principal Amount (P): ");
	scanf("%f",&p);
	printf("\nEnter Rate of Interest (R): ");
	scanf("%f",&r);
	printf("\nEnter Time (T): ");
	scanf("%d",&t);
	amount = p * pow(1 + r/100,t);
	ci = amount - p;
	printf("\n-----------------------------");
	printf("\n Total Amount  = %f",amount);
	printf("\n-----------------------------");
	printf("\n Compound interest = %f%",ci);
	printf("\n-----------------------------");
}
int Quad_roots()
{
    int a,b,c;
    float d,x,y,real_part,imag_part;
    printf("enter the value of coefficients:");
    scanf("%d%d%d",&a,&b,&c);
    d = (b*b)-(4*a*c);
    if (d >=0)
    {
        x = (-b + sqrt(d))/(2*a);
        y = (-b - sqrt(d))/(2*a);
        printf("\n-----------------------------");
        printf("\nthe Roots of the eqn are distinct and equal to: %f and %f",x,y);
        printf("\n-----------------------------");
    }
    else
    {
        real_part = (-b)/(2*a);
        imag_part = sqrt(-d)/(2*a);
        printf("\n-----------------------------");
        printf("\n Roots are %f + %f i ",real_part,imag_part);
        printf("\n and %f - %f i ",real_part,imag_part);
        printf("\n-----------------------------");
    }
}
