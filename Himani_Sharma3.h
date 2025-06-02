
#define SIZE 10
#define MX 10
// Define real functions instead of macros for better type safety
double f1(double x) { return cos(x) - x * exp(x); }
double df1(double x) { return -sin(x) - exp(x) * (x + 1); }
double ddf1(double x) { return -cos(x) - exp(x) * (x + 2); }

double f2(double x) { return x * x - 3; }
double df2(double x) { return 2 * x; }
double ddf2(double x) { return 2; }

double f3(double x) { return x * x - 5 * x + 1; }
double df3(double x) { return 2 * x - 5; }
double ddf3(double x) { return 2; }

double f4(double x) { return 1/(1+x) ; }


double (*f)(double);
double (*df)(double);
double (*ddf)(double);
int select_fn(double (**f)(double), double (**df)(double), double (**ddf)(double));
int Bisection_Method();
int Secant_Method();
int Newton_Raphson_Method();
int Bairstow_Method();
int Regula_Falsi();
int Simpson3by8();
int Simpson1by3();
int Trapezoidal();
int Jacobi_Iteration();
int Birge_Vieta();
int Lagrange_interpolation();
int Hermite_IP();
int Backward_Subsitution();
int Forward_Subtitution();
// Updated select_fn to use function pointers
int select_fn(double (**f)(double), double (**df)(double), double (**ddf)(double)) {
    int k;
    
    printf("\n 1) cos(x)-x*exp(x)");
    printf("\n 2) x^2 - 3");
	printf("\n 3) x^2 - 5*x + 1");
	printf("\n 4) 1/(1+x)");
	printf("\nChoose a function (1, 2, 3 or 4): ");
    scanf("%d", &k);

    switch (k) {
        case 1:
        	printf("\n 1) cos(x)-x*exp(x)\n");
            *f = f1;
            *df = df1;
            *ddf = ddf1;
            break;
        case 2:
        	printf("\n 2) x^2 - 3\n");
            *f = f2;
            *df = df2;
            *ddf = ddf2;
            break;
        case 3:
        	printf("\n 3) x^2 - 5*x + 1\n");
            *f = f3;
            *df = df3;
            *ddf = ddf3;
            break;
        case 4:
        	printf("\n 4) 1/(1+x) \n");
        	*f  = f4;
        	break;
        default:
            printf("\nInvalid input. Try again.\n");
            return 0;
    }
    return 1;
}

// Sample implementation for Bisection Method using the selected function
int Bisection_Method() {
    //double (*f)(double), (*df)(double), (*ddf)(double);
    if (!select_fn(&f, &df, &ddf)) return 1;

    double a, b, m, e;
    int k = 0;

    while (1) {
        printf("\nEnter interval [a, b]: ");
        scanf("%lf %lf", &a, &b);
        if ((f(a)) * (f(b)) < 0.0) break;
        printf("Invalid interval. Try again.\n");
    }

    printf("Enter error tolerance: ");
    scanf("%lf", &e);

    while (k < 6) {
        m = (a + b) / 2.0;
        printf("Iteration %d: a = %lf, b = %lf, m = %lf, f(m) = %lf\n", k + 1, a, b, m, f(m));
        if (fabs(f(m)) < e) {
            printf("Root found: %lf\n", m);
            break;
        }
        if ((f(a)) * (f(m)) < 0.0) {
            b = m;
        } else {
            a = m;
        }
        k++;
    }
    return 0;
}
int NR(double x0) {
	if (!select_fn(&f, &df, &ddf)) return 1;
    double x1,e;
    int iteration = 1;
	printf("\n Enter Error tolerance:");
	scanf("%lf",&e);
    while (1) {
        double fx = f(x0);
        double dfx = df(x0);

        if (dfx == 0.0) {
            printf("Math Error: Derivative is zero. No solution.\n");
            return;
        }

        x1 = x0 - fx / dfx;

        printf("Iteration %d: x = %.6lf\n", iteration, x1);

        if (fabs(x1 - x0) < e) {
            printf("\nRoot found: %.6lf\n", x1);
            break;
        }

        x0 = x1;
        iteration++;
    }
}
int Newton_Raphson_Method() {
    double x0;
    printf("Enter initial guess: ");
    scanf("%lf", &x0);
    NR(x0);
}
int Bairstow_Method()
{
     float a[10], b[10], c[10], p, q, cc, den, delp, delq;
     int i, n, m, j, k, l;
     printf("Input p, q, degree, number of iterations:\n");
     scanf("%f%f%d%d",&p,&q,&n,&m);
     printf("Input coefficients of polynomial in decreasing order of n \n");
     for(i = 0; i <= n; i++)
      {
       scanf("%f", &a[i]);
      }
//  generate b[k] & c[k]
     for (j = 0;j <= m;j++)
      {
	b[0] = a[0];
	b[1] = a[1] - p * b[0];
	b[2] = a[2] - p * b[1] - q;
	for (k = 3; k <= n; k++)
	  b[k] = a[k] - p * b[k-1] - q * b[k-2];
	printf("\n b[0] = %f b[1] = %f b[2] = %f b[%d] = %f \n",b[0],b[1],b[2],k,b[--k]);

	c[0] = b[0];
	c[1] = b[1] - p*c[0];
	c[2] = b[2] - p * c[1] - q;
	l = n - 1;
	for (k = 3; k <= l; k++)
	  c[k] = b[k] - p * c[k-1] - q * c[k-2];
	cc = c[n-1] - b[n-1];
	den = c[n-2] * c[n-2] - cc * c[n-3];
	delp = -(b[n] * c[n-3] - b[n-1] * c[n-2]) / den;
	delq = -(b[n-1] * cc - b[n] * c[n-2]) / den;
	printf("\n den = %f \n delp = %f \n delq = %f \n",den,delp,delq);
	p = p + delp;
	q = q + delq;
	printf("ITERATIONS = %d, P = %f, Q = %f \n", j, p, q);
	}
}
int check(double x1,double x2)
	{
//	if (!select_fn(&f, &df, &ddf)) return 1;
		if ((f(x1))*(f(x2)) < 0.0)
			return(1);
		else{
			return(0);
		}
	}
double secant(double x1,double x2,double tolerance){
	double x3;
	int k = 1;
//	if (!select_fn(&f, &df, &ddf)) return 1;
	printf("\n Iteration\t  x1\t\t  x2\t\t  x3\t\t f(x3)\n");
    printf("-----------------------------------------------------------------------\n");
    
	while(k<=6){
		x3 = (x1*(f(x2))-x2*(f(x1)))/((f(x2))-(f(x1))); 
		printf("k = %d\t\t%.6lf\t%.6lf\t%.6lf\t%.6lf\n", k, x1, x2, x3,(f(x3)));
		if ((f(x3))==0.0)
		{
			return(x3);
		}
		else{
			x1 = x2;
			x2 = x3;
		}
			
	k++;
    }
	return (x3);
}
int Secant_Method(){
	if (!select_fn(&f, &df, &ddf)) return 1;
	double x1,x2,e,root;
	printf("\n Enter Initial Approximations :");
	scanf("%lf %lf",&x1,&x2);
	if (check(x1,x2)) {
        printf("\nThe interval [%.6lf, %.6lf] brackets a root.\n", x1,x2);
    } else {
        printf("\nThe interval [%.6lf, %.6lf] does NOT bracket a root.\n", x1,x2);
    }
    
    printf("\n enter error tolerance:");
    scanf("%lf",&e);
    
    root = secant(x1,x2,e);
    printf("\n Approximated Root is %lf",root);
}
double RegulaFalsi(double x1,double x2,double tolerance){
	double x3;
	int k = 1;
	
	printf("\n Iteration\t  x1\t\t  x2\t\t  x3\t\t f(x3)\n");
    printf("-----------------------------------------------------------------------\n");
    
	while(k<=6){
		x3 = (x1*(f(x2))-x2*(f(x1)))/((f(x2))-(f(x1))); 
		printf("k = %d\t\t%.6lf\t%.6lf\t%.6lf\t%.6lf\n", k, x1, x2, x3,(f(x3)));
		if ((f(x3))==0.0)
		{
			return(x3);
		}
		else{
			if ((f(x1)*(f(x3)))<0.0){
					x1 = x1;
					x2 = x3;
			}
			else {
				x1 = x3;
				x2 = x2;
			}
		}
			
	k++;
    }
	return (x3);
}
int Regula_Falsi(){
	if (!select_fn(&f, &df, &ddf)) return 1;
	double x1,x2,e,root;
	printf("\n Enter Initial Approximations :");
	scanf("%lf %lf",&x1,&x2);
	if (check(x1,x2)) {
        printf("\nThe interval [%.6lf, %.6lf] brackets a root.\n", x1,x2);
    } else {
        printf("\nThe interval [%.6lf, %.6lf] does NOT bracket a root.\n", x1,x2);
    }
    printf("\n enter error tolerance:");
    scanf("%lf",&e);
    
    root = RegulaFalsi(x1,x2,e);
    printf("\n Approximated Root is %lf",root);
}
double calcx4(double x1,double x2,double x3)
{
    double f1,f2,f3,h1,h2,lambda2,s2,g2,a1,c2,lambda3,x4;
//    printf("\nf1=%lf f2=%lf f3=%lf",x1,x2,x3);
    f1=f(x1);
    f2=f(x2);
    f3=f(x3);
    h1=x2-x1;
    h2=x3-x2;
    lambda2=h2/h1;
    s2=1+lambda2;
    g2=(lambda2*lambda2*f1)-(s2*s2*f2)+((lambda2+s2)*f3);
    c2=lambda2*(lambda2*f1-s2*f2+f3);
    if(g2<0.00)
      lambda3=(-2*s2*f3)/(g2-sqrt(g2*g2-4*s2*f3*c2));
    else if(g2>0.00)
      lambda3=(-2*s2*f3)/(g2+sqrt(g2*g2-4*s2*f3*c2));
    else
      printf("ERROR:g2 cannot be negative and it is giving negative ");
    x4=x3+(x3-x2)*lambda3;
//printf("\nh1=%lf   h2=%lf  \nlemda2 =%lf   s2=%lf\n  g2=%lf c2=%lf\n lambda3=%lf \n x4=%lf",h1,h2,lemda2,s2,g2,c2,lemda3,x4);
return x4;
}
int Muller_Method()
{
int k=1;
double a,d,b,x1,x2,x3,x4,f4,e;
//code for initial approximation
printf("\nenter the three initial approxiation ");
scanf("%lf%lf%lf",&x1,&x2,&x3);
printf("\nf(x1)=%lf",f(x1));
printf("\nf(x2)=%lf",f(x2));
printf("\nf(x3)=%lf",f(x3));
printf("\n x1,x2 and x3 are %lf,%lf and %lf",x1,x2,x3);
//code for error tolerance
printf("\n enter the error tolerance");
scanf("%lf",&e);
printf("\n -------------------Muller-Method-Table---------------------");
printf("\n|  k |   x(n-1)  |    xn   | x(n+1)   |  x(n+2) | f(x(n+2) |  x(n+2)-x(n+1)| ");
     while(k<=5)
	    {
	    x4=calcx4(x1,x2,x3);
	   // d=x4-x3;
	  //  printf("\nd=%lf",d);
	    printf("\n|  %d | %8.4lf  |%8.4lf | %8.4lf  |%8.4lf | %8.4lf | %8.4lf |",k,x1,x2,x3,x4,f(x4),d);
	    k++;
	    double c=fabs(f(x4));
		  if(c<e)
			{
			printf("\n -----------------------------------------------------------" );
			printf("\nRoot is %lf",x4);
			break;
			}
		  else
			{
			x1=x2;
			x2=x3;
			x3=x4;
			}
	    }
}
double Cheby(double x0,double tolerance){
	double x1;
	int k = 1;
	
	printf("\nIteration\t  x0\t\t  x1\t\t  f(x0)\n");
    printf("-----------------------------------------------------------------------\n");
    
	while(k<=6){
		double fx = f(x0);
        double dfx = df(x0);
        double ddfx = ddf(x0);
        
        //Iteration function
        x1 = x0 - (fx / dfx) - (1.0/2)*((fx*fx)/(dfx*dfx*dfx))*(ddfx);
		
		printf("%d\t\t%.6lf\t%.6lf\t%.6lf\n", k, x0, x1, fx);
		if (fabs(x1-x0)<tolerance)	{
			return(x1);
	    }
		x0 = x1;
		k++;
    }
	return (x1);
}
int Chebyshev(){
	if (!select_fn(&f, &df, &ddf)) return 1;
	double x0,e,root;
	printf("\n Enter Initial Approximation :");
	scanf("%lf",&x0);
	
    printf("\n enter error tolerance:");
    scanf("%lf",&e);
    
    root = Cheby(x0,e);
    printf("\n----------------------------------------------");
    printf("\n Approximated Root is %lf",root);
       printf("\n----------------------------------------------");
}
int Trapezoidal()
{
	if (!select_fn(&f, &df, &ddf)) return 1;
	float a,b,trapezoidal;
	//printf("f(x)=integration from a to b (1/(1+x))\n");
	printf("Enter the value of a and b:\n");
	scanf("%f%f",&a,&b);
	trapezoidal=((f(a)+f(b))/2.0)*(b-a);
	printf("\n----------------------------------------------");
	printf("\nValue of f(x) using trapezoidal rule is %f",trapezoidal);
   printf("\n----------------------------------------------");
}

int Simpson1by3()
{
	if (!select_fn(&f, &df, &ddf)) return 1;
	float a,b,simpson1_3;
	//printf("f(x)=integration from a to b (1/(1+x))\n");
	printf("Enter the value of a and b:\n");
	scanf("%f%f",&a,&b);
	simpson1_3=(f(a)+(4*f((a+b)/2))+f(b))*((b-a)/6.0);
	   printf("\n----------------------------------------------");
	printf("\nValue of f(x) using simpson one-third rule is %f",simpson1_3);
	   printf("\n----------------------------------------------");
}
int Simpson3by8()
{
	if (!select_fn(&f, &df, &ddf)) return 1;
	float a,b,simpson3_8;
	//printf("f(x)=integration from a to b (1/(1+x))\n");
	printf("Enter the value of a and b:\n");
	scanf("%f%f",&a,&b);
	simpson3_8=(f(a)+(3*f((2*a+b)/3.0))+(3*f((a+2*b)/3.0))+f(b))*((3*(b-a)/(8.0*3.0)));
	   printf("\n----------------------------------------------");
	printf("\nValue of f(x) using simpson three-eighth rule is %f",simpson3_8);
	   printf("\n----------------------------------------------");
}


int Jacobi_Iteration() {
    int n, i, j,k, iterations;
    float a[SIZE][SIZE], b[SIZE], x[SIZE], x_new[SIZE];
    printf("Enter the number of equations: ");
    scanf("%d", &n);
    printf("Enter the augmented matrix coefficients (a[i][j]) and constant terms (b[i]):\n");
    for (i = 0; i < n; i++) {
        printf("Equation %d:\n", i + 1);
        for (j = 0; j < n; j++) {
            printf("a[%d][%d] = ", i + 1, j + 1);
            scanf("%f", &a[i][j]);
        }
        printf("b[%d] = ", i + 1);
        scanf("%f", &b[i]);
    }

    printf("Enter initial guesses for the variables:\n");
    for (i = 0; i < n; i++) {
        printf("x[%d] = ", i + 1);
        scanf("%f", &x[i]);
    }

    printf("Enter number of iterations: ");
    scanf("%d", &iterations);
    // Jacobi Iteration
    for ( k = 1; k <= iterations; k++) {
        for (i = 0; i < n; i++) {
            float sum = 0;
            for (j = 0; j < n; j++) {
                if (j != i)
                    sum += a[i][j] * x[j];
            }
            x_new[i] = (b[i] - sum) / a[i][i];
        }
        // Update x values
        for (i = 0; i < n; i++) {
            x[i] = x_new[i];
        }
        // Display current iteration results
        printf("\nAfter iteration %d:\n", k);
        for (i = 0; i < n; i++) {
            printf("x[%d] = %.6f\n", i + 1, x[i]);
        }
    }
    return 0;
}


typedef struct {
    int degree;
    double coeff[10];  
} Polynomial;

// Function to evaluate polynomial using Birge-Vieta method
double birge_vieta(Polynomial poly, double x0) {
    int i, iter = 0;
    const double EPSILON = 0.0001;
    double b[10], c[10];
    double x = x0, fx, fpx;

    do {
        // Synthetic division to compute b coefficients
        b[poly.degree] = poly.coeff[poly.degree];
        for (i = poly.degree - 1; i >= 0; i--) {
            b[i] = poly.coeff[i] + x * b[i + 1];
        }

        fx = b[0];

        // Compute c coefficients for derivative f'(x)
        c[poly.degree] = b[poly.degree];
        for (i = poly.degree - 1; i > 0; i--) {
            c[i] = b[i] + x * c[i + 1];
        }

        fpx = c[1];
        if (fpx == 0) {
            printf("Derivative became zero. Method fails.\n");
            return x;
        }

        x0 = x;
        x = x - (fx / fpx);
        iter++;
    } while (fabs(x - x0) > EPSILON && iter < 100);

    return x;
}

int Birge_Vieta() {
    Polynomial poly;
    int i;
    double x0, root;

    printf("Enter the degree of the polynomial: ");
    scanf("%d", &poly.degree);

    if (poly.degree > 9) {
        printf("Degree too high (max supported: 9)\n");
        return 1;
    }

    printf("Enter the coefficients (from highest to lowest degree):\n");
    for (i = poly.degree; i >= 0; i--) {
        printf("Coefficient of x^%d: ", i);
        scanf("%lf", &poly.coeff[i]);
    }

    printf("Enter the initial guess: ");
    scanf("%lf", &x0);

    root = birge_vieta(poly, x0);
    printf("\nApproximate root: %.6f\n", root);

    return 0;
}
// Function to calculate the Lagrange Interpolation
double lagrange_interpolation(double x[], double y[], int n, double x_value) {
    double result = 0.0;
	int i,j;
    // Loop over each term in the summation
    for (i = 0; i < n; i++) {
        double term = y[i];  // Initialize the term with the y value
        
        // Compute the Lagrange basis polynomial L_i(x) for each term
        for (j = 0; j < n; j++) {
            if (j != i) {
                term *= (x_value - x[j]) / (x[i] - x[j]);
            }
        }
        
        // Add the current term to the result
        result += term;
    }
    
    return result;
}

int Lagrange_interpolation() {
    int n,i;
    printf("Enter the number of data points (n): ");
    scanf("%d", &n);

    double x[n], y[n];

    // Input the x and y values of the data points
    printf("Enter the x values:\n");
    for (i = 0; i < n; i++) {
        scanf("%lf", &x[i]);
    }

    printf("Enter the y values:\n");
    for (i = 0; i < n; i++) {
        scanf("%lf", &y[i]);
    }

    // Get the x_value where we want to interpolate
    double x_value;
    printf("Enter the value of x to interpolate: ");
    scanf("%lf", &x_value);

    // Calculate the interpolated value at x_value
    double result = lagrange_interpolation(x, y, n, x_value);

    // Output the result
    printf("The interpolated value at x = %.2lf is: %.2lf\n", x_value, result);

    return 0;
}
int computeDividedDifference(double x[], double y[][10], int n)
 {
 	int i,j;
    for (i = 1; i < 2 * n; i++)
     {
        for (j = 2 * n - 1; j >= i; j--)
         {
            if (x[j] == x[j - i])
             {
                y[j][i] = y[j][1]; // Using derivative values for repeated nodes
            } else
             {
                y[j][i] = (y[j][i - 1] - y[j - 1][i - 1]) / (x[j] - x[j - i]);
            }
        }
    }
}

// Function to interpolate using Hermite interpolation formula
double hermiteInterpolation(double x[], double y[][10], int n, double x_interp) 
{
    double result = y[0][0];
    double product = 1.0;
	int i;
    for (i = 1; i < 2 * n; i++)
     {
       product = product * (x_interp - x[i - 1]);
       result =  result + product * y[i][i];
    }

    return result;
}

int Hermite_IP() {
    int n,i;
    printf("Enter number of data points: ");
    scanf("%d", &n);

    double x[20], y[20][10];
    printf("Enter the values of x, y, and y':\n");
    
    for (i = 0; i < n; i++) 
    {
        printf("x[%d]: ", i);
        scanf("%lf", &x[2 * i]);
        x[2 * i + 1] = x[2 * i]; // Duplicate x values
        
        printf("y[%d]: ", i);
        scanf("%lf", &y[2 * i][0]);
        y[2 * i + 1][0] = y[2 * i][0]; // Duplicate y values
        
        printf("y'[%d]: ", i);
        scanf("%lf", &y[2 * i + 1][1]); // First derivative
    }

    computeDividedDifference(x, y, n);

    double x_interp;
    printf("Enter the value of x to interpolate: ");
    scanf("%lf", &x_interp);

    double interpolated_value = hermiteInterpolation(x, y, n, x_interp);
    printf("Interpolated value at x = %.5lf is %.5lf\n", x_interp, interpolated_value);

   getch();
}

int forwardSubstitution(float L[MX][MX], float b[MX], float x[MX], int n) 
{
	int i,j;
    for (i = 0; i < n; i++) 
    {
        x[i] = b[i];  
        for (j = 0; j < i; j++) 
        {
            x[i] -= L[i][j] * x[j];  
        }
        x[i] /= L[i][i];
    }
}

int displayEquations(float L[MX][MX], float b[MX], int n) 
{
	int i,j;
    printf("\nSystem of Equations:\n");
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++) 
        {
            printf("%.2fx%d ", L[i][j], j + 1);
            if (j < n - 1) 
            {
                printf("+ ");
            }
        }
        printf(" = %.2f\n", b[i]);
    }
}

void displayLowerTriangularMatrix(float L[MX][MX], int n) 
{
	int i,j;
    printf("\nLower Triangular Matrix (L):\n");
    for (i = 0; i < n; i++) 
    {
        for (j = 0; j < n; j++) 
        {
            printf("%.2f ", L[i][j]);
        }
        printf("\n");
    }
}

int Forward_Subtitution() 
{
    float L[MX][MX], b[MX], x[MX];
    int n;
	int i,j;
    printf("Enter the number of equations: ");
    scanf("%d", &n);

    printf("Enter the lower triangular matrix L (%d x %d):\n", n, n);
    for (i = 0; i < n; i++) 
    {
        for (j = 0; j < n; j++) 
        {
            scanf("%f", &L[i][j]);
        }
    }

    printf("Enter the right-hand side vector b:\n");
    for (i = 0; i < n; i++) {
        scanf("%f", &b[i]);
    }
    displayLowerTriangularMatrix(L, n);
    displayEquations(L, b, n);

    forwardSubstitution(L, b, x, n);

    printf("\nSolution vector x:\n");
    for (i = 0; i < n; i++) 
    {
        printf("x[%d] = %.4f\n", i + 1, x[i]);
    }

    return 0;
}

int Backward_Subsitution() 
{
    float L[MX][MX], b[MX], x[MX];
    int n,i,j;
    printf("Enter the number of equations: ");
    scanf("%d", &n);

    printf("Enter the lower triangular matrix L (%d x %d):\n", n, n);
    for (i = 0; i < n; i++) 
    {
        for (j = 0; j < n; j++) 
        {
            scanf("%f", &L[i][j]);
        }
    }

    printf("Enter the right-hand side vector b:\n");
    for (i = 0; i < n; i++) 
    {
        scanf("%f", &b[i]);
    }

    displayLowerTriangularMatrix(L, n);
    displayEquations(L, b, n);
    forwardSubstitution(L, b, x, n);

    printf("\nSolution vector x:\n");
    for (i = 0; i < n; i++) 
    {
        printf("x[%d] = %.4f\n", i + 1, x[i]);
    }

    return 0;
}

