
int Mean(){
	int n,i;
	double sum=0,mean;
	printf("\n Enter the value of N:");
	scanf("%d",&n);
	double arr[n];
	if (n<=0)
	printf("\n Invalid Input , take N > 0");
	printf("\n Enter %d values",n);
	for (i=0;i<n;i++){
		scanf("%lf",&arr[i]);
	}
	for (i=0;i<n;i++){
				sum = sum+arr[i];
		}
	mean = sum/n;
	printf("\n The Mean of given Observations is : %lf",mean);
}
int Variance(){
	int N,i;
	double sum=0,mean,var = 0;
	printf("\n Enter the value of N:");
	scanf("%d",&N);
	double arr[N];
	if (N<=0)
	printf("\n Invalid Input , take N > 0");
	printf("\n Enter %d values",N);
	for (i=0;i<N;i++){
		scanf("%lf",&arr[i]);

	}
		for (i=0;i<N;i++){
		sum = sum+arr[i];
	}
	mean = sum/N;
	for (i=0;i<N;i++){
		var = var + pow(arr[i]-mean,2);
	}
	var = var/N;
	printf("\n The Variance of given Observations is : %lf",var);
	return var;
}
int Std_deviation(){
	double sd,v;
	v = Variance();
	sd = sqrt(v);
	printf("\n The Standard Deviation of given Observations is : %lf",sd);	
}
int CV(){
	int n,i;
	double sum=0,mean,sd,cv;;
	printf("\n Enter the value of N:");
	scanf("%d",&n);
	double arr[n];
	if (n<=0)
	printf("\n Invalid Input , take N > 0");
	printf("\n Enter %d values",n);
	for (i=0;i<n;i++){
		scanf("%lf",&arr[i]);
	}
	for (i=0;i<n;i++){
				sum = sum+arr[i];
		}
	mean = sum/n;
	printf("\n The Mean of given Observations is : %lf",mean);
	sd = Std_deviation();
	cv = sd/mean;
	printf("\n The Coefficient of Variation is :%lf",cv);
	getch();
}
int Range(){
	int n,i,j,k;
	double rng;
	printf("\n Enter the size of the array :");
	scanf("%d",&n);
	double array[n];
	printf("\n enter observation :");
	for (i = 0;i<n;i++){
		scanf("%lf",&array[i]);	
	}
	j = array[0];
	for (i =1;i<n;i++){
	
		if (array[i]>j)
		j = array[i];
}
	printf("\n max observation is %lf",j);
	k = array[0];
	for (i =1;i<n;i++)
		if (array[i]<k)
	k = array[i];
	printf("\n min observation is %lf",k);
	rng = j-k;
	printf("\n----------------------------");
	printf("\n Range = %lf\n",rng);
	printf("-----------------------------");
}


int Skewness_Kurtosis(){
	int n,i;
	double mean,sum=0,var=0,dev = 0,devsq = 0,mu3,mu2,mu4,beta1,beta2,gamma1,gamma2;
	printf("\n Enter the value of n:");
	scanf("%d",&n);
	double arr[n];
	printf("\n Enter the observations one by one:");
	for(i=0;i<n;i++){
		scanf("%lf",&arr[i]);
	}
	for (i=0;i<n;i++){
		sum = sum + arr[i];
	}
	mean = sum/n;
	for(i=0;i<n;i++){
		var = var + pow(arr[i]-mean,2);
	}
	mu2 = var/n;
	for(i=0;i<n;i++){
		dev = dev + pow(arr[i]-mean,3);
	}
	mu3 = dev/n;
		for(i=0;i<n;i++){
		devsq = devsq + pow(arr[i]-mean,4);
	}
	mu4 = devsq/n;
        beta1 = ((mu3 * mu3) / (mu2 * mu2 * mu2));
        gamma1 = sqrt(beta1);
        printf("\n---------------------------------------------");
        printf("\nSkewness = %lf \n or Skewness Coefficient = %lf\n", beta1, gamma1);
     	printf("\n----------------------------------------------");
        if (gamma1 > 0)
    		printf("\n The distribution is positively skewed (tail to the right).");
		else if (gamma1 < 0)
		    printf("\n The distribution is negatively skewed (tail to the left).");
		else
    		printf("\n The distribution is symmetric.");
    printf("\n----------------------------------------------");
	beta2 = (mu4 / (mu2 * mu2));
	gamma2 = beta2 - 3;

	printf("\nKurtosis = %lf\nor Kurtosis Coefficient = %lf\n", beta2, gamma2);
	printf("\n----------------------------------------------");
	if (gamma2 > 3)
    printf("\n? The distribution is leptokurtic (peaked, heavy tails).");
	else if (gamma2 < 3)
    printf("\n The distribution is platykurtic (flat, light tails).");
	else
    printf("\n The distribution is mesokurtic (normal shape).");
	printf("\n----------------------------------------------");
}
int GeometricMean()
{
    int n,i;
	double prod=0,gm;
	printf("\n Enter the value of N:");
	scanf("%d",&n);
	double arr[n];
	if (n<=0)
	printf("\n Invalid Input , take N > 0");
	printf("\n Enter %d values",n);
	for (i=0;i<n;i++){
		scanf("%lf",&arr[i]);
	}	 
    for (i=0;i<n;i++){
        prod = prod * arr[i];
    }
    n = i;
    printf("\nGeometric Mean = %Lf\n", gm);
}
int HarmonicMean()
{
 	int n,i;
	double sum=0,hm;
	printf("\n Enter the value of N:");
	scanf("%d",&n);
	double arr[n];
	if (n<=0)
	printf("\n Invalid Input , take N > 0");
	printf("\n Enter %d values",n);
	for (i=0;i<n;i++){
		scanf("%lf",&arr[i]);
	}	 
    for (i=0;i<n;i++){
        sum = sum + (1 / arr[i]);
    }
    n = i;
    hm = n / sum;
    printf("\nHarmonic Mean = %.3Lf\n", hm);
}


double findMode(double arr[], int size) {
    int i, j;
    int maxCount = 0;
    double maxValue = arr[0];

    for (i = 0; i < size; ++i) {
        int count = 0;

        for (j = 0; j < size; ++j) {
            if (arr[j] == arr[i])
                ++count;
        }

        if (count > maxCount) {
            maxCount = count;
            maxValue = arr[i];
        }
    }

    return maxValue;
}

int Mode() {
    int n, i;
    printf("\nEnter the size of the array: ");
    scanf("%d", &n);

    double data[n];
    printf("\nEnter observations:\n");
    for (i = 0; i < n; i++) {
        scanf("%lf", &data[i]);
    }

    double mode = findMode(data, n);
    printf("Mode = %.2lf\n", mode);

    return 0;
}

