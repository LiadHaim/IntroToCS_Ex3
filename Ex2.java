package exe.ex2;

import java.util.Arrays;

/**
 * Introduction to Computer Science 2023, Ariel University,
 * Ex2: arrays, static functions and JUnit
 *
 * This class represents a set of functions on a polynom - represented as array of doubles.
 * The array {0.1, 0, -3, 0.2} represents the following polynom: 0.2x^3-3x^2+0.1
 * This is the main Class you should implement (see "add your code here")
 *
 * @author boaz.benmoshe
 */
public class Ex2 {
	/** Epsilon value for numerical computation, it serves as a "close enough" threshold. */
	public static final double EPS = 0.001; // the epsilon to be used for the root approximation.
	/** The zero polynom is represented as an array with a single (0) entry. */
	public static final double[] ZERO = {0};
	/**
	 * Computes the f(x) value of the polynom at x.
	 * @param poly
	 * @param x
	 * @return f(x) - the polynom value at x.
	 */
	public static double f(double[] poly, double x) {
		double ans = 0;
		for(int i=0;i<poly.length;i++) {
			double c = Math.pow(x, i);
			ans +=c*poly[i];
		}
		return ans;
	}
	/** Given a polynom (p), a range [x1,x2] and an epsilon eps. 
	 * This function computes an x value (x1<=x<=x2) for which |p(x)| < eps, 
	 * assuming p(x1)*p(x1) <= 0. 
	 * This function should be implemented recursively.
	 * @param p - the polynom
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
	 * @return an x value (x1<=x<=x2) for which |p(x)| < eps.
	 */
	public static double root_rec(double[] p, double x1, double x2, double eps) {
		double f1 = f(p,x1);
		double f2 = f(p,x2);
		double x12 = (x1+x2)/2;
		double f12 = f(p,x12);
		if (f1*f2<=0 && Math.abs(f12)<eps) {return x12;}
		if(f12*f1<=0) {return root_rec(p, x1, x12, eps);}
		else {return root_rec(p, x12, x2, eps);}
	}

	/**
	 * This function computes a polynomial representation from a set of 2D points on the polynom.
	 * The solution is based on: //	http://stackoverflow.com/questions/717762/how-to-calculate-the-vertex-of-a-parabola-given-three-points
	 * Note: this function only works for a set of points containing up to 3 points, else returns null.
	 * @param xx
	 * @param yy
	 * @return an array of doubles representing the coefficients of the polynom.
	 */
	public static double[] PolynomFromPoints(double[] xx, double[] yy) {
		double [] ans = null;
		int lx = xx.length;
		int ly = yy.length;
		if(xx!=null && yy!=null && lx==ly && lx>1 && lx<4) {
			double x1 = xx[0], x2 = xx[1], x3 = xx[2];
			double y1 = yy[0], y2 = yy[1], y3 = yy[2];
			double a = (x1 * (y3 - y2) + x2 * (y1 - y3) + x3 * (y2 - y1)) / ((x1 - x2) * (x1 - x3) * (x2 - x3));
			double b = (y2 - y1) / (x2 - x1) - (a * (x1 + x2));
			double c = y1 - a * (x1 * x1) - b * x1;
			double[] ans1 = {c, b, a};
			ans = ans1;
		}
		return ans;
	}

	/** Two polynoms are equal if and only if the have the same values f(x) for 1+n values of x, 
	 * where n is the max degree (over p1, p2) - up to an epsilon (aka EPS) value.
	 * @param p1 first polynom
	 * @param p2 second polynom
	 * @return true iff p1 represents the same polynom as p2.
	 */

	public static boolean equals(double[] p1, double[] p2) {
		//checking if p1 or p2 has an 0 for a constant and if so, remove it
		int i = p1.length - 1; 
		while (i >= 0 && p1[i] == 0) {
			i--;
		}
		int j = p2.length - 1;
		while (j >= 0 && p2[j] == 0) {
			j--;
		}

		if (i != j) {		//if their sizes are different, return false
			return false;
		}
		for (int k = 0; k <= i; k++) {
			if (Math.abs(f(p1, k) - f(p2, k)) > EPS) {		//comparing f(x) values in range [0,max power of p1/p2]
				return false;
			}
		}
		return true;
	}

	private static double f(int[] p, int x) {
		double result = 0;
		for (int i = p.length - 1; i >= 0; i--) {
			result = result * x + p[i];
		}
		return result;
	}


	/** 
	 * Computes a String representing the polynom.
	 * For example the array {2,0,3.1,-1.2} will be presented as the following String  "-1.2x^3 +3.1x^2 +2.0"
	 * @param poly the polynom represented as an array of doubles
	 * @return String representing the polynom: 
	 */
	public static String poly(double[] poly) {
		// get poly 
		// return the string of poly
		// add it to string and return val of string
		if(poly.length==0)
		{
			return "";
			//return if length is 0  
		}
		String ans = "";
		int counter=0; //init counter for assigning the array size of the polynomial
		for(int i=poly.length-1;i>=0;i--)
		{

			if(poly[i]==0)
			{
				counter++;
				continue ;
			}

			if(i!=0)
			{
				ans=ans+poly[i]+"x^"+i;
				if(poly[i-1]>0)
				{
					ans=ans+"+";
				}	

			}
			if(i==0&&poly[0]>0)
			{
				ans=ans+poly[i];
			}
			if(i==0&&poly[0]<0)
			{
				ans=ans+poly[i];
			}		
		}
		if(counter==poly.length)
		{
			return "0";
		}

		return ans;
	}
	/**
	 * Given two polynoms (p1,p2), a range [x1,x2] and an epsilon eps. This function computes an x value (x1<=x<=x2)
	 * for which |p1(x) -p2(x)| < eps, assuming (p1(x1)-p2(x1)) * (p1(x2)-p2(x2)) <= 0.
	 * @param p1 - first polynom
	 * @param p2 - second polynom
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
	 * @return an x value (x1<=x<=x2) for which |p1(x) - p2(x)| < eps.
	 */
	public static double sameValue(double[] p1, double[] p2, double x1, double x2, double eps) {
		double x12 = (x1+x2)/2;
		// *** add your code here ***
		while (Math.abs(f(p1, x12) - f(p2, x12)) > eps) {// the abs value of the function (p1,x12) and (p2,x12) is bigger then EPS
			if ((f(p1, x1) < f(p2, x1) && f(p1, x12) < f(p2, x12)) // check when the value of x is the same in both functions when the value of y is bigger
					|| (f(p1, x1) > f(p2, x1) && f(p1, x12) > f(p2, x12))) {
				x1 = x12;//x1 val take the val of x12
			} else {
				x2 = x12;//x2 val take the val of x12
			}
			x12 = (x1 + x2) / 2;//x12 equals to average of x1 + x2
		}
		// **************************
		return x12;
	}
	/**
	 * Given a polynom (p), a range [x1,x2] and an integer with the number (n) of sample points. 
	 * This function computes an approximation of the length of the function between f(x1) and f(x2) 
	 * using n inner sample points and computing the segment-path between them.
	 * assuming x1 < x2. 
	 * This function should be implemented iteratively (none recursive).
	 * @param p - the polynom
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param numberOfSegments - (A positive integer value (1,2,...).
	 * @return the length approximation of the function between f(x1) and f(x2).
	 */
	public static double length(double[] p, double x1, double x2, int numberOfSegments) {
		double ans = 0;
		double segLen = (x2 - x1) / numberOfSegments; //calculate the width of the segment
		for (int i = 0; i < numberOfSegments; i++) {
			double x1Seg = (x1 + i )* segLen;
			double x2Seg = x1 + (i + 1) * segLen;
			double y1 = f(p, x1Seg);
			double y2 = f(p, x2Seg);
			double dx = x2Seg - x1Seg;
			double dy = y2 - y1;
			double segLenApprox = Math.sqrt(dx * dx + dy * dy); //calculating each segment's polynomuial's length using pythagoras
			ans += segLenApprox;
		}
		return ans;
	}



	/**
	 * Given two polynoms (p1,p2), a range [x1,x2] and an integer representing the number of Trapezoids between the functions (number of samples in on each polynom). 
	 * This function computes an approximation of the area between the polynoms within the x-range.
	 * The area is computed using Riemann's like integral (https://en.wikipedia.org/wiki/Riemann_integral)
	 * @param p1 - first polynom
	 * @param p2 - second polynom
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param numberOfTrapezoid - a natural number representing the number of Trapezoids between x1 and x2.
	 * @return the approximated area between the two polynoms within the [x1,x2] range.
	 */

	public static double area(double[] p1, double[] p2, double x1, double x2, int numberOfBoxes) {
		double ans = 0;
		double x0=x1;
		double y01 = Ex2.f(p1, x0);
		double y02 = Ex2.f(p2, x0);
		double delta = (x2-x1)/numberOfBoxes;	//calculate the width of the segment
		for(double x = x1+delta;x<=x2;x+=delta) {
			double y11 = Ex2.f(p1, x);
			double y12 = Ex2.f(p2, x);
			double y = Math.abs(y01-y02);
			double y1 = Math.abs(y11-y12);
			ans+=((y+y1)*delta)*0.5; 	//calculating each Trapezoids's area and adding each one to their sum
			x0=x;
			y01 = y11;
			y02=y12;
		} 
		return ans;
	}
	/**
	 * This function computes the array representation of a polynom from a String
	 * representation. Note:given a polynom represented as a double array,  
	 * getPolynomFromString(poly(p)) should return an array equals to p.
	 * 
	 * @param p - a String representing polynom.
	 * @return
	 */
		
		public static double[] getPolynomFromString(String p) {
			String h = p.replace(" ", "");
			double[]arr11= {};
			if(h.length()==0)
			{
				return arr11;
			}
			String s=h;
			s.replace("X", "x");// replace all occurrences of "X" with "x" in input string
			//	for (int i=0;i<s.length();i++)
			//	{
			//		System.out.println(s.charAt(i));
			//	}
			int index=0;
			int counter=0;
			int Size=ArrSize(s);
			double[]arr=new double[Size+1];
			String[]temp=s.split("((?=-)|(?=-)|(?>\\+)|(?=\\+))");// split input string into an array of monomials
			//System.out.println(Arrays.toString(temp));
			String[]t = new String[temp.length];
			for(int i=0;i<temp.length;i++)
			{
				t=temp[i].split("\\^");// split each monomial into its coefficient and exponent
				for(int j=0;j<t.length;j++)
				{
					if(t.length==1)// if the monomial has no exponent
					{
						if(Checker(t[j])==1)// if the monomial has an "x" term but no exponent
						{
							String []y=t[j].split("x");// split the coefficient and "x" term
							double x=Double.parseDouble(y[0]); // parse the coefficient to a double
							System.out.println(x);
							arr[1]=x;// set the coefficient of the x term in the array to the parsed coefficient
						}
						if(Checker(t[j])==0)// if the monomial has no "x" term or exponent
						{
							double x=Double.parseDouble(t[j]);// parse the coefficient to a double
							System.out.println(x);// set the constant term in the array to the parsed coefficient
							arr[0]=x;
						}

					}
					if(t.length==2)// if the monomial has an exponent
					{
						int  index1=Integer.parseInt(t[1]);// parse the exponent to an integer
						String []y=t[0].split("x"); // split the coefficient and "x" term
						if(y.length==0)//if the coefficient is missing
						{
							double x=1; // assume coefficient is 1
							arr[index1]=1; // set the coefficient of the x term in the array to 1
						}
						else
						{

							double x=Double.parseDouble(y[0]);// parse the coefficient to a double
							arr[index1]=x; // set the coefficient of the x term in the array to the parsed coefficient
						}

					}
				}
				System.out.print(Arrays.toString(t));
			}
			System.out.print(Arrays.toString(arr));
			return arr;
		}
		public static int Checker(String s)
		{
			int counter=0;
			for(int i=0;i<s.length();i++)
			{
				if(s.charAt(i)=='X'|s.charAt(i)=='x')
				{
					counter++;
				}
			}
			if(counter==0)
			{
				return 0;
			}

			return counter;

		}

	public static int ArrSize(String s) // method to get the size of the new array



	{
		String poly=s;
		String help=poly;
		int counter=0;
		int  tempo=0;
		int max=0;
		String[]temp=poly.split("((?=-)|(?=-)|(?>\\+)|(?=\\+))");
		for(int i=0;i<temp.length;i++)
		{
			String []mono=temp[i].split("x");
			for(int j=0;j<mono.length;j++)
			{
				if(mono.length>1)
				{
					if(mono[j].contains("^"))
					{
						String x=mono[j].substring(1);
						x.replace(" ", "");
						tempo=Integer.parseInt(x);

						if(tempo>max)
						{
							max=tempo;
						}

					}
				}
			}
		}
		if(max!=0)
			return max;
		{
		}
		for(int i=0;i<help.length();i++)
		{
			if(poly.charAt(i)=='X'|poly.charAt(i)=='x')
			{
				counter++;
			}
		}
		if(counter==0)
		{
			return max;
		}
		if(counter==1)
		{
			for(int i=0;i<help.length();i++)
			{
				if(poly.charAt(i)=='+'|poly.charAt(i)=='-')
				{
					counter++;
				}
			}

		}
		return counter;

	}

	/**
	 * This function computes the polynom which is the sum of two polynoms (p1,p2)
	 * @param p1
	 * @param p2
	 * @return
	 * 
	 */

	public static double[] add(double[] p1, double[] p2) {
		// add the p1 and p2 
		// do it by the value
		// take the value of the p1 and p2
		// enter p1[i]+p2[i] to the index in the new array l,
		double[]arr11= {};
		if(p2.length==0)// If one of the polynomials is empty, return the other polynomi
		{
			return p1;
		}
		if(p1.length==0)
		{
			return p2;
		}
		if(p2.length==0&&p1.length==0)// If both polynomials are empty, return an empty array
		{
			return arr11;
		}
		int min=Math.min(p1.length, p2.length);// Determine the minimum and maximum lengths of p1 and p2
		int max=Math.max(p1.length, p2.length);
		double[] arr=new double[max];// Initialize an array of length max to store the sum of p1 and p2
		if(p1.length<=p2.length)// Add the coefficients of p1 and p2 for the first min elements
		{
			for(int i=0;i<min;i++)
			{
				arr[i]=p2[i]+p1[i];
			}
			for(int i=min;i<max;i++)// Copy the remaining coefficients of p2 to the result array
			{
				arr[i]=p2[i];
			}
		}
		else
		{
			for(int i=0;i<min;i++)
			{
				arr[i]=p2[i]+p1[i];
			}
			for(int i=min;i<max;i++)// Copy the remaining coefficients of p1 to the result array
			{
				arr[i]=p1[i];
			}
		}
		return arr;
	}
	/**
	 * This function computes the polynom which is the multiplication of two polynoms (p1,p2)
	 * @param p1
	 * @param p2
	 * @return
	 */

	public static double[] mul(double[] p1, double[] p2) {
		double[]arr11= {};// creates an empty array to be returned in case both input arrays are empty
		if(p2.length==0)// checks if p2 is an empty array
		{
			return p1;// returns p1 if p2 is empty
		}
		if(p1.length==0) // checks if p1 is an empty array
		{
			return p2; // returns p2 if p1 is empty
		}
		if(p2.length==0&&p1.length==0)// checks if both p1 and p2 are empty arrays
		{
			return arr11;// returns an empty array if both p1 and p2 are empty
		}

		double[]arr=new double[p1.length+p2.length-1]; // creates an array of size (p1.length+p2.length-1) to store the result
		for(int i=0;i<arr.length;i++)
		{
			arr[i]=0.0;// initializes all elements of arr to 0.0
		}

		for(int i=0;i<p1.length;i++)// loops through all elements of p1
		{
			for(int j=0;j<p2.length;j++) // loops through all elements of p2
			{
				double temp=p1[i]*p2[j]; // multiplies the current element of p1 with the current element of p2
				arr[i+j]+=temp;// adds the result of the multiplication to the appropriate index of arr
			}
		}
		return arr;
	}
	/**
	 * This function computes the derivative polynom:.
	 * @param po
	 * @return
	 */
	public static double[] derivative (double[] po) {
		double[]arr11= {}; // Initialize empty array
		double[]arr12= {0};// Initialize array with single zero element

		if(po.length==0)// If the input array is empty, return the empty array
		{
			return arr11;
		}
		if(po.length==1)// If the input array has only one element, return an array with a single zero element
		{
			return arr12;
		}
		double[]arr=new double[po.length-1];// Initialize a new array with length equal to the input array minus one
		for(int i=1;i<=arr.length;i++)		//calculating derivative 
		{
			arr[i-1]=po[i]*i;// The derivative of the term with exponent i is i times the coefficient of that term
		}

		return arr; // Return the derivative array
	}
}