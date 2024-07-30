package exe.ex2;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

/**
 *  * Introduction to Computer Science 2023, Ariel University,
 *  * Ex2: arrays, static functions and JUnit
 *
 * This JUnit class represents a JUnit (unit testing) for Ex2 - 
 * It contains few testing functions for the polynum functions as define in Ex2.
 * Note: you should add additional JUnit testing functions to this class.
 *
 * @author boaz.ben-moshe
 */

class Ex2Test {
	static final double[] P1 ={2,0,3, -1,0}, P2 = {0.1,0,1, 0.1,3};
	static double[] po1 = {2,2}, po2 = {-3, 0.61, 0.2};;
	static double[] po3 = {2,1,-0.7, -0.02,0.02};
	static double[] po4 = {-3, 0.61, 0.2};

	
	
	@Test
	public void testroot_rec() {
		try {
			double empty = Ex2.root_rec(po2,-1,1,Ex2.EPS);
			// code that may cause stack overflow
		} catch (StackOverflowError e) {
			// handle stack overflow exception
			System.err.println("Stack overflow detected: " + e.getMessage());
			// do something to free up stack space or increase stack size
		}
	}

	@Test
	/**
	 * Tests that f(x) == poly(x).
	 */
	public void testF() {
		double fx0 = Ex2.f(po1, 0);
		double fx1 = Ex2.f(po1, 1);
		double fx2 = Ex2.f(po1, 2);
		assertEquals(fx0, 2, Ex2.EPS);
		assertEquals(fx1, 4, Ex2.EPS);
		assertEquals(fx2, 6, Ex2.EPS);
	}
	@Test
	public void  TestFo()
	{	// create two polynomial arrays with different values

		double[] p2= {6.8,8,1,4};
		double[] p1= {6.8001,8.0002,1.001,4.1};
		boolean t=Ex2.equals(p1, p2);	// compare the two arrays using the Ex2.equals method and save the result in a boolean variable

		assertEquals(false,Ex2.equals(p1, p2));
	}
	@Test
	public void  testf()// Define the input polynomial as an array of coefficients
	// Call the function f with input polynomial x and value -3

	{
		double[] x= {12,0,-1};
		assertEquals(3,Ex2.f(x,-3));
	}
	@Test
	/**
	 * Tests that p1(x) + p2(x) == (p1+p2)(x)
	 */
	public void  testF2() {
		double x = Math.PI;
		double[] po12 = Ex2.add(po1, po2);
		double f1x = Ex2.f(po1, x);
		double f2x = Ex2.f(po2, x);
		double f12x = Ex2.f(po12, x);
		assertEquals(f1x + f2x, f12x, Ex2.EPS);
	}
	@Test
	/**
	 * Tests that p1+p2+ (-1*p2) == p1
	 */
	public void  test8Add() {
		double[] p12 = Ex2.add(po1, po2);
		double[] minus1 = {-1};
		double[] pp2 = Ex2.mul(po2, minus1);
		double[] p1 = Ex2.add(p12, pp2);
		assertTrue(Ex2.equals(p1, po1));
	}
	@Test
	/**

	This method tests the add method of the Ex2 class
	*/
	public void  testAdd() {
		double[] p12 = Ex2.add(po1, po2);
		double[] minus1 = {-1};
		double[] pp2 = Ex2.mul(po2, minus1);
		double[] p1 = Ex2.add(p12, pp2);
		assertEquals(Ex2.poly(po1), Ex2.poly(p1));
	}
	@Test
	/**
	 * Tests that p1+p2 == p2+p1
	 */
	public void  testAdd2() {// Add the two given polynomials

		double[] p12 = Ex2.add(po1, po2);// Create a one-term polynomial with coefficient -1

		double[] p21 = Ex2.add(po2, po1);// Ensure that the result polynomial p1 is equal to the expected polynomial (which is po1)

		assertTrue(Ex2.equals(p12, p21));
	}
	@Test
	/**
	 * Tests that p1+0 == p1
	 */
	public void  testAdd3() {
		double[] p1 = Ex2.add(po1, Ex2.ZERO);
		assertTrue(Ex2.equals(p1, po1));
	}
	
	@Test
	public void  TestAddNew()// This test is to check the correctness of the 'add' function in the Ex2 class

	{    // create arrays for testing

		double[] x= {5,6,7};
		double[] x1= {5,2};
		double[]x3= {9,8,8,9,1};	
		double[]x2=Ex2.add(x, x1);
		assertEquals(false,Ex2.equals(x3,x2)); // check if the result of adding x and x1 equals x3 using the Ex2.equals() method

	}
	@Test
	/**
	 * Tests that p1*0 == 0
	 */
	public void  testMul1() {
		double[] p1 = Ex2.mul(po1, Ex2.ZERO);
		assertTrue(Ex2.equals(p1, Ex2.ZERO));
	}
	@Test
	/**
	 * Tests that p1*p2 == p2*p1
	 */
	public void  testMul2() {
		double[] p12 = Ex2.mul(po1, po2);
		double[] p21 = Ex2.mul(po2, po1);
		assertTrue(Ex2.equals(p12, p21));
	}
	@Test
	public void  mul()
	{
		double[] p2= {10,1,5};// create an array of doubles for polynomial 2
		double[] p1= {2,-5,4};// create an array of doubles for polynomial 1
		double []p3={20,-48,45,-21,20};// create an array of doubles for the result of multiplication of polynomial 1 and 2
		assertEquals(true, Ex2.equals(p3,Ex2.mul(p1, p2)));;;
	}
	@Test
	public void  testSameValue() {
		double[] po3 = {2,1,1};// the array
		double test = Ex2.sameValue(po3,po1,-1,2,Ex2.EPS);// the function
		assertEquals(0,test,Ex2.EPS);// we expect to have this
	}
	@Test
	/**
	 * Tests that p1(x) * p2(x) = (p1*p2)(x),
	 */
	public void  testMulDoubleArrayDoubleArray() {
		double[] xx = {0,1,2,3,4.1,-15.2222};
		double[] p12 = Ex2.mul(po1, po2);
		for(int i = 0;i<xx.length;i=i+1) {
			double x = xx[i];
			double f1x = Ex2.f(po1, x);
			double f2x = Ex2.f(po2, x);
			double f12x = Ex2.f(p12, x);
			assertEquals(f12x, f1x*f2x, Ex2.EPS);
		}
	}


	@Test
	public void  Computes() {
		double [] f= {5,2,3,4,9};// define an array of double values
		double x=2;// define the value of x
		assertEquals(197, Ex2.f(f, x));// compute f(x) and compare it to the expected 
	}

	@Test
	public void  equals0()
	{
		double[] p12= {6.8,8,-1,4,0,0,0};// define an array of double values
		double[] p1= {};// define an empty array of double values
		boolean ans=Ex2.equals(p12, p1);// check if the two arrays are equal
		assertEquals(false,Ex2.equals(p12, p1));// compare the result to the expected 
	}
	@Test
	public void  equalsFalse()
	{
		double[] p12= {6.8,8,-1,4,0,0,0};// define an array of double values
		double[] p1= {6.8,8,-1,4.1}; // define another array of double values
		boolean ans=Ex2.equals(p12, p1);// check if the two arrays are equal
		assertEquals(false,Ex2.equals(p12, p1));// compare the result to the expected 
	}
	@Test
	public void  testDerivativeArrayDoubleArray() {
		double[] p = {1,2,3}; // 3X^2+2x+1
		double[] dp1 = {2,6}; // 6x+2
		double[] dp2 = Ex2.derivative(p); // compute the derivative of the given polynomial
		assertEquals(dp1[0], dp2[0],Ex2.EPS);// compare the first value of the computed derivative to the expected value
		assertEquals(dp1[1], dp2[1],Ex2.EPS);// compare the second value of the computed derivative to the expected value

		assertEquals(dp1.length, dp2.length);// compare the lengths of the two arrays
	}
		@Test
		void AreaIntegal()
		{
			double[] x= {0,3};// define an array of double values
			double[] x1= {0,1,0};// define another array of double values
			assertEquals(9, Ex2.area(x, x1, 0 ,3, 1000000),Ex2.EPS); // compute the area under the curve and compare it to the expected value
		}
	@Test
	/**
	 * Tests a simple derivative examples - till ZERO.
	 */
	public void  testD9erivativeArrayDoubleArray() {
		double[] p = {1,2,3}; // 3X^2+2x+1
		double[] pt = {2,6}; // 6x+2
		double[] dp1 = Ex2.derivative(p); // 2x + 6
		double[] dp2 = Ex2.derivative(dp1); // 2
		double[] dp3 = Ex2.derivative(dp2); // 0
		double[] dp4 = Ex2.derivative(dp3); // 0
		assertTrue(Ex2.equals(dp1, pt));
		assertTrue(Ex2.equals(Ex2.ZERO, dp3));
		assertTrue(Ex2.equals(dp4, dp3));
	}
	@Test
	/** 
	 * Tests the parsing of a polynom in a String like form.
	 */
	public void testFromString() {
		double[] p = {-1.1,2.3,3.1}; // 3.1X^2+ 2.3x -1.1
		String sp2 = "3.1x^2 +2.3x -1.1";
		String sp = Ex2.poly(p);
		System.out.println("sp:  "+sp);
		double[] p1 = Ex2.getPolynomFromString(sp);
		double[] p2 = Ex2.getPolynomFromString(sp2);
		boolean isSame1 = Ex2.equals(p1, p);
		boolean isSame2 = Ex2.equals(p2, p);
		if(!isSame1) {fail();}
		if(!isSame2) {fail();}
		assertEquals(sp, Ex2.poly(p1));
	}
	@Test
	/**
	 * Tests the equality of pairs of arrays.
	 */
	public void testEquals() {
		double[][] d1 = {{0}, {1}, {1,2,0,0}};
		double[][] d2 = {Ex2.ZERO, {1+Ex2.EPS/2}, {1,2}};
		double[][] xx = {{-2*Ex2.EPS}, {1+Ex2.EPS*1.2}, {1,2,Ex2.EPS*2}};
		for(int i=0;i<d1.length;i=i+1) {
			assertTrue(Ex2.equals(d1[i], d2[i]));
		}
		for(int i=0;i<d1.length;i=i+1) {
			assertFalse(Ex2.equals(d1[i], xx[i]));
		}
	}

	@Test
	/**
	 * Tests is the sameValue function is symmetric.
	 */
	public void testSameValue2() {
		double x1=-4, x2=0;
		double rs1 = Ex2.sameValue(po1,po2, x1, x2, Ex2.EPS);
		double rs2 = Ex2.sameValue(po2,po1, x1, x2, Ex2.EPS);
		assertEquals(rs1,rs2,Ex2.EPS);
	}
	@Test
	/**
	 * Test the area function - it should be symmetric.
	 */
	public void testArea() {
		double x1=0, x2=4;
		double a1 = Ex2.area(po1, po2, x1, x2, 10);
		double a2 = Ex2.area(po2, po1, x1, x2, 10);
		assertEquals(a1,a2,Ex2.EPS);
	}
	@Test
	public void testlength() {
		double x1=0, x2=4;// define range
		double a1 = Ex2.length(po1,  x1, x2, 10);// calculate the length of po1 in the range [0,4] with 10 segments
		double a2 = Ex2.length(po2,  x1, x2, 10); // calculate the length of po2 in the range [0,4] with 10 segments
		System.out.println(a1);// print the length of po1
		System.out.println(a2);
		assertEquals(a1,8.94427190999916,Ex2.EPS);// compare the expected length of po1 with the calculated one
		assertEquals(a2,7.002347386697557,Ex2.EPS);

	}
	@Test
	public void testLength() {
		// Test 1
		int numOfSegment = 1000;// define the number of segments
		double x1 = 0, x2 = 15;// define the range
		double[] arr1 = {0,-8.2}; // define the array of coefficients
		double res1 = Ex2.length(arr1,x1,x2,numOfSegment); // calculate the length of the curve with the given array of coefficients
		double expect1 = 123.91125856837955;// expected length for the given array of coefficients and range
		
	
		assertEquals(expect1, res1, Ex2.EPS);
	}
}