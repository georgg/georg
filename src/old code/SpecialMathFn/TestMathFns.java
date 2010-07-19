package edu.mit.csail.psrg.georg.SpecialMathFn;

public class TestMathFns {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
	//	testDigamma();
		testTrigamma();
	}
	
	public static void testDigamma() {
		double d = 0.0;
		d = Digamma.eval(9) - 2.140641477955609996536345;
		System.out.println(d);
		d = Digamma.eval(2.5) - 0.7031566406452431872257;
		System.out.println(d);
		d = Digamma.eval(0.1) - (-10.42375494041107679516822);
		System.out.println(d);
		d = Digamma.eval(7e-4) - (-1429.147493371120205005198);
		System.out.println(d);
		d = Digamma.eval(7e-5) - (-14286.29138623969227538398);
		System.out.println(d);
		d = Digamma.eval(7e-6) - (-142857.7200612932791081972);
		System.out.println(d);
		d = Digamma.eval(2e-6) - (-500000.5772123750382073831);
		System.out.println(d);
		d = Digamma.eval(1e-6) - (-1000000.577214019968668068);
		System.out.println(d);
		d = Digamma.eval(7e-7) - (-1428572.005785942019703646);
		System.out.println(d);
		d = Digamma.eval(-0.5) - (.03648997397857652055902367);
		System.out.println(d);
		d = Digamma.eval(-1.1) - 10.15416395914385769902271;
		System.out.println(d);
		d = Digamma.eval(-1);
		System.out.println(d);
		d = Digamma.eval(0);
		System.out.println(d);
	}
	
	public static void testTrigamma() {
		double d = 0.0;
		d = Trigamma.eval(9)    - .117512014694031425134547;
		System.out.println(d);
		d = Trigamma.eval(2.5)  - .4903577561002348649728011;
		System.out.println(d);
		d = Trigamma.eval(0.1)  - 101.4332991507927588172155;
		System.out.println(d);
		d = Trigamma.eval(7e-4) - 2040817.969783389022409943;
		System.out.println(d);
		d = Trigamma.eval(7e-5) - 204081634.2978270192803090;
		System.out.println(d);
		d = Trigamma.eval(7e-6) - 20408163266.95103968719027;
		System.out.println(d);
		d = Trigamma.eval(7e-7) - 2040816326532.257177281929;
		System.out.println(d);
		d = Trigamma.eval(-0.5) - 8.934802200544679309417246;
		System.out.println(d);
		d = Trigamma.eval(-1.1) - 102.7489862404689390536693;
		System.out.println(d);
		d = Trigamma.eval(-1);
		System.out.println(d);
		d = Trigamma.eval(0);
		System.out.println(d);
	}

}
