use num::bigint::BigInt;
use num::rational::{Ratio,BigRational};
use num::traits::{One,Zero};

use std::ops::{Div,Rem,Mul};
use std::cmp::{Ordering,Eq};

use super::tableau::Tableau;
use super::tableau_vars::{TableauVariable,TableauVariables};
use super::lex_min_ratio::lexminratio;

#[cfg(test)] use num::traits::{FromPrimitive,ToPrimitive};

// LCP (aka. Linear Complementarity Problem)
// =============================================================================
// (1) Mz + q >= 0
// (2) z >= 0
// (3) z'(Mz + q) = 0
//
// (1) and (2) are feasibility conditions.
// (3) is complementarity condition (also written as w = Mz + q where w and z are orthogonal)
// Lemke algorithm takes this (M, q) and a covering vector (d) and outputs a solution
//
// TODO: Does M have to be square?  We seem to make that assumption?



pub struct LCP {
	m: Vec<BigRational>,
	q: Vec<BigRational>,
	pub n: usize,
    d: Vec<BigRational>,

	pub tableau: Tableau,
    pub vars: TableauVariables,

	/* scale factors for variables z
	 * scfa[Z(0)]   for  d,  scfa[RHS] for  q
	 * scfa[Z(1..n)] for cols of  M
	 * result variables to be multiplied with these
	 */
	scale_factors: Vec<BigInt>,
}

impl LCP {

    fn new(m: Vec<BigRational>, q: Vec<BigRational>) -> LCP {

    	if m.len()%q.len() != 0 {
    		panic!("M and q are not right dimensions");  // TODO: return Result instead
    	}

    	let ncols = m.len() / q.len();
    	let nrows = q.len();
    	println!("Creating LCP with matrix M [{}x{}] and vector q [{}x1]\n", nrows, ncols, nrows);

    	if ncols != nrows {
    		panic!("M must be a square matrix but was {}x{}", nrows, ncols);
    	}

        // TODO: most of this state is just for the lemke algorithm...
    	let mut lcp = LCP {
            m: m,
            q: q,
            d: vec![BigRational::zero(); nrows],
            n: nrows,
            vars: TableauVariables::new(nrows),
            tableau: Tableau::new(nrows),
        	scale_factors: vec![BigInt::zero(); nrows+2],
        };
    	lcp.init_tableau();

    	lcp
    }

    fn init_tableau(&mut self) {

    	for j in 1..self.tableau.ncols {

            let scale_factor = {
        		self.compute_scale_factor(self.n, |i: usize| {
                    if j == self.n+1 {
                        &self.q[i]
                    } else {
                        &self.m[i*self.n+(j-1)]
                    }
        		})
            };

    		for i in 0..self.tableau.nrows {
                let value = {
        			let rat = if j == self.n+1 {
        				&self.q[i]
        			} else {
        			    &self.m[i*self.n+(j-1)]
                    };

        			/* cols 0..n of  A  contain LHS cobasic cols of  Ax = b     */
        			/* where the system is here         -Iw + dz_0 + Mz = -q    */
        			/* cols of  q  will be negated after first min ratio test   */
        			/* A[i][j] = num * (scfa[j] / den),  fraction is integral       */
                    rat.numer().mul(&scale_factor).div(rat.denom())
                };
    			self.tableau.set(i, j, value);
    		}
            self.scale_factors[j] = scale_factor;
    	}
    }

    /*
     * compute lcm  of denominators for  col  j  of  A
     * Necessary for converting fractions to integers and back again
     */
    fn compute_scale_factor<'a, F>(&'a self, n: usize, vec: F) -> BigInt where F : Fn(usize) -> &'a BigRational {

    	let mut lcm = BigInt::one();
    	for i in 0..n {
            let rat = vec(i);
            let gcd = LCP::euclid_gcd(&lcm, rat.denom());
    		lcm = lcm.div(&gcd).mul(rat.denom());
    	}
    	lcm
    }

    fn euclid_gcd(a: &BigInt, b: &BigInt) -> BigInt {
        let mut a = a.clone();
        let mut b = b.clone();
        while BigInt::zero().eq(&b) {
            let t = b;
            b = a.rem(&t);
            a = t;
        }
        a
    }

    // TODO: convert this to a builder pattern?
    pub fn add_covering_vector(&mut self, d: Vec<BigRational>) {

        self.d = d;

    	let scale_factor = self.compute_scale_factor(self.n, |i| &self.d[i]);

    	for i in 0..self.tableau.nrows {
    		let rat = &self.d[i];
    		let value = rat.numer().mul(&scale_factor).div(rat.denom());
    		self.tableau.set(i, 0, value)
    	}

        self.scale_factors[0] = scale_factor;
    }
}

/*
 * asserts that  d >= 0  and not  q >= 0  (o/w trivial sol)
 * and that q[i] < 0  implies  d[i] > 0
 */
 // TODO: don't error on trivial solution... just return it
fn validate_inputs(q: &Vec<BigRational>, d: &Vec<BigRational>) {

	let mut is_q_pos = true;
	for i in 0..q.len() {
		if d[i].cmp(&Ratio::zero()) == Ordering::Less {
			panic!("Covering vector  d[{}] = {} negative. Cannot start Lemke.", i+1, d[i]);
		} else if q[i].cmp(&Ratio::zero()) == Ordering::Less {
			is_q_pos = false;
			if d[i].is_zero() {
				panic!("Covering vector  d[{}] = 0  where  q[{}] = {}  is negative. Cannot start Lemke.", i+1, i+1, q[i]);
			}
		}
	}

	if is_q_pos {
		panic!("No need to start Lemke since  q>=0. Trivial solution  z=0.");
	}
}

fn lemke(m: Vec<BigRational>, q: Vec<BigRational>, d: Vec<BigRational>) -> Vec<BigRational> {
	lemke_with_pivot_max(m, q, d, 0)
}

// LemkeWithPivotMax solves the linear complementarity probelm via Lemke's algorithm.
// It will only perform up to maxCount pivots before exiting.
// TODO: ray termination...
fn lemke_with_pivot_max(m: Vec<BigRational>, q: Vec<BigRational>, d: Vec<BigRational>, pivot_max: usize) -> Vec<BigRational> {

	validate_inputs(&q, &d);

	// TODO: better way?
    let mut lcp = LCP::new(m, q);
	lcp.add_covering_vector(d);

	let mut enter = lcp.vars.z(0); // z0 enters the basis to obtain lex-feasible solution
	let (mut leave, mut z0_can_leave) = lexminratio(&lcp.tableau, &lcp.vars, &enter);

	lcp.vars.negate_rhs(&mut lcp.tableau); // now give the entering q-col its correct sign

	let mut pivot_count = 1;
	loop {

		println!("{}. entering: {}, leaving: {}", pivot_count, enter.to_string(), leave.to_string());

		lcp.vars.pivot(&mut lcp.tableau, &leave, &enter);

		if z0_can_leave {
			break; // z0 will have a value of zero but may still be basic... amend?
		}

		enter = leave.complement();  // select pivot

		// FIXME: better way?  Maybe rust will fix this...
		let (next_leave, next_z0_can_leave) = lexminratio(&lcp.tableau, &lcp.vars, &enter);
        leave = next_leave;
        z0_can_leave = next_z0_can_leave;

		if pivot_count == pivot_max {
			// maxcount == 0 is equivalent to infinity since pivotcount starts at 1
			// TODO: negative result...
			break;
		}

		pivot_count += 1;
	}

	lcp.vars.solution(&lcp.tableau, &lcp.scale_factors)
}

#[cfg(test)]
fn into_bigrats(ints: Vec<i64>) -> Vec<BigRational> {
    ints.into_iter()
        .map(|val| { Ratio::from_integer(BigInt::from_i64(val).unwrap()) })
        .collect()
}

#[test]
fn lemke2() {

	let m = into_bigrats(vec![2, 1, 1, 3]);
	let q = into_bigrats(vec![-1, -1]);
	let d = into_bigrats(vec![2, 1]);

	let z = lemke(m, q, d);

	assert_eq!(2, z.len());
	assert_eq!(false, z[0].is_integer());
	assert_eq!(2, z[0].numer().to_i64().unwrap());
	assert_eq!(5, z[0].denom().to_i64().unwrap());
	assert_eq!(false, z[1].is_integer());
	assert_eq!(1, z[1].numer().to_i64().unwrap());
	assert_eq!(5, z[1].denom().to_i64().unwrap());
}

#[test]
fn lemke3() {

	let m = into_bigrats(vec![0, -1, 2, 2, 0, -2, -1, 1, 0]);
	let q = into_bigrats(vec![-3, 6, -1]);
	let d = into_bigrats(vec![1, 1, 1]);

	let z = lemke(m, q, d);

	assert_eq!(3, z.len());
	assert_eq!(true, z[0].is_integer());
	assert_eq!(0, z[0].numer().to_i64().unwrap());
	assert_eq!(true, z[1].is_integer());
	assert_eq!(1, z[1].numer().to_i64().unwrap());
	assert_eq!(true, z[2].is_integer());
	assert_eq!(3, z[2].numer().to_i64().unwrap());
}
