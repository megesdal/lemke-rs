use num::bigint::BigInt;
use num::rational::BigRational;
use num::traits::{One,Zero};

use std::ops::{Div,Rem,Mul};

use tableau::Tableau;

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

struct TableauVariable {
    value: usize,
    n: usize,
}

impl TableauVariable {

    fn is_z(&self) -> bool {
    	self.value <= self.n
    }

    fn is_w(&self) -> bool {
    	self.value > self.n
    }

    /*
     * complement of  v  in VARS, error if  v==Z(0).
     * this is  W(i) for Z(i)  and vice versa, i=1...n
     */
    fn complement(&self) -> TableauVariable {

    	if self.value == 0 {
    		panic!("Attempt to find complement of z0.");
    	}

    	if self.is_z() {
            TableauVariable{ value: self.value + self.n, n: self.n }
    	} else {
            TableauVariable{ value: self.value - self.n, n: self.n }
        }
    }
}

struct TableauVariables {
    /*  v in VARS, v cobasic:  TABCOL(v) is v's tableau col */
	/*  v  basic:  TABCOL(v) < 0,  TABCOL(v)+n   is v's row */
	/* VARS   = 0..2n = Z(0) .. Z(n) W(1) .. W(n)           */
	/* ROWCOL = 0..2n,  0 .. n-1: tabl rows (basic vars)    */
	/*                  n .. 2n:  tabl cols  0..n (cobasic) */
	vars2rowcol: Vec<usize>,
	rowcol2vars: Vec<usize>,
    n: usize,
}

impl TableauVariables {
    /*
     * init tableau variables:
     * Z(0)...Z(n)  nonbasic,  W(1)...W(n) basic
     * This is for setting up a complementary basis/cobasis
     */
    fn new(n: usize) -> TableauVariables {

        let mut vars = TableauVariables {
            n: n,
            vars2rowcol: vec![0; 2*n+1],
            rowcol2vars: vec![0; 2*n+1],
        };

        // w[1..n]
        for i in 1..n+1 {
            vars.rowcol2vars[i-1] = i + n;
            vars.vars2rowcol[i+n] = i - 1;
        }

        // z[0..n]
        for i in 0..n+1 {
    		vars.rowcol2vars[i+n] = i;
            vars.vars2rowcol[i] = n + i;
    	}

        vars
    }

    // TODO: check bounds with panic?
    fn z(&self, idx: usize) -> TableauVariable {
        TableauVariable{ value: idx, n: self.n }
    }

    fn w(&self, idx: usize) -> TableauVariable {
        TableauVariable{ value: idx + self.n, n: self.n }
    }

    fn from_row(&self, row: usize) -> TableauVariable {
        TableauVariable{ value: self.rowcol2vars[row], n: self.n }
    }

    fn from_col(&self, col: usize) -> TableauVariable {
    	TableauVariable{ value: self.rowcol2vars[col+self.n], n: self.n }
    }

    fn to_row(&self, var: &TableauVariable) -> usize {
        self.vars2rowcol[var.value]
    }

    fn to_col(&self, var: &TableauVariable) -> usize {
        self.vars2rowcol[var.value] - self.n
    }

    fn is_basic(&self, var: &TableauVariable) -> bool {
    	self.vars2rowcol[var.value] < self.n
    }

    fn swap(&mut self, enter: &TableauVariable, leave: &TableauVariable) -> (usize, usize) {

    	let leaveRow = self.to_row(leave); // basic var is leaving
    	let enterCol = self.to_col(enter); // cobasic var is entering

    	self.vars2rowcol[leave.value] = enterCol + self.n;
    	self.rowcol2vars[enterCol+self.n] = leave.value;

    	self.vars2rowcol[enter.value] = leaveRow;
    	self.rowcol2vars[leaveRow] = enter.value;

    	(leaveRow, enterCol)
    }
}

struct LCP {
	m: Vec<BigRational>,
	q: Vec<BigRational>,
	pub n: usize,

	a: Tableau,
    vars: TableauVariables,

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

    	let mut lcp = LCP {
            m: m,
            q: q,
            n: nrows,
            vars: TableauVariables::new(nrows),
            a: Tableau::new(nrows, nrows+2),
        	scale_factors: vec![BigInt::zero(); nrows+2],
        };
    	lcp.init_tableau();

    	lcp
    }

    fn init_tableau(&mut self) {

    	for j in 1..self.a.ncols {

            let scale_factor = {
        		self.compute_scale_factor(self.n, |i: usize| {
                    if j == self.n+1 {
                        &self.q[i]
                    } else {
                        &self.m[i*self.n+(j-1)]
                    }
        		})
            };

    		for i in 0..self.a.nrows {
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
    			self.a.set(i, j, value);
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
}

#[test]
fn variable_assignments_work() {

	let vars = TableauVariables::new(4);

	for i in 0..vars.n+1 {
		let zi = vars.z(i);
		assert_eq!(i, zi.value);  // z(i) should equal i
		assert_eq!(i, vars.to_col(&zi));  // col(z(i)) should equal i
		assert_eq!(zi.value, vars.from_col(i).value);  // var(col(z(i))) = var(i) should equal z(i)
		assert_eq!(false, vars.is_basic(&zi));  // z(i) should NOT be basic
	}

	for i in 1..vars.n+1 {
		let wi = vars.w(i);
		assert_eq!(i+vars.n, wi.value);  // w(i) should equal i + n
		assert_eq!(i-1, vars.to_row(&wi));  // row(w(i)) should equal i - 1 (w is 1-indexed)
		assert_eq!(wi.value, vars.from_row(i-1).value);  // var(row(w(i - 1))) = var(i - 1) should equal w(i)
		assert_eq!(true, vars.is_basic(&wi));  // w(i) should be basic
	}
}

#[test]
fn complement_works() {

    let vars = TableauVariables::new(4);

	for i in 1..vars.n+1 {
		let compVar = vars.z(i).complement();
		assert_eq!(vars.w(i).value, compVar.value);  // w(i) should be complement of z(i)
	}

	//_, err := lcp.complement(lcp.z(0))
	//assert.NotNil(t, err, "Should not be able to get complement of z0")
}

#[test]
fn swap_works() {

	let mut vars = TableauVariables::new(4);

	let leave_var = vars.w(1);                    // first row, w(1) = 4
	assert_eq!(true, vars.is_basic(&leave_var));    // leaving var should be basic
	assert_eq!(0, vars.to_row(&leave_var));        // w(1) points to row 0

	let enter_var = vars.z(0);                    // first col, z(0) = 0
	assert_eq!(false, vars.is_basic(&enter_var));  // entering var should be cobasic
	assert_eq!(0, vars.to_col(&enter_var));        // z(0) points to col 0

	let (row, col) = vars.swap(&enter_var, &leave_var);

	assert_eq!(0, col);  // col should be original col(w1)
	assert_eq!(0, row);  // row should be original row(z0)

	assert_eq!(leave_var.value, vars.from_col(col).value);  // w1 should be var(col)
	assert_eq!(enter_var.value, vars.from_row(row).value);  // z0 should be var(row)
}
